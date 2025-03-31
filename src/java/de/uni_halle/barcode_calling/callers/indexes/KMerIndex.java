package de.uni_halle.barcode_calling.callers.indexes;

import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.Sorting;
import de.uni_halle.barcode_calling.util.KMer;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;

/**
 * Implements the heuristic index from Section 3.4.4
 * @author Riko Uphoff
 *
 */
public class KMerIndex implements OrderedIndex {
	
	protected static final double DEFAULT_BARCODES_FRACTION_PROBED = 0.01;
	protected static final int MAX_POSITION_DIFFERENCE = 5;
	
	protected int[][][] index;
	protected int barcodeLength, barcodesCount;
	protected int k;
	protected double barcodesFractionProbed;
	protected BarcodeDataset barcodes;
	
	public KMerIndex(BarcodeDataset barcodes, int k, double barcodesFractionProbed) {
		this.setK(k);
		this.setBarcodes(barcodes);
		this.barcodesFractionProbed = barcodesFractionProbed;
	}
	
	public KMerIndex(BarcodeDataset barcodes, int k) {
		this(barcodes, k, DEFAULT_BARCODES_FRACTION_PROBED);
	}
	
	public KMerIndex(int k) {
		this(new BarcodeDataset(), k, DEFAULT_BARCODES_FRACTION_PROBED);
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public void setBarcodes(BarcodeDataset barcodes) {
		this.barcodeLength = barcodes.getSequenceLength();
		this.barcodesCount = barcodes.getSize();
		this.barcodes = barcodes;
		
		if (this.barcodesCount > 0)
			this.insert(barcodes);
	}
	
	@Override
	public BarcodeDataset getBarcodes() {
		return this.barcodes;
	}
	
	public int getK() {
		return this.k;
	}
	
	public void setK(int k) {
		this.k = k;
	}
	
	@Override
	public String getName() {
		StringBuilder sb = new StringBuilder(this.getClass().getSimpleName());
		
		sb.append("_k");
		sb.append(this.k);
		
		return sb.toString();
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(this.getName());
		
		sb.append("_pos");
		sb.append(this.barcodesFractionProbed);
		
		return sb.toString();
	}
	
	@Override
	public double getBarcodesFractionProbed() {
		return this.barcodesFractionProbed;
	}
	
	@Override
	public void setBarcodesFractionProbed(double barcodesFractionProbed) {
		this.barcodesFractionProbed = barcodesFractionProbed;
	}
	
	public void insert(DNASequence barcode, int id, List<Integer>[][] helpIndex) {
		int pos = 0;
		
		for (KMer kMer : KMer.iterateKMers(barcode, this.k)) {
			pos = kMer.getPosition();
			int kMerHash = kMer.hashCode();
			
			if (helpIndex[kMerHash][pos] == null) {
				helpIndex[kMerHash][pos] = new ArrayList<Integer>();
			}
			
			helpIndex[kMerHash][pos].add(id);
			
		}
	}
	
	public void insert(BarcodeDataset barcodes) {
		@SuppressWarnings("unchecked")
		List<Integer>[][] helpIndex = new List[(int)Math.pow(4, k)][this.barcodeLength - k + 1];
		
		for (int i = 0; i < barcodes.getSize(); i++) {
			this.insert(barcodes.getSequence(i), i, helpIndex);
		}
		
		// Get max list length
		int maxSize = 0;
		for (int kMerHash = 0; kMerHash < (int)Math.pow(4, k); kMerHash++) {
			for (int pos = 0; pos < this.barcodeLength - k + 1; pos++) {
				if (helpIndex[kMerHash][pos] != null) {
					maxSize = maxSize < helpIndex[kMerHash][pos].size() ? helpIndex[kMerHash][pos].size() : maxSize;
				}
			}
		}
		
		// Create static index
		this.index = new int[(int)Math.pow(4, k)][this.barcodeLength - k + 1][maxSize];
		for (int kMerHash = 0; kMerHash < (int)Math.pow(4, k); kMerHash++) {
			for (int pos = 0; pos < this.barcodeLength - k + 1; pos++) {
				if (helpIndex[kMerHash][pos] == null) {
					this.index[kMerHash][pos][0] = -1;
				}
				else {
					for (int i = 0; i < helpIndex[kMerHash][pos].size(); i++) {
						this.index[kMerHash][pos][i] = helpIndex[kMerHash][pos].get(i);
					}
					
					if (helpIndex[kMerHash][pos].size() < maxSize)
						this.index[kMerHash][pos][helpIndex[kMerHash][pos].size()] = -1;
				}
			}
		}
	}
	
	public int insertHitDistances(
			DNASequence seq,
			int[] candidatesScore,
			int[] candidatesId
	) {
		int pos, lookupPos, kMerHash;
		int candidatesCount = 0;
		
		for (KMer kMer : KMer.iterateKMers(seq, this.k)) {
			
			pos = kMer.getPosition();
			kMerHash = kMer.hashCode();
			
			for (int offset = 0; offset <= MAX_POSITION_DIFFERENCE; offset++) {
				for (int sign = -1; sign <= 1 && offset != 0 || sign == -1; sign += 2) {
					lookupPos = pos + sign*offset;
					
					if (lookupPos >= 0 && lookupPos < this.barcodeLength - k + 1) {
						int[] barcodeEntries = index[kMerHash][lookupPos];
						int i = 0;
						
						while (i < barcodeEntries.length && barcodeEntries[i] != -1) {
							// Get k-mer distance
							// Count number of kMer matches per barcode
							int candidate = barcodeEntries[i];
							
							if (candidatesScore[candidate] == 0) {
								candidatesId[candidatesCount] = candidate;
								candidatesCount++;
							}
							
							// Increment pseudo distance
							candidatesScore[candidate] += offset - this.barcodeLength;
							
							i++;
						}
					}
				}
			}
		}
		
		return candidatesCount;
	}
	
	public int[] getOrderedCandidates(DNASequence seq) {
		int[] candidatesScore = new int[this.barcodesCount];
		int[] candidatesId = new int[this.barcodesCount];
		
		// Score
		int candidatesCount = this.insertHitDistances(seq, candidatesScore, candidatesId);
		
		// Sort
		Sorting.mergeSort(candidatesId, candidatesCount, (int a, int b) -> candidatesScore[a] - candidatesScore[b]);
		
		int indexSize = (int)Math.ceil(this.barcodesCount * this.barcodesFractionProbed);
		indexSize = Math.min(indexSize, candidatesCount);
		int[] candidates = new int[indexSize];
		
		// Shorten
		for (int i = 0; i < indexSize; i++)
			candidates[i] = candidatesId[i];
		
		return candidates;
	}

}
