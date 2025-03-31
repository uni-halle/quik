package de.uni_halle.barcode_calling.callers.indexes;

import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.Sorting;
import de.uni_halle.barcode_calling.util.KMer;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;

/**
 * Implementation of KMerIndex that counts the number of inner iterations. Only to be used in experiments.
 * @author Riko Uphoff
 *
 */
public class KMerIndexCounter extends KMerIndex {
	// TODO Store all runs, not just average
	
	private static final double DEFAULT_BARCODES_FRACTION_PROBED = 0.01;
	private static final int MAX_POSITION_DIFFERENCE = 5;
	
	private int readCounter = 0;
	public List<Integer> operationCounter = new ArrayList<Integer>();
	public List<Integer> candidateCounter = new ArrayList<Integer>();
	public List<Double> timeInsertHitDistance = new ArrayList<Double>();
	public List<Double> timeInsertHitDistancePos = new ArrayList<Double>();
	public List<Double> timeInsertHitDistanceScore = new ArrayList<Double>();
	public List<Double> timeSorting = new ArrayList<Double>();
	public List<Double> timeTotal = new ArrayList<Double>();
	
	public KMerIndexCounter(BarcodeDataset barcodes, int k, double barcodesFractionProbed) {
		super(barcodes, k, barcodesFractionProbed);
	}
	
	public KMerIndexCounter(BarcodeDataset barcodes, int k) {
		this(barcodes, k, DEFAULT_BARCODES_FRACTION_PROBED);
	}
	
	public KMerIndexCounter(int k) {
		this(new BarcodeDataset(), k, DEFAULT_BARCODES_FRACTION_PROBED);
	}
	
	public void setBarcodes(BarcodeDataset barcodes) {
		super.setBarcodes(barcodes);
		this.resetCounter();
	}
	
	public void resetCounter() {
		this.readCounter = 0;
		this.operationCounter = new ArrayList<Integer>();
		this.candidateCounter = new ArrayList<Integer>();
		this.timeInsertHitDistance = new ArrayList<Double>();
		this.timeInsertHitDistancePos = new ArrayList<Double>();
		this.timeInsertHitDistanceScore = new ArrayList<Double>();
		this.timeSorting = new ArrayList<Double>();
		this.timeTotal = new ArrayList<Double>();
	}
	
	public synchronized void incrementReadCounter() {
		this.readCounter++;
	}
	
	public synchronized void addCandidateCounter(int candidateCount) {
		this.candidateCounter.add(candidateCount);
	}
	
	public synchronized void addOperationCounter(int operationCounter) {
		this.operationCounter.add(operationCounter);
	}
	
	public synchronized void addTimeInsertHitDistance(double time) {
		this.timeInsertHitDistance.add(time);
	}
	
	public synchronized void addTimeInsertHitDistancePos(double time) {
		this.timeInsertHitDistancePos.add(time);
	}
	
	public synchronized void addTimeInsertHitDistanceScore(double time) {
		this.timeInsertHitDistanceScore.add(time);
	}
	
	public synchronized void addTimeSorting(double time) {
		this.timeSorting.add(time);
	}
	
	public synchronized void addTimeTotal(double time) {
		this.timeTotal.add(time);
	}
	
	public synchronized List<Integer> getCandidateCounter() {
		return this.candidateCounter;
	}
	
	public synchronized List<Integer> getOperationCounter() {
		return this.operationCounter;
	}
	
	public double getAverageOperationCounter() {
		return this.operationCounter.stream().mapToLong(Integer::longValue).sum() / (double)this.readCounter;
	}
	
	public double getAverageCandidateCounter() {
		return this.candidateCounter.stream().mapToLong(Integer::longValue).sum() / (double)this.readCounter;
	}
	
	@Override
	public String getName() {
		StringBuilder sb = new StringBuilder(KMerIndex.class.getSimpleName());
		
		sb.append("_k");
		sb.append(this.k);
		
		return sb.toString();
	}
	
	public int[] insertHitDistancesCounter(
			DNASequence seq,
			int[] candidatesScore,
			int[] candidatesId
	) {
		double timeTotal = 0, timePos = 0, timeCandidatesScore = 0;
		long startTime, finishTime, startTimeTotal, finishTimeTotal;
		
		startTimeTotal = System.nanoTime();
		
		// Start
		int pos, lookupPos, kMerHash;
		int candidatesCount = 0;
		int operationCounter = 0;
		
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
							operationCounter++;
							
							i++;
						}
					}
				}
			}
		}
		
		finishTimeTotal = System.nanoTime();
		timeTotal = (double)(finishTimeTotal - startTimeTotal)/1000000;  // in ms
		
		
		this.addTimeInsertHitDistance(timeTotal);
		//this.addTimeInsertHitDistancePos(timePos);
		//this.addTimeInsertHitDistanceScore(timeCandidatesScore);
		
		return new int[] {operationCounter, candidatesCount};
	}
	
	@Override
	public int[] getOrderedCandidates(DNASequence read) {
		long startTimeTotal = System.nanoTime();
		
		int[] candidatesScore = new int[this.barcodesCount];
		int[] candidatesId = new int[this.barcodesCount];
		
		// Score
		int[] res = this.insertHitDistancesCounter(read, candidatesScore, candidatesId);
		int hitCounter = res[0];
		int candidatesCount = res[1];
		
		// Sort
		long startTime = System.nanoTime();
		
		// TODO Different sorting algorithm
		Sorting.mergeSort(candidatesId, candidatesCount, (int a, int b) -> candidatesScore[a] - candidatesScore[b]);
		
		long finishTime = System.nanoTime();
		double timeSorting = (double)(finishTime - startTime)/1000000;  // in ms
		this.addTimeSorting(timeSorting);
		
		// Shorten
		int indexSize = (int)Math.ceil(this.barcodesCount * this.barcodesFractionProbed);
		indexSize = Math.min(indexSize, candidatesCount);
		int[] candidates = new int[indexSize];
		
		// Shorten
		for (int i = 0; i < indexSize; i++)
			candidates[i] = candidatesId[i];
		
		long finishTimeTotal = System.nanoTime();
		double timeTotal = (double)(finishTimeTotal - startTimeTotal)/1000000;  // in ms
		this.addTimeTotal(timeTotal);
		
		// Count
		this.incrementReadCounter();
		this.addOperationCounter(hitCounter);
		this.addCandidateCounter(indexSize);
		
		return candidates;
	}
}