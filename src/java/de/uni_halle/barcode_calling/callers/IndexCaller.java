package de.uni_halle.barcode_calling.callers;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.indexes.OrderedIndex;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

/**
 * Implements the idea of an index-based calling algorithm
 * @author Riko Uphoff
 *
 */
public class IndexCaller extends MaxDistanceCaller {
	
	public static DistanceMeasureBase DEFAULT_DISTANCE_MEASURE = new SequenceLevenshteinBase();
	
	protected long executionTimeIndex = 0;
	protected long executionTimeNN = 0;
	protected long executionTimeTotal = 0;
	protected BarcodeDataset barcodes;
	protected OrderedIndex index;

	/*
	 * Constructors
	 */
	public IndexCaller(BarcodeDataset barcodes, OrderedIndex index, int maxDist, double indexFractionProbed) {
		super(barcodes, DEFAULT_DISTANCE_MEASURE, maxDist);
		
		this.index = index;
		this.index.setBarcodesFractionProbed(indexFractionProbed);
		this.barcodes = barcodes;
		this.index.setBarcodes(barcodes);
		this.maxDist = maxDist;
	}
	
	public IndexCaller(BarcodeDataset barcodes, OrderedIndex index, int maxDist) {
		this(barcodes, index, maxDist, 1);
	}
	
	public IndexCaller(BarcodeDataset barcodes, OrderedIndex index) {
		this(
				barcodes, 
				index, 
				DEFAULT_DISTANCE_MEASURE.getDefaultMaxDistance(barcodes.getSequenceLength()), 
				1.0
		);
	}
	
	public IndexCaller(OrderedIndex index) {
		this(new BarcodeDataset(), index);
	}
	
	@Override
	public String getName() {
		StringBuilder sb = new StringBuilder(this.getClass().getSimpleName());
		
		sb.append("_index(");
		sb.append(this.index.getName());
		sb.append(")");
		
		return sb.toString();
	}
	
	public void setIndex(OrderedIndex index) {
		this.index = index;
	}
	
	public OrderedIndex getIndex() {
		return this.index;
	}
	
	public void setDistanceModel(DistanceMeasureBase distanceMeasure) {
		if (this.barcodes != null) {
			this.setMaxDist(
					distanceMeasure.getDefaultMaxDistance(this.barcodes.getSequenceLength())
			);
		}
		else {
			this.setMaxDist(0);
		}
		
		this.distanceMeasure = distanceMeasure;
	}
	
	public void setBarcodesFractionProbed(double barcodesFractionProbed) {
		this.index.setBarcodesFractionProbed(barcodesFractionProbed);
	}
	
	public double getBarcodesFractionProbed() {
		return this.index.getBarcodesFractionProbed();
	}
	
	public synchronized void resetExecutionTimes() {
		this.executionTimeIndex = 0;
		this.executionTimeNN = 0;
		this.executionTimeTotal = 0;
	}
	
	public synchronized void addExecTimeIndex(long t) {
		this.executionTimeIndex += t;
	}
	
	public synchronized void addExecTimeNN(long t) {
		this.executionTimeNN += t;
	}
	
	public synchronized void addExecTime(long t) {
		this.executionTimeTotal += t;
	}
	
	public synchronized void setExecTimeIndex(long t) {
		this.executionTimeIndex = t;
	}
	
	public synchronized void setExecTimeNN(long t) {
		this.executionTimeNN = t;
	}
	
	public synchronized void setExecTime(long t) {
		this.executionTimeTotal = t;
	}
	
	public long getExecTimeIndex() {
		return this.executionTimeIndex;
	}
	
	public long getExecTimeNN() {
		return this.executionTimeNN;
	}
	
	public long getExecTime() {
		return this.executionTimeTotal;
	}
	

	@Override
	public void setBarcodes(BarcodeDataset barcodes) {
		this.barcodes = barcodes;
		this.index.setBarcodes(barcodes);
	}
	
	@Override
	public DNASequence getBarcode(int id) {
		return this.barcodes.getSequences()[id];
	}
	
	@Override
	public BarcodeDataset getBarcodes() {
		return this.barcodes;
	}
	
	/**
	 * Calls the read and returns the presumed barcode id and the SL distance
	 * @param read
	 * @return [callId, dist]
	 */
	public int[] callDetails(DNASequence read) {
		int callId = -1, callIdSecond = -1;
		int minDist = this.maxDist + 1, minDistSecond = this.maxDist + 1;
		int dist;
		DNASequence barcode;
		
		long startTime = System.nanoTime();	// Track performance
		int[] candidates = index.getOrderedCandidates(read);
		long finishTime = System.nanoTime();	// Track performance
		this.addExecTimeIndex(finishTime - startTime);	// Track performance
		
		startTime = System.nanoTime();	// Track performance
		for (int barcodeId : candidates) {
			barcode = this.getBarcode(barcodeId);
			
			dist = distanceMeasure.inDistanceResult(read, barcode, minDistSecond-1);
			
			if (dist < minDist) {
				if (minDist < minDistSecond) {
					minDistSecond = minDist;
					callIdSecond = callId;
				}
				minDist = dist;
				callId = barcodeId;
			}
			else if (dist < minDistSecond) {
				minDistSecond = dist;
				callIdSecond = barcodeId;
			}
		}
		finishTime = System.nanoTime();	// Track performance
		this.addExecTimeNN(finishTime - startTime);	// Track performance
		
		int[] out = {callId, minDist, callIdSecond, minDistSecond};
		return out;
	}
	
	/**
	 * Result of calling read i is written into callingResults[i] and in case of call, the distance is written into callingDistances[i]
	 * @param reads
	 * @param callingResults
	 * @param callingDistances
	 */
	public void callParallel(
			DNASequence[] reads, 
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond
		) {
		super.callParallel(
				reads,
				callingResults, 
				callingDistances,
				callingResultsSecond, 
				callingDistancesSecond
			);
		
		double indexExecTimeFrac = (double)this.getExecTimeIndex() / (double)(this.getExecTimeNN() + this.getExecTimeIndex());
		this.setExecTimeIndex((long)(this.getExecTime() * indexExecTimeFrac));
		this.setExecTimeNN((long)(this.getExecTime() * (1-indexExecTimeFrac)));
	};
	
	@Override
	public double[] tuneParameters(LabeledDataset readsArtif, double minPres, double minRec) {
		
		// Index param. pos
		double recIndex = Math.sqrt(minRec);
		this.index.tuneParameters(readsArtif, recIndex);
		
		// Caller param. maxDist (L)
		double[] out = super.tuneParameters((LabeledDataset)readsArtif, minPres, minRec);
		
		return out;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(this.getClass().getSimpleName());
		
		sb.append("_index(");
		sb.append(this.index.getName());
		sb.append(")");
		
		sb.append("_L");
		sb.append(maxDist);
		sb.append("_pos");
		sb.append(this.index.getBarcodesFractionProbed());
		
		return sb.toString();
	}
}
