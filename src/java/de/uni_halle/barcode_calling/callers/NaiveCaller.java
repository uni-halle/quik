package de.uni_halle.barcode_calling.callers;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

/**
 * Provides an abstract base for methods implementing the first baseline from the paper with SL distance (see Section 3.4.1)
 * @author Riko Uphoff
 *
 */
public abstract class NaiveCaller implements BarcodeCaller {
	
	private DistanceMeasureBase distanceMeasure;
	private BarcodeDataset barcodes;
	private long executionTime;
	
	public NaiveCaller(BarcodeDataset barcodes, DistanceMeasureBase distanceMeasure) {
		this.barcodes = barcodes;
		this.distanceMeasure = distanceMeasure;
	}
	
	public NaiveCaller(DistanceMeasureBase distanceMeasure) {
		this(new BarcodeDataset(), distanceMeasure);
	}
	
	@Override
	public BarcodeDataset getBarcodes() {
		return this.barcodes;
	}
	
	public void setBarcodes(BarcodeDataset barcodes) {
		this.barcodes = barcodes;
	}
	
	@Override
	public double[] tuneParameters(LabeledDataset readsArtif, double minPrec, double minRec) {
		Main.writeToLog("No tuning needed for caller " + this.getName());
		return new double[3];
	}
	
	@Override
	public synchronized void resetExecutionTimes() {
		this.executionTime = 0;
	}
	
	@Override
	public synchronized void addExecTime(long t) {
		this.executionTime += t;
		//System.out.println("New exec time NN: " + this.executionTimeNN);	// TODO
	}
	
	@Override
	public synchronized void setExecTime(long t) {
		this.executionTime = t;
	}
	
	@Override
	public long getExecTime() {
		return this.executionTime;
	}
	
	@Override
	public int[] callDetails(DNASequence read) {
		long startTime = System.nanoTime();	// Track performance
		
		int callId = -1, callIdSecond = -1;
		int callDist = Integer.MAX_VALUE, callDistSecond = Integer.MAX_VALUE;
		int dist;
		
		for (int i = 0; i < this.barcodes.getSize(); i++) {
			dist = distanceMeasure.distance(read, this.barcodes.getSequence(i));
			if (dist < callDist) {
				// Remember second best call
				if (callDist < callDistSecond) {
					callIdSecond = callId;
					callDistSecond = callDist;
				}
				
				callId = i;
				callDist = dist;
			}
			// Remember second best call
			else if (dist < callDistSecond) {
				callIdSecond = i;
				callDistSecond = dist;
			}
		}
		
		long finishTime = System.nanoTime();	// Track performance
		this.addExecTime(finishTime - startTime);	// Track performance
		
		return new int[] {callId, callDist, callIdSecond, callDistSecond};
	}

	@Override
	public DNASequence getBarcode(int id) {
		return this.barcodes.getSequences()[id];
	}
	
	@Override
	public String toString() {
		return this.getName();
	}
	
}
