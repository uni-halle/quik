package de.uni_halle.barcode_calling.callers;

import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;

/**
 * Base for calling algorithms that stop once they find a barcode within distance maxDist.
 * @author Riko Uphoff
 *
 */
public class MaxDistanceCaller implements BarcodeCaller{
	
	protected DistanceMeasureBase distanceMeasure;
	int maxDist;
	private BarcodeDataset barcodes;
	long executionTime;
	
	public MaxDistanceCaller(BarcodeDataset barcodes, DistanceMeasureBase distanceMeasure, int maxDist) {
		this.distanceMeasure = distanceMeasure;
		this.barcodes = barcodes;
		this.maxDist = maxDist;
	}
	
	public MaxDistanceCaller(DistanceMeasureBase distanceMeasure, int maxDist) {
		this(new BarcodeDataset(), distanceMeasure, maxDist);
	}
	
	public MaxDistanceCaller(BarcodeDataset barcodes, DistanceMeasureBase distanceMeasure) {
		this(
				barcodes, 
				distanceMeasure, 
				distanceMeasure.getDefaultMaxDistance(barcodes.getSequenceLength())
		);
	}
	
	@Override
	public BarcodeDataset getBarcodes() {
		return this.barcodes;
	}
	
	public void setBarcodes(BarcodeDataset barcodes, int maxDist) {
		this.barcodes = barcodes;
		this.setMaxDist(maxDist);
	}
	
	public void setBarcodes(BarcodeDataset barcodes) {
		this.barcodes = barcodes;
	}
	
	/**
	 * Sets the maximum distance to which barcodes are accepted.
	 * @param maxDist Maximum accepted distance.
	 */
	public void setMaxDist(int maxDist) {
		this.maxDist = Math.max(0, maxDist);
	}
	
	public int getMaxDist() {
		return this.maxDist;
	}

	public synchronized void resetExecutionTimes() {
		this.executionTime = 0;
	}
	
	public synchronized void addExecTime(long t) {
		this.executionTime += t;
	}
	
	public synchronized void setExecTime(long t) {
		this.executionTime = t;
	}
	
	public long getExecTime() {
		return this.executionTime;
	}
	
	public int[] callDetails(DNASequence read) {
		int callId = -1, callIdSecond = -1;
		int minDist = this.maxDist + 1, minDistSecond = this.maxDist + 1;
		int dist;
		DNASequence barcode;
		
		long startTime = System.nanoTime();	// Track performance
		for (int barcodeId = 0; barcodeId < this.barcodes.getSize(); barcodeId++) {
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
		long finishTime = System.nanoTime();	// Track performance
		this.addExecTime(finishTime - startTime);	// Track performance
		
		int[] out = {callId, minDist, callIdSecond, minDistSecond};
		return out;
	}

	@Override
	public DNASequence getBarcode(int id) {
		return this.barcodes.getSequences()[id];
	}
	
	@Override
	public String getName() {
		return this.getClass().getSimpleName();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder(this.getName());
		
		//sb.append("_tau");
		//sb.append(this.getBarcodes().getMinMutualEditDist() / 2);
		sb.append("_L");
		sb.append(maxDist);
		
		return sb.toString();
	}
	
}
