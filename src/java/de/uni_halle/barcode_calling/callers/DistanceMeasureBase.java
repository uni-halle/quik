package de.uni_halle.barcode_calling.callers;

import de.uni_halle.barcode_calling.util.DNASequence;

public interface DistanceMeasureBase {
	
	public int distance(DNASequence a, DNASequence b);
	
	/**
	 * @param a Sequence A
	 * @param b Sequence B
	 * @param maxDist Maximum computable distance
	 * @return Are A and B within a distance of maxDist?
	 */
	public boolean inDistance(DNASequence a, DNASequence b, int maxDist);
	
	/**
	 * @param a Sequence A
	 * @param b Sequence B
	 * @param maxDist Maximum computable distance
	 * @return Distance between a and b up to maxDist
	 */
	public int inDistanceResult(DNASequence a, DNASequence b, int maxDist);
	
	/**
	 * Default maximum accepted distance.
	 * @param length Length of barcodes.
	 * @return Default L.
	 */
	public int getDefaultMaxDistance(int length);
	
	public String getName();
	public String toString();
}
