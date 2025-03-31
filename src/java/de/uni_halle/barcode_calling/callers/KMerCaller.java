package de.uni_halle.barcode_calling.callers;

import java.util.Map;

import de.uni_halle.barcode_calling.callers.indexes.KMerIndex;
import de.uni_halle.barcode_calling.exceptions.MissingParameterException;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;

/**
 * Implements the calling algorithm from Section 3.4.4
 * @author Riko Uphoff
 *
 */
public class KMerCaller extends IndexCaller {
	
	private static double DEFAULT_BARCODES_FRACTION_PROBED = 0.01;
	
	private int k;
	
	public KMerCaller(BarcodeDataset barcodes, int k, int maxDist, double indexFractionProbed) {
		super(barcodes, new KMerIndex(k), maxDist, indexFractionProbed);
		
		this.k = k;
	}
	
	public KMerCaller(BarcodeDataset barcodes, int k, int maxDist) {
		this(barcodes, k, maxDist, DEFAULT_BARCODES_FRACTION_PROBED);
	}
	
	public KMerCaller(BarcodeDataset barcodes, int k) {
		this(
				barcodes,
				k, 
				IndexCaller.DEFAULT_DISTANCE_MEASURE.getDefaultMaxDistance(barcodes.getSequenceLength())
		);
	}
	
	public KMerCaller(int k, int maxDist) {
		this(new BarcodeDataset(), k, maxDist);
	}
	
	public KMerCaller(int k) {
		this(k, 0);
	}
	
	public int getK() {
		return this.k;
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
		
		sb.append("_L");
		sb.append(maxDist);
		sb.append("_pos");
		sb.append(this.getBarcodesFractionProbed());
		
		return sb.toString();
	}
	
	public static KMerCaller fromParameters(Map<String, Object> params, BarcodeDataset barcodes) {
		int k;
		
		if (params.containsKey("k")) {
			k = (Integer)params.get("k");
		}
		else {
			throw new MissingParameterException("Parameter 'k' missing!");
		}
		
		KMerCaller alg = new KMerCaller(barcodes, k);
		
		if (params.containsKey("L")) {
			int maxDist = (Integer)params.get("L");
			alg.setMaxDist(maxDist);
		}
		else {
			System.out.println("Notice: Parameter 'L' will be inferred.");
		}
		
		if (params.containsKey("p")) {
			int pos = (Integer)params.get("p");
			alg.setBarcodesFractionProbed((double)pos / (double)barcodes.getSize());
		}
		else {
			System.out.println("Notice: Parameter 'p' will be inferred.");
		}
		
		return alg;
	}

}
