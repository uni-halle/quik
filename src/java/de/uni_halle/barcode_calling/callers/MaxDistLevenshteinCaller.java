package de.uni_halle.barcode_calling.callers;

import java.util.Map;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

/**
 * Implements the second baseline from the paper with Levenshtein distance (see Section 3.4.1)
 * @author Riko Uphoff
 *
 */
public class MaxDistLevenshteinCaller extends MaxDistanceCaller {
	
	/*
	 * Constructors
	 */
	public MaxDistLevenshteinCaller(BarcodeDataset barcodes) {
		super(
			barcodes,
			new LevenshteinBase()
		);
	}
	
	@Override
	public double[] tuneParameters(LabeledDataset readsArtif, double minPrec, double minRec) {
		Main.writeToLog("No tuning needed for matcher " + this.getName());
		return new double[3];
	}
	
	public static MaxDistLevenshteinCaller fromParameters(Map<String, Object> params, BarcodeDataset barcodes) {
		MaxDistLevenshteinCaller alg = new MaxDistLevenshteinCaller(barcodes);
		
		if (params.containsKey("L")) {
			int maxDist = (Integer)params.get("L");
			alg.setMaxDist(maxDist);
		}
		
		return alg;
	}
}
