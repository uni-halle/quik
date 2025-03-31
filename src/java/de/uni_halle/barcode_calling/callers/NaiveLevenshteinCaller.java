package de.uni_halle.barcode_calling.callers;

import java.util.Map;

import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;

/**
 * Implements the first baseline from the paper with Levenshtein distance (see Section 3.4.1)
 * @author Riko Uphoff
 *
 */
public class NaiveLevenshteinCaller extends NaiveCaller {
	
	/*
	 * Constructors
	 */
	public NaiveLevenshteinCaller(BarcodeDataset barcodes) {
		super(barcodes, new LevenshteinBase());
	}
	
	public NaiveLevenshteinCaller() {
		super(new LevenshteinBase());
	}
	
	public static NaiveLevenshteinCaller fromParameters(Map<String, Object> params, BarcodeDataset barcodes) {
		NaiveLevenshteinCaller alg = new NaiveLevenshteinCaller(barcodes);
		
		return alg;
	}

}
