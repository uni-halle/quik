package de.uni_halle.barcode_calling.callers;

import java.util.Map;

import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;

/**
 * Implements the first baseline from the paper with SL distance (see Section 3.4.1)
 * @author Riko Uphoff
 *
 */
public class NaiveSequenceLevenshteinCaller extends NaiveCaller {

	/*
	 * Constructors
	 */
	public NaiveSequenceLevenshteinCaller(BarcodeDataset barcodes) {
		super(barcodes, new SequenceLevenshteinBase());
	}
	
	public NaiveSequenceLevenshteinCaller() {
		super(new SequenceLevenshteinBase());
	}
	
	public static NaiveSequenceLevenshteinCaller fromParameters(Map<String, Object> params, BarcodeDataset barcodes) {
		NaiveSequenceLevenshteinCaller alg = new NaiveSequenceLevenshteinCaller(barcodes);
		
		return alg;
	}
	
}
