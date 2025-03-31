package de.uni_halle.barcode_calling.util.datasets;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.util.DNASequence;

public class LabeledDataset extends Dataset {
	
	/*
	 * Attribute declarations
	 */
	// Valid label must be a non-negative int. If negative, label is unknown (miss).
	private int[] labels;
	private double pSub = 0, pIns = 0, pDel = 0;
	
	
	/*
	 * Constructors
	 */
	public LabeledDataset(DNASequence[] sequences, int[] labels) {
		super(sequences);
		this.labels = labels;
	}
	
	public LabeledDataset(DNASequence[] sequences) {
		super(sequences);
		this.labels = new int[sequences.length];
		for (int i = 0; i < this.labels.length; i++) {
			this.labels[i] = i;
		}
	}
	
	
	/*
	 * Methods
	 */
	public int[] getLabels() {
		return this.labels;
	}
	
	public int getLabel(int idx) {
		return this.labels[idx];
	}
	
	@Override
	public LabeledDataset getSubset(int subsetSize) {
		DNASequence[] subset = new DNASequence[subsetSize];
		int[] subsetLabels = new int[subsetSize];
		
		for (int i = 0; i < subsetSize; i++) {
			subset[i] = this.getSequence(i);
			subsetLabels[i] = this.getLabel(i);
		}
		
		return new LabeledDataset(subset, subsetLabels);
	}
	
	public void setErrorProbabilities(double pSub, double pIns, double pDel) {
		this.pSub = pSub;
		this.pIns = pIns;
		this.pDel = pDel;
	}
	
	public double[] getErrorProbabilities() {
		double[] out = {pSub, pIns, pDel};
		return out;
	}
	
	public void toFile(Path readsFile, Path labelsFile) {
		this.toFile(readsFile);
		
		// Write .txt file
		List<String> labelsList = new ArrayList<String>();
		for (Integer label : this.labels) {
			labelsList.add(label.toString());
		}
		
		try {
			Files.write(labelsFile, labelsList);
		} catch (IOException e) {
			Main.writeToLog(e);
		}
	}
	
	public void toFile(String readsFile, String labelsFile) {
		this.toFile(Path.of(readsFile), Path.of(labelsFile));
	}
	
//	@Override
//	public int[] sort() {
//		int[] mapping = super.sort();
//		
//		int[] newLabels = new int[this.getSize()];
//		
//		for (int i = 0; i < this.getSize(); i++) {
//			newLabels[i] = this.getLabel(mapping[i]);
//		}
//		
//		return mapping;
//	}
	
	/*
	 * Static methods
	 */
	public static LabeledDataset fromFile(Path readsFile, Path labelsFile) {
		// TODO Handle file format
		if (readsFile.toFile().isDirectory() || labelsFile.toFile().isDirectory()) {
			throw new IllegalArgumentException("One given path is a directory!");
		}
		
		// Read .txt file
		LabeledDataset out = null;
		
		try {
			List<String> sequences = Files.readAllLines(readsFile);
			List<String> labels = Files.readAllLines(labelsFile);
			DNASequence[] sequenceStrings = new DNASequence[sequences.size()];
			int[] labelsInts = new int[labels.size()];
			
			for(int i = 0; i < sequenceStrings.length; i++) {
				sequenceStrings[i] = new DNASequence(sequences.get(i));
				labelsInts[i] = Integer.parseInt(labels.get(i));
			}
			
			out = new LabeledDataset(sequenceStrings, labelsInts);
			
		} catch (FileNotFoundException e) {
			Main.writeToLog(e);
		} catch (IOException e) {
			Main.writeToLog(e);
		}
		
		return out;
	}
	
	public static LabeledDataset fromFile(String readsFile, String labelsFile) {
		Path readsPath = Path.of(readsFile);
		Path labelsPath = Path.of(labelsFile);
		return fromFile(readsPath, labelsPath);
	}
	
}
