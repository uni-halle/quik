package de.uni_halle.barcode_calling.util.datasets;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Supplier;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.experiments.data.DNAErrorModel;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.Tuples.SequenceIndexTuple;

public class Dataset {
	
	/*
	 * Constants
	 */
	public enum RandomDistribution {
		CONSTANT, UNIFORM, NORMAL;
	}
	
	
	/*
	 * Attribute declarations
	 */
	private DNASequence[] sequences;
	private String[] tags;
	
	
	/*
	 * Constructors
	 */
	public Dataset(DNASequence[] sequences, String[] tags) {
		this.sequences = sequences;
		this.tags = tags;
	}
	
	public Dataset(DNASequence[] sequences) {
		this(sequences, dummyTags(sequences.length));
	}
	
	public Dataset(String[] sequenceStrings, String[] tags) {
		this(stringsToSequences(sequenceStrings), tags);
	}
	
	public Dataset(String[] sequenceStrings) {
		this(stringsToSequences(sequenceStrings));
	}
	
	public Dataset() {
		this(new DNASequence[0]);
	}
	
	
	/*
	 * Methods
	 */
	public DNASequence[] getSequences() {
		return this.sequences;
	}
	
	public void setSequences(DNASequence[] sequences) {
		this.sequences = sequences;
	}
	
	public DNASequence getSequence(int idx) {
		return this.getSequences()[idx];
	}
	
	public String[] getTags() {
		return this.tags;
	}
	
	public String getTag(int idx) {
		return idx >= 0 ? this.tags[idx] : "-";
	}
	
	public int getSize() {
		return this.sequences != null ? this.sequences.length : 0;
	}
	
	public Dataset getSubset(int subsetSize) {
		DNASequence[] subset = new DNASequence[subsetSize];
		String[] tagsSub = new String[subsetSize];
		
		for (int i = 0; i < subsetSize; i++) {
			subset[i] = this.getSequence(i);
			tagsSub[i] = this.getTag(i);
		}
		
		return new Dataset(subset, tagsSub);
	}
	
	public Dataset getSubset(int[] subsetIds) {
		int subsetSize = subsetIds.length;
		DNASequence[] subset = new DNASequence[subsetSize];
		String[] tagsSub = new String[subsetSize];
		
		for (int i = 0; i < subsetSize; i++) {
			subset[i] = this.getSequence(subsetIds[i]);
			tagsSub[i] = this.getTag(subsetIds[i]);
		}
		
		return new Dataset(subset, tagsSub);
	}
	
	public Dataset getSubset(List<Integer> subsetIds) {
		int subsetSize = subsetIds.size();
		DNASequence[] subset = new DNASequence[subsetSize];
		String[] tagsSub = new String[subsetSize];
		
		for (int i = 0; i < subsetSize; i++) {
			subset[i] = this.getSequence(subsetIds.get(i));
			tagsSub[i] = this.getTag(subsetIds.get(i));
		}
		
		return new Dataset(subset, tagsSub);
	}
	
	public void toFile(Path filePath) {
		// TODO Handle file format
		String filePathString = filePath.toString();
		if(!filePathString.endsWith(".txt")) {
			throw new IllegalArgumentException("Wrong file extension");
		}
		
		// Write .txt file
		List<String> seqStrings = new ArrayList<String>();
		for (DNASequence seq : this.sequences) {
			seqStrings.add(seq.toString());
		}
		
		try {
			Files.write(filePath, seqStrings);
		} catch (IOException e) {
			Main.writeToLog(e);
		}
	}
	
	public void toFile(String filePath) {
		Path path = Path.of(filePath);
		this.toFile(path);
	}
	
	public LabeledDataset corruptPress(
			double pSub,
			double pIns,
			double pDel,
			int readsCount
	) {
		return this.corrupt(pSub, pIns, pDel, readsCount, RandomDistribution.CONSTANT, DNAErrorModel.PRESS);
	}
	
	public LabeledDataset corrupt(
			double pSub,
			double pIns,
			double pDel,
			int readsCount,
			RandomDistribution abundanceDistribution,
			DNAErrorModel errorModel
	) {
		// TODO Refactor
		int expectedReadsPerBarcode = (int)Math.ceil((double)readsCount / (double)this.getSize());
		
		List<DNASequence> readsList = new ArrayList<>();
		List<Integer> labelsList = new ArrayList<>();
		// Generates (random) numbers of corrupted reads
		Supplier<Integer> sizeGenerator = getSizeGenerator(expectedReadsPerBarcode, abundanceDistribution);
		int readCount;
		
		//Corrupt
		for (int i = 0; i < this.sequences.length; i++) {
			DNASequence sequence = this.sequences[i];
			readCount = sizeGenerator.get();
			
			for(int j = 0; j < readCount; j++) {
				readsList.add(sequence.corrupt(pSub, pIns, pDel, errorModel));
				labelsList.add(i);
			}
		}
		
		// Shuffle dataset and labels coherently
		long seed = System.currentTimeMillis();
		Collections.shuffle(readsList, new Random(seed));
		Collections.shuffle(labelsList, new Random(seed));
		// Lists to arrays
		DNASequence[] reads = readsList.toArray(new DNASequence[0]);
		int[] labels = labelsList.stream().mapToInt(i->i).toArray();
		LabeledDataset dataset = new LabeledDataset(reads, labels);
		dataset.setErrorProbabilities(pSub, pIns, pDel);
		
		// TODO Only works for CONSTANT abundance
		return dataset.getSubset(readsCount);
	}
	
	public LabeledDataset corrupt(
			double pSub,
			double pIns,
			double pDel,
			int readsCount,
			RandomDistribution abundanceDistribution
	) {
		return this.corrupt(pSub, pIns, pDel, readsCount, abundanceDistribution, DNAErrorModel.UPHOFF);
	}
	
	public LabeledDataset corrupt(
			double pSub,
			double pIns,
			double pDel,
			int readsCount
	) {
		return this.corrupt(pSub, pIns, pDel, readsCount, RandomDistribution.CONSTANT);
	}
	
	public LabeledDataset corrupt(
			double p,
			int readsCount,
			RandomDistribution abundanceDistribution
	) {
		return this.corrupt(p/3.0, p/3.0, p/3.0, readsCount, abundanceDistribution);
	}
	
	public LabeledDataset corrupt(
			double p,
			int readsCount
	) {
		return this.corrupt(p, readsCount, RandomDistribution.CONSTANT);
	}
	
	public LabeledDataset corrupt(
			double pSub,
			double pIns,
			double pDel,
			int readsCount,
			DNAErrorModel errorModel
	) {
		switch (errorModel) {
		case PRESS:
			return this.corruptPress(pSub, pIns, pDel, readsCount);
		case UPHOFF:
			return this.corrupt(pSub, pIns, pDel, readsCount);
		default:
			return null;
		}
	}
	
//	/**
//	 * Sorts the underlying sequences lexicographically and return a mapping to the old ordering.
//	 * @return
//	 */
//	public int[] sort() {
//		int seqCount = this.getSize();
//		
//		SequenceIndexTuple[] seqTuples = new SequenceIndexTuple[seqCount];
//		for (int i = 0; i < seqCount; i++) {
//			seqTuples[i] = new SequenceIndexTuple(this.getSequence(i), i);
//		}
//		
//		Arrays.parallelSort(seqTuples);
//		
//		int[] mapping = new int[seqCount];
//		DNASequence[] seqs = new DNASequence[seqCount];
//		String[] tags = new String[seqCount];
//		
//		for (int i = 0; i < seqCount; i++) {
//			seqs[i] = seqTuples[i].seq;
//			tags[i] = this.tags[seqTuples[i].idx];
//			mapping[i] = seqTuples[i].idx;
//		}
//		
//		this.sequences = seqs;
//		this.tags = tags;
//		
//		return mapping;
//	}
	
	public void shuffle() {
		Random rnd = new Random();
		DNASequence tempS;
		String tempT;
		
	    for (int i = 1; i < this.getSize(); i++)
	    {
	      int index = rnd.nextInt(i + 1);
	      // Simple swap
	      tempS = this.sequences[index];
	      this.sequences[index] = this.sequences[i];
	      this.sequences[i] = tempS;
	      
	      tempT = this.tags[index];
	      this.tags[index] = this.tags[i];
	      this.tags[i] = tempT;
	    }
	}
	
	@Override
	public String toString() {
		StringBuffer b = new StringBuffer(this.sequences[0].toString());
		
		for (int i = 1; i < this.sequences.length; i++) {
			b.append("\n");
			b.append(this.sequences[i].toString());
		}
		
		return b.toString();
	}
	

	/*
	 * Static methods
	 */
	
	public static DNASequence[] stringsToSequences(String[] sequenceStrings) {
		List<DNASequence> sequenceList = new ArrayList<>();
		
		for (int i = 0; i < sequenceStrings.length; i++) {
			try {
				DNASequence seq = new DNASequence(sequenceStrings[i]);
				sequenceList.add(seq);
			}
			catch (IllegalArgumentException e) {
				// Unknown base / only Ns contained / empty string
				// TODO This may break automatic labeling
				Main.writeToLog(e);
				Main.writeToLog("Skipping sequence at position " + i);
			}
		}
		
		return sequenceList.toArray(new DNASequence[0]);
	}
	
	public static String[] dummyTags(int size) {
		String[] tags = new String[size];
		
		for (int i = 0; i < size; i++)
			tags[i] = "" + i;
		
		return tags;
	}
	
	/**
	 * Returns a set of n random sequences with length l
	 * @param n Number of random sequences to be returned
	 * @param l Length of resulting sequences
	 * @return Set on random sequences
	 */
	public static Dataset random(int n, int l) {
		DNASequence[] barcodeSequences = new DNASequence[n];
		
		for (int i = 0; i < n; i++)
			barcodeSequences[i] = DNASequence.random(l);
		
		Dataset barcodes = new Dataset(barcodeSequences);
		
		return barcodes;
	}

	public static Dataset[] datasetsFromPaths(Path[] files) {
		Dataset[] datasets = new Dataset[files.length];
		Dataset curDataset;
		
		for (int i = 0; i < files.length; i++) {
			curDataset = fromFile(files[i]);
			datasets[i] = curDataset;
		}
		
		return datasets;
	}
	
	public static Dataset fromFile(Path filePath) {
		String extension = IOMethods.getFileExtension(filePath);
		
		// Read .txt file
		Dataset out = null;
		
		try {
			if(extension == null || extension.equals(".txt")) {
				/*
				 * One sequence per line expected
				 */
				List<String> sequences = Files.readAllLines(filePath);
				String[] sequenceStrings = new String[sequences.size()];
				
				out = new Dataset(sequences.toArray(sequenceStrings));
			}
			else if (extension.equals(".fastq")) {
				/*
				 * Typical fastq format:
				 * 1. @TITLE
				 * 2. DNA_STRING
				 * 3. +
				 * 4. QUALITY_ASSESMENT
				 */
				List<String> sequences = new ArrayList<>();
				
				try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()))) {
					for (String line = reader.readLine(); line != null; line = reader.readLine()) {
						if (line.charAt(0) == '@') {
							line = reader.readLine();	// Remove title
							sequences.add(line);		// Keep DNA string
							line = reader.readLine();	// Remove +
							line = reader.readLine();	// Remove quality
						}
					}
					
					String[] sequenceStrings = new String[sequences.size()];
					
					out = new Dataset(sequences.toArray(sequenceStrings));
				}
				
				
			}
			else {
				throw new IllegalArgumentException("Wrong file extension " + extension + " for path " + filePath);
			}
			
			
		} catch (FileNotFoundException e) {
			Main.writeToLog(e);
		} catch (IOException e) {
			Main.writeToLog(e);
		}
		
		return out;
	}
	
	public static Dataset fromFile(String filePath) {
		Path path = Path.of(filePath);
		return fromFile(path);
	}
	
	private static Supplier<Integer> getSizeGenerator(
			int expectedReadsPerBarcode,
			RandomDistribution abundanceDistribution
	) {
		// Init random generator for count of corrupted reads for a barcode
		switch (abundanceDistribution) {
		case CONSTANT:
			return () -> expectedReadsPerBarcode;
		case UNIFORM:
			return () -> (int)(Math.random() * (2*expectedReadsPerBarcode + 1));  // TODO correct?
		case NORMAL:
			// TODO
			throw new IllegalArgumentException("Normal distribution not yet implemented!");
		default:
			throw new IllegalArgumentException("Invalid RandomDistribution!");
		}
	}
	
}
