package de.uni_halle.barcode_calling.util.datasets;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.IOMethods;

public class ReadDataset extends Dataset {
	
	int barcodeLength;
	
	/*
	 * Constructors
	 */
	
	public ReadDataset(DNASequence[] sequences, String[] tags, int barcodeLength) {
		super(sequences, tags);
		
		this.barcodeLength = barcodeLength;
		
		// Pre checks
		for (DNASequence seq : sequences) {
			if (seq.getLength() < barcodeLength) {
				throw new RuntimeException("Read " + seq + " is shorter that barcode length " + barcodeLength + "!");
			}
		}
	}
	
	public ReadDataset(DNASequence[] sequences, String[] tags) {
		this(sequences, tags, 0);
	}
	
	public ReadDataset(DNASequence[] sequences) {
		this(sequences, dummyTags(sequences.length));
	}
	
	public ReadDataset(DNASequence[] sequences, int barcodeLength) {
		this(sequences, dummyTags(sequences.length), barcodeLength);
	}
	
	public ReadDataset() {
		this(new DNASequence[0]);
	}
	
	public ReadDataset(String[] sequences, String[] tags, int barcodeLength) {
		this(stringsToSequences(sequences), tags, barcodeLength);
	}
	
	public ReadDataset(String[] sequences, String[] tags) {
		this(sequences, tags, 0);
	}
	
	public ReadDataset(String[] sequences) {
		this(sequences, dummyTags(sequences.length));
	}
	
	public ReadDataset(String[] sequences, int barcodeLength) {
		this(sequences, dummyTags(sequences.length), barcodeLength);
	}
	
	public ReadDataset(Dataset data, int barcodeLength) {
		this(data.getSequences(), data.getTags(), barcodeLength);
	}
	
	public ReadDataset(Dataset data) {
		this(data.getSequences(), data.getTags());
	}
	
	/*
	 * Methods
	 */
	
	public int getBarcodeLength() {
		return this.barcodeLength;
	}
	
	public ReadDataset getSubset(int subsetSize) {
		Dataset rawSub = super.getSubset(subsetSize);
		return new ReadDataset(rawSub, this.getBarcodeLength());
	}
	
	public ReadDataset getSubset(int[] subsetIds) {
		Dataset rawSub = super.getSubset(subsetIds);
		return new ReadDataset(rawSub, this.getBarcodeLength());
	}
	
	public ReadDataset getSubset(List<Integer> subsetIds) {
		Dataset rawSub = super.getSubset(subsetIds);
		return new ReadDataset(rawSub, this.getBarcodeLength());
	}
	
	/*
	 * Static methods
	 */
	
	public static ReadDataset fromFile(Path filePath, int barcodeLength) {
		String extension = IOMethods.getFileExtension(filePath);
		
		// Read file
		ReadDataset out = null;
		
		try {
			if(extension == null || extension.equals(".txt")) {
				/*
				 * One sequence per line expected
				 */
				List<String> sequences = Files.readAllLines(filePath);
				String[] sequenceStrings = new String[sequences.size()];
				
				out = new ReadDataset(sequences.toArray(sequenceStrings));
			}
			else if (extension.equals(".fastq")) {
				/*
				 * Typical fastq format:
				 * 1. @TITLE
				 * 2. DNA_STRING
				 * 3. +
				 * 4. QUALITY_ASSESMENT
				 */
				List<DNASequence> sequences = new ArrayList<>();
				List<String> tags = new ArrayList<>();
				
				try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()))) {
					for (String line = reader.readLine(); line != null; line = reader.readLine()) {
						if (line.charAt(0) == '@') {
							String tag = line;				// Keep title
							line = reader.readLine();
							DNASequence seq = new DNASequence(line);		// Keep DNA string
							seq = seq.subSequence(0, barcodeLength);
							
							if (seq.getLength() == barcodeLength) {
								tags.add(tag);
								sequences.add(seq);
							}
							else {
								System.out.println(String.format("Skipping read %s as it is too short after removing Ns...", tag));
							}
							
							line = reader.readLine();	// Remove +
							line = reader.readLine();	// Remove quality
						}
					}
					
					String[] refString = new String[0];
					DNASequence[] refSeq = new DNASequence[0];
					
					out = new ReadDataset(sequences.toArray(refSeq), tags.toArray(refString), barcodeLength);
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
	
	public static Dataset fromFile(String filePath, int barcodeLength) {
		Path path = Path.of(filePath);
		return fromFile(path, barcodeLength);
	}
	
	public static Dataset fromFile(Path filePath) {
		return fromFile(filePath, 0);
	}
	
	public static Dataset fromFile(String filePath) {
		return fromFile(filePath, 0);
	}
}
