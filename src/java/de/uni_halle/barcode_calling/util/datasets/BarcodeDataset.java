package de.uni_halle.barcode_calling.util.datasets;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.MissingResourceException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.exceptions.MissingParameterException;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.IOMethods;

public class BarcodeDataset extends Dataset {
	
	private int minMutualEditDist;
	private int length;
	
	/*
	 * Constructors
	 */
	
	public BarcodeDataset(DNASequence[] sequences, String[] tags, int minMutualEditDist) {
		super(sequences, tags);
		this.minMutualEditDist = minMutualEditDist;
		
		// All barcodes are expected to have the same length
		this.length = this.getSize() > 0 ? this.getSequence(0).getLength() : 0;
		
		// Pre checks
		for (DNASequence seq : sequences) {
			if (seq.getLength() != this.length) {
				throw new RuntimeException("Sequence " + seq + " does not have expected length " + this.length + "!");
			}
		}
	}
	
	public BarcodeDataset(DNASequence[] sequences, String[] tags) {
		this(sequences, tags, 0);
	}
	
	public BarcodeDataset(DNASequence[] sequences) {
		this(sequences, dummyTags(sequences.length));
	}
	
	public BarcodeDataset(DNASequence[] sequences, int minMutualEditDist) {
		this(sequences, dummyTags(sequences.length), minMutualEditDist);
	}
	
	public BarcodeDataset() {
		this(new DNASequence[0]);
	}
	
	public BarcodeDataset(String[] sequences, String[] tags, int minMutualEditDist) {
		this(stringsToSequences(sequences), tags, minMutualEditDist);
	}
	
	public BarcodeDataset(String[] sequences, String[] tags) {
		this(sequences, tags, 0);
	}
	
	public BarcodeDataset(String[] sequences) {
		this(sequences, dummyTags(sequences.length));
	}
	
	public BarcodeDataset(String[] sequences, int minMutualEditDist) {
		this(sequences, dummyTags(sequences.length), minMutualEditDist);
	}
	
	public BarcodeDataset(Dataset data, int minMutualEditDist) {
		this(data.getSequences(), data.getTags(), minMutualEditDist);
	}
	
	public BarcodeDataset(Dataset data) {
		this(data.getSequences(), data.getTags());
	}
	
	/*
	 * Methods
	 */
	
	public int getMinMutualEditDist() {
		return this.minMutualEditDist;
	}
	
	public void setMinMutualEditDist(int dist) {
		this.minMutualEditDist = dist;
	}
	
	public int getSequenceLength() {
		return this.length;
	}
	
	public BarcodeDataset getSubset(int subsetSize) {
		Dataset rawSub = super.getSubset(subsetSize);
		return new BarcodeDataset(rawSub, this.getMinMutualEditDist());
	}
	
	public BarcodeDataset getSubset(int[] subsetIds) {
		Dataset rawSub = super.getSubset(subsetIds);
		return new BarcodeDataset(rawSub, this.getMinMutualEditDist());
	}
	
	public BarcodeDataset getSubset(List<Integer> subsetIds) {
		Dataset rawSub = super.getSubset(subsetIds);
		return new BarcodeDataset(rawSub, this.getMinMutualEditDist());
	}
	
	
	/*
	 * Static methods
	 */
	public static String cleanToken(String token) {
		token = token.strip();
		// TODO Handle quotes properly
		token = token.replace("\"", "");
		token = token.replace("'", "");
		return token;
	}
	
	public static BarcodeDataset fromFile(Path filePath) {
		String extension = IOMethods.getFileExtension(filePath);
		
		// Read file
		BarcodeDataset out = null;
		
		try {
			if(extension == null || extension.equals(".txt")) {
				/*
				 * One sequence per line expected
				 */
				List<String> sequences = Files.readAllLines(filePath);
				String[] sequenceStrings = new String[sequences.size()];
				
				out = new BarcodeDataset(sequences.toArray(sequenceStrings));
			}
			else if (extension.equals(".list")) {
				/*
				 * Expected .list format:
				 * bc.ID	pattern	bc	freq
				 */
				List<String> sequences = new ArrayList<>();
				List<String> tags = new ArrayList<>();
				
				try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()))) {
					reader.readLine();	// Skip header
					
					for (String line = reader.readLine(); line != null; line = reader.readLine()) {
						String[] lineSplit = line.split("\t");
						tags.add(lineSplit[0]);
						sequences.add(lineSplit[2]);
					}
					
					String[] ref = new String[0];
					
					out = new BarcodeDataset(sequences.toArray(ref), tags.toArray(ref));
				}
				
				
			}
			else if (extension.equals(".csv")) {
				/*
				 * Expects separator ",". Quotes are not yet supported. Expected column names:
				 * id: ID of barcode (optional).
				 * sequence: Sequence of barcode.
				 * TODO Make more flexible.
				 */
				List<String> sequences = new ArrayList<>();
				List<String> tags = new ArrayList<>();
				
				try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()))) {
					int idColumn = -1, sequenceColumn = -1;
					
					// Read header
					String line = reader.readLine();
					String[] tokens = line.split(",");
					
					for (int i = 0; i < tokens.length; i++) {
						String token = tokens[i];
						token = cleanToken(token);
						
						
						if (token.equals("id")) {
							idColumn = i;
						}
						else if (token.equals("sequence")) {
							sequenceColumn = i;
						}
					}
					
					if (sequenceColumn < 0)
						throw new IOException("Required column 'sequence' not found in input .csv file!");
					
					// Read lines
					for (line = reader.readLine(); line != null; line = reader.readLine()) {
						tokens = line.split(",");
						
						sequences.add(cleanToken(tokens[sequenceColumn]));
						
						if (idColumn >= 0)
							tags.add(cleanToken(tokens[idColumn]));
					}
					
					// Init data set
					String[] ref = new String[0];
					
					if (idColumn >= 0)
						out = new BarcodeDataset(sequences.toArray(ref), tags.toArray(ref));
					else
						out = new BarcodeDataset(sequences.toArray(ref));
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
	
	public static BarcodeDataset fromFile(String filePath) {
		Path path = Path.of(filePath);
		return fromFile(path);
	}
	
	public static int getEditDistFromPath(Path file) {
		int editDist = 0;
		Pattern distPattern = Pattern.compile("_d(\\d+)");
		Matcher distMatcher;
		
		String filename = file.toFile().getName();
			
		// Get mutual edit distance from file name
		distMatcher = distPattern.matcher(filename);
		
		if (distMatcher.find()) {
		    editDist = Integer.parseInt(distMatcher.group(1));
		}
		
		return editDist;
	}
	
	public static BarcodeDataset[] datasetsFromPaths(Path[] files) {
		BarcodeDataset[] datasets = new BarcodeDataset[files.length];
		BarcodeDataset curDataset;
		int editDist;
		
		for (int i = 0; i < files.length; i++) {
			editDist = getEditDistFromPath(files[i]);
			
			curDataset = fromFile(files[i]);
			curDataset.setMinMutualEditDist(editDist);
			datasets[i] = curDataset;
		}
		
		return datasets;
	}
}
