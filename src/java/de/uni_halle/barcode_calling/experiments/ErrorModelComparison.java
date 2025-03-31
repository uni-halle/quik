package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.LevenshteinBase;
import de.uni_halle.barcode_calling.callers.SequenceLevenshteinBase;
import de.uni_halle.barcode_calling.experiments.data.DNAErrorModel;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.DatasetAnalyzer;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

/**
 * Implements experiment 1 from the thesis.
 * @author Riko Uphoff
 *
 */
public class ErrorModelComparison implements Experiment {
	// TODO
	
	/*
	 * Constants
	 */
	static final DNAErrorModel[] DEFAULT_ERROR_MODELS = {DNAErrorModel.UPHOFF, DNAErrorModel.PRESS};
	static final int DEFAULT_MUTATION_RATE_STEPS = 5;
	static final double DEFAULT_MUTATION_RATE = 0.1;
	static final double DEFAULT_MAX_MUTATION_RATE = 0.25;
	static final int DEFAULT_READS_COUNT = 10000;
	static final int DEFAULT_SEQUENCE_LENGTH_STEPS = 5;
	static final int DEFAULT_SEQUENCE_LENGTH = 34;
	static final int DEFAULT_MAX_SEQUENCE_LENGTH = 50;
	static final String[] CSV_HEADER = {
		"error model", "distance model", "sequence length", "mutation probability", "distance", "average substitutions", "average insertions", "average deletions", "fraction"
	};
	static final String DEFAULT_EXPERIMENT_NAME = "error-model-comparison";
	
	/*
	 * Attribute declarations
	 */
	String name;
	
	/*
	 * Constructors
	 */
	public ErrorModelComparison(
	) {
		this.setExperimentName(DEFAULT_EXPERIMENT_NAME);
	}
	
	
	/*
	 * Methods
	 */
	@Override
	public String getExperimentName() {
		return this.name;
	}

	@Override
	public void setExperimentName(String newName) {
		this.name = newName;
	}
	
	@Override
	public void doExperiment(Path outPath) {
		
		Main.writeToLog("Starting experiment " + this.getExperimentName());
		
		Path outFilePath = this.resolveOutPath(outPath);
		
		// "error model", "distance model", "sequence length", "mutation probability", "distance", "average substitutions", "average insertions", "average deletions", "fraction"
		List<Object[]> csvEntries = new ArrayList<>();
		
		// Iterate over mutation probability
		double stepSizeP = DEFAULT_MAX_MUTATION_RATE/(double)DEFAULT_MUTATION_RATE_STEPS;
		double p;
		
		for (int i = 1; i <= DEFAULT_MUTATION_RATE_STEPS; i++) {
			p = i*stepSizeP;
			this.analyze(p, DEFAULT_SEQUENCE_LENGTH, csvEntries);
		}
		
		// Iterate over sequence length
		int stepSizeLength = DEFAULT_MAX_SEQUENCE_LENGTH/DEFAULT_SEQUENCE_LENGTH_STEPS;
		int sequenceLength;
		
		for (int i = 1; i <= DEFAULT_SEQUENCE_LENGTH_STEPS; i++) {
			sequenceLength = i*stepSizeLength;
			this.analyze(DEFAULT_MUTATION_RATE, sequenceLength, csvEntries);
		}
		
		
		IOMethods.writeCSV(outFilePath, CSV_HEADER, csvEntries);
		Main.writeToLog("Experiment ran successfully!");
	}
	
	@Override
	public String getDefaultResultFileName() {
		StringBuilder sb = new StringBuilder(this.getExperimentName());
		sb.append("_");
		sb.append(getTimestampString());
		sb.append(".csv");
		
		return sb.toString();
	}
	
	private void analyze(double p, int sequenceLength, List<Object[]> csvEntries) {
		DNASequence barcode, read;
		
		for (DNAErrorModel errorModel : DEFAULT_ERROR_MODELS) {
			Map<Integer, double[]> levDistances = new HashMap<Integer, double[]>();
			Map<Integer, double[]> seqLevDistances = new HashMap<Integer, double[]>();
			
			Dataset barcodes = Dataset.random(DEFAULT_READS_COUNT, sequenceLength);
			LabeledDataset reads = barcodes.corrupt(p/3, p/3, p/3, DEFAULT_READS_COUNT, errorModel);
			
			for (int i = 0; i < DEFAULT_READS_COUNT; i++) {
				barcode = barcodes.getSequence(reads.getLabel(i));
				read = reads.getSequence(i);
				
				double[][][] detailedDistMatrix = LevenshteinBase.computeDetailedDistMatrix(read, barcode);
				
				// Levenshtein
				double[] levDistDetail = LevenshteinBase.distanceDetails(detailedDistMatrix);
				int levDist = (int)levDistDetail[0];
				
				if (!levDistances.containsKey(levDist)) {
					levDistances.put(levDist, new double[4]);
				}
				levDistances.get(levDist)[0] += 1;
				levDistances.get(levDist)[1] += levDistDetail[1];
				levDistances.get(levDist)[2] += levDistDetail[2];
				levDistances.get(levDist)[3] += levDistDetail[3];
				
				// Sequence-Levenshtein
				double[] seqLevDistDetail = SequenceLevenshteinBase.distanceDetails(detailedDistMatrix);
				int seqLevDist = (int)seqLevDistDetail[0];
				
				if (!seqLevDistances.containsKey(seqLevDist)) {
					seqLevDistances.put(seqLevDist, new double[4]);
				}
				seqLevDistances.get(seqLevDist)[0] += 1;
				seqLevDistances.get(seqLevDist)[1] += seqLevDistDetail[1];
				seqLevDistances.get(seqLevDist)[2] += seqLevDistDetail[2];
				seqLevDistances.get(seqLevDist)[3] += seqLevDistDetail[3];
			}
			
			// Write Levenshtein
			for (int dist : levDistances.keySet()) {
				int count = (int)levDistances.get(dist)[0];
				
				Object[] entry = {
					errorModel.toString(),
					"Levenshtein",
					sequenceLength,
					p,
					dist,
					(double)levDistances.get(dist)[1] / (double)count,
					(double)levDistances.get(dist)[2] / (double)count,
					(double)levDistances.get(dist)[3] / (double)count,
					(double)count / (double)DEFAULT_READS_COUNT
				};
				csvEntries.add(entry);
			}
			
			// Write Sequence-Levenshtein
			for (int dist : seqLevDistances.keySet()) {
				int count = (int)seqLevDistances.get(dist)[0];
				
				Object[] entry = {
					errorModel.toString(),
					"Sequence-Levenshtein",
					sequenceLength,
					p,
					dist,
					(double)seqLevDistances.get(dist)[1] / (double)count,
					(double)seqLevDistances.get(dist)[2] / (double)count,
					(double)seqLevDistances.get(dist)[3] / (double)count,
					(double)count / (double)DEFAULT_READS_COUNT
				};
				csvEntries.add(entry);
			}
			
		}
		
		
	
		// TODO
		/*
		int barcodeCount = barcodes.getSize();
		int barcodeLength = barcodes.getSequence(0).getLength();
		int editDist = barcodes.getMinMutualEditDist();
		
		// Analyze
		System.out.println("Analyzing index with n = " + barcodes.getSize() + ", l = " + barcodes.getFirstSequenceLength() + ", d = " + barcodes.getMinMutualEditDist() + "...");
		int dist;
		
		for (int a = 0; a < barcodeCount; a++) {
			for (int b = a+1; b < barcodeCount; b++) {
				dist = distance(barcodes.getSequence(a), barcodes.getSequence(b));
				
				Object[] entry = {
					barcodeCount,
					barcodeLength,
					editDist,
					a,
					b,
					dist
				};
				csvEntries.add(entry);
			}
		}
		*/
		
	}
}