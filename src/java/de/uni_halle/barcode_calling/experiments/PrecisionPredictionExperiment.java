package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.DatasetAnalyzer;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;
import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.IndexCaller;

public class PrecisionPredictionExperiment extends CallingExperiment implements Experiment {

	/*
	 * Constants
	 */
	static final int DEFAULT_READS_COUNT = 10000;
	static final double DEFAULT_MAX_MUTATION_RATE = 0.3;
	static final int DEFAULT_REPEATS = 100;
	
	static final String[] CSV_HEADER = {
		"matcher", "p_sub", "p_ins", "p_del", "predicted p_sub", "predicted p_ins", "predicted p_del", "predicted precision", "predicted recall", "precision", "recall", "execution time tuning", "execution time analysis"
	};
	static final String DEFAULT_EXPERIMENT_NAME = "precision-prediction-analysis";
	
	/*
	 * Attribute declarations
	 */
	BarcodeDataset barcodeSet;
	IndexCaller callingAlgorithm;
	String name;
	
	/*
	 * Constructors
	 */
	public PrecisionPredictionExperiment(
			BarcodeDataset barcodeSet,
			IndexCaller callingAlgorithm
	) {
		this.barcodeSet = barcodeSet;
		this.callingAlgorithm = callingAlgorithm;

		this.setExperimentName(DEFAULT_EXPERIMENT_NAME);
	}
	
	public PrecisionPredictionExperiment(
			Path barcodeFile,
			IndexCaller callingAlgorithm
	) {
		this(BarcodeDataset.fromFile(barcodeFile), callingAlgorithm);
	}
	
	/*
	 * Methods
	 */
	@Override
	public String getDefaultResultFileName() {
		StringBuilder sb = new StringBuilder(this.getExperimentName());
		sb.append("_r");
		sb.append(DEFAULT_READS_COUNT);
		sb.append("_n");
		sb.append(this.barcodeSet.getSize());
		sb.append("_");
		sb.append(getTimestampString());
		sb.append(".csv");
		
		return sb.toString();
	}
	
	@Override
	public void doExperiment(Path outPath) {
		
		Main.writeToLog("Starting experiment " + this.getExperimentName());
		
		Path outFilePath = this.resolveOutPath(outPath);
		
		// "matcher", "p_sub", "p_ins", "p_del", "predicted precision", "predicted recall", "precision", "recall", "execution time tuning", "execution time analysis"
		List<Object[]> csvEntries = new ArrayList<>();
		
		for (int i = 0; i < DEFAULT_REPEATS; i++) {
			double pSub = DEFAULT_MAX_MUTATION_RATE/3 * Math.random();
			double pIns = DEFAULT_MAX_MUTATION_RATE/3 * Math.random();
			double pDel = DEFAULT_MAX_MUTATION_RATE/3 * Math.random();
			
			this.analyze(barcodeSet, pSub, pIns, pDel, csvEntries);
			
		}
		
		
		IOMethods.writeCSV(outFilePath, CSV_HEADER, csvEntries);
		Main.writeToLog("Experiment ran successfully!");
	}
	
	private void analyze(BarcodeDataset barcodes, double pSub, double pIns, double pDel, List<Object[]> csvEntries) {
		
		// Generate reads
		LabeledDataset reads = this.barcodeSet.corrupt(pSub, pIns, pDel, DEFAULT_READS_COUNT);
		
		Main.writeToLog("Reads with p_sub = " + pSub + ", p_ins = " + pIns + ", p_del = " + pDel + " created...");
		
		// Measure execution time of analysis
		long start = System.nanoTime();
		
		double[] estimatedErrorProbs = DatasetAnalyzer.estimateErrorProbabilities(reads, barcodes);
		LabeledDataset readsArtif = barcodes.corrupt(pSub, pIns, pDel, DatasetAnalyzer.DEFAULT_K);
		
		long finish = System.nanoTime();
		double timeElapsedAnalysis = (double)(finish - start)/1000000;
		
		Main.writeToLog("> Analysis finished after " + timeElapsedAnalysis + " ms");
		
		// Measure execution time of tuning
		start = System.nanoTime();
		
		double[] matchingRatesPredicted = this.callingAlgorithm.tuneParameters(readsArtif, 0.0, 0.0);
		
		finish = System.nanoTime();
		double timeElapsedTuning = (double)(finish - start)/1000000;
		
		Main.writeToLog("> Tuning finished after " + timeElapsedTuning + " ms");
		
		// Statistics
		int[] matches = this.callingAlgorithm.call(reads);
		double[] matchingRates = CallingExperiment.getMatchingRates(reads.getLabels(), matches);
		
		Main.writeToLog("> Matching rates predicted: " + Arrays.toString(matchingRatesPredicted));
		Main.writeToLog("> Matching rates: " + Arrays.toString(matchingRates));
		
		Object[] entry = {
			this.callingAlgorithm.getName(),
			pSub,
			pIns,
			pDel,
			estimatedErrorProbs[0],
			estimatedErrorProbs[1],
			estimatedErrorProbs[2],
			matchingRatesPredicted[0],
			matchingRatesPredicted[1],
			matchingRates[0],
			matchingRates[1],
			timeElapsedTuning,
			timeElapsedAnalysis,
		};
		
		csvEntries.add(entry);
		
		System.out.println();
	}
	
	
	/*
	 * Factory methods
	 */
	
}
