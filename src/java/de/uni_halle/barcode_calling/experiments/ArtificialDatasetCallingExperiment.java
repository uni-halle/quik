package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.experiments.CallingExperiment.getRepresentativeBarcodeSet;
import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.BarcodeCaller;
import de.uni_halle.barcode_calling.callers.IndexCaller;
import de.uni_halle.barcode_calling.callers.KMerCaller;
import de.uni_halle.barcode_calling.callers.MaxDistLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.MaxDistSequenceLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.NaiveLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.NaiveSequenceLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.TrimerCaller;
import de.uni_halle.barcode_calling.callers.indexes.KMerIndex;
import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.DatasetAnalyzer;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset.RandomDistribution;

/**
 * Implements experiment 3 of the thesis.
 * @author Riko Uphoff
 *
 */
public class ArtificialDatasetCallingExperiment extends CallingExperiment {
	/*
	 * Constants
	 */
	static final RandomDistribution DEFAULT_ABUNDANCE_DISTRIBUTION = RandomDistribution.CONSTANT;
	static final int DEFAULT_MUTATION_RATE_STEPS = 6;
	static final double DEFAULT_MUTATION_RATE = 0.1;
	static final double DEFAULT_MAX_MUTATION_RATE = 0.3;
	static final int DEFAULT_READS_COUNT = 1000;
	static final int DEFAULT_BARCODE_COUNT = 1000;
	static final int DEFAULT_BARCODE_LENGTH = 34;
	static final int DEFAULT_BARCODE_DIST = 9;
	static final double DEFAULT_MIN_PRECISION = 0.8;
	static final double DEFAULT_MIN_RECALL = 0.8;
	static final String[] CSV_HEADER = {
		"caller", "barcode count", "barcode length", "mutual edit distance", "reads count", "mutation probability", "execution time calling", "execution time tuning", "precision", "recall", "f1", "predicted precision", "predicted recall", "predicted f1", "L", "pos"
	};
	static final String DEFAULT_EXPERIMENT_NAME = "artificial-dataset-calling-comparison";
	
	/*
	 * Attribute declarations
	 */
	protected BarcodeDataset[] barcodeSets;
	protected double[] mutationRates;
	protected int readsCount;
	protected RandomDistribution abundanceDistribution;
	protected BarcodeCaller[] candidates;
	
	/*
	 * Constructors
	 */
	public ArtificialDatasetCallingExperiment(
			BarcodeDataset[] barcodeSets,
			BarcodeCaller[] candidates,
			int readsCount,
			double[] mutationRates,
			RandomDistribution abundanceDistribution
	) {
		this.barcodeSets = barcodeSets;
		this.candidates = candidates;
		this.mutationRates = mutationRates;
		this.readsCount = readsCount;
		this.abundanceDistribution = abundanceDistribution;
		this.setExperimentName(DEFAULT_EXPERIMENT_NAME);
	}
	
	public ArtificialDatasetCallingExperiment(
			BarcodeDataset[] barcodeSets,
			BarcodeCaller[] candidates,
			int readsCount,
			int mutationRateSteps,
			double maxMutationRate
	) {
		this(barcodeSets, candidates, readsCount, null, DEFAULT_ABUNDANCE_DISTRIBUTION);
		
		// Set mutation rates
		this.mutationRates = new double[mutationRateSteps+1];
		double stepSize = maxMutationRate/(double)mutationRateSteps;
		
		for (int i = 0; i <= mutationRateSteps; i++) {
			this.mutationRates[i] = i*stepSize;
		}
	}
	
	public ArtificialDatasetCallingExperiment(
			BarcodeDataset[] barcodeSets,
			BarcodeCaller[] candidates
	) {
		this(
				barcodeSets, 
				candidates, 
				DEFAULT_READS_COUNT,
				DEFAULT_MUTATION_RATE_STEPS, 
				DEFAULT_MAX_MUTATION_RATE
		);
	}
	
	
	/*
	 * Methods
	 */
	
	@Override
	public String getDefaultResultFileName() {
		StringBuilder sb = new StringBuilder(this.getExperimentName());
		sb.append("_r");
		sb.append(this.readsCount);
		sb.append("_");
		sb.append(getTimestampString());
		sb.append(".csv");
		
		return sb.toString();
	}
	
	@Override
	public void doExperiment(Path outPath) {
		
		Main.writeToLog("Starting experiment " + this.getExperimentName());
		
		Path outFilePath = this.resolveOutPath(outPath);
		
		// "caller", "barcode count", "barcode length", "mutual edit distance", "reads count", "mutation probability", "execution time", "recall", "precision", "f1"
		List<Object[]> csvEntries = new ArrayList<>();
		
		// Iterate over barcode-specific parameters
		if (this.barcodeSets.length > 1) {
			for (BarcodeDataset barcodes : this.barcodeSets) {
				for (BarcodeCaller candidate : this.candidates) {
					candidate.setBarcodes(barcodes);
				}
					
				this.doCalling(barcodes, DEFAULT_MUTATION_RATE, DEFAULT_READS_COUNT, csvEntries);
			}
		}
		
		// Set default barcode set
		BarcodeDataset defaultBarcodeSet = getRepresentativeBarcodeSet(this.barcodeSets, DEFAULT_BARCODE_COUNT, DEFAULT_BARCODE_LENGTH, DEFAULT_BARCODE_DIST);
		for (BarcodeCaller candidate : this.candidates) {
			candidate.setBarcodes(defaultBarcodeSet);
		}
		
		// Iterate over read-specific parameters
		for (double p : this.mutationRates) {
			this.doCalling(defaultBarcodeSet, p, this.readsCount, csvEntries);
		}
		
		IOMethods.writeCSV(outFilePath, CSV_HEADER, csvEntries);
		Main.writeToLog("Experiment ran successfully!");
	}
	
	protected void doCalling(
			BarcodeDataset barcodes, double p, int readsCount, List<Object[]> csvEntries
		) {
		
		int barcodeCount = barcodes.getSize();
		int barcodeLength = barcodes.getSequenceLength();
		int editDist = barcodes.getMinMutualEditDist();
		
		// Generate reads
		LabeledDataset reads = barcodes.corrupt(p, readsCount, this.abundanceDistribution);
		LabeledDataset readsArtif = DatasetAnalyzer.generateSimilarReads(reads, barcodes);
		Main.writeToLog("Corrupted dataset with " + readsCount + " reads and p = " + p + " created...");
		
		// Call
		Main.writeToLog("Calling with n = " + barcodes.getSize() + ", l = " + barcodes.getSequenceLength() + ", d = " + barcodes.getMinMutualEditDist() + "...");
		
		for (BarcodeCaller candidate : this.candidates) {
			// Call
			long start = System.nanoTime();
			
			double[] predictions = candidate.tuneParameters(readsArtif, DEFAULT_MIN_PRECISION, DEFAULT_MIN_RECALL);
			
			long finish = System.nanoTime();
			double timeTuning = (double)(finish - start)/1000000;
			
			Main.writeToLog("> Caller \"" + candidate.getName() + "\" finished tuning in " + timeTuning + " ms...");
			
			// Call
			start = System.nanoTime();
			
			int[] calles = candidate.callParallel(reads.getSequences());
			
			finish = System.nanoTime();
			double timeCalling = (double)(finish - start)/1000000;
			
			Main.writeToLog("> Caller \"" + candidate.getName() + "\" finished calling in " + timeCalling + " ms...");
			
			// Statistics
			double[] callingRates = getMatchingRates(reads.getLabels(), calles);
			
			Object[] entry = {
				candidate.getName(),
				barcodeCount,
				barcodeLength,
				editDist,
				readsCount,
				p,
				timeCalling / (double)readsCount,	// Millisecond per read
				timeTuning,
				callingRates[0],
				callingRates[1],
				callingRates[2],
				predictions[0],
				predictions[1],
				predictions[2],
				candidate.getMaxDist(),
				candidate.getBarcodesFractionProbed(),
			};
			
			csvEntries.add(entry);
		}
		
		System.out.println();
	}
	
	
	/*
	 * Factory methods
	 */
	public static ArtificialDatasetCallingExperiment createBaselineComparison(
			Path[] barcodesFiles
	) {
		BarcodeDataset dummy = new BarcodeDataset();
		BarcodeCaller[] candidates = {
				new NaiveLevenshteinCaller(dummy),
				new NaiveSequenceLevenshteinCaller(dummy),
				new MaxDistLevenshteinCaller(dummy),
				new MaxDistSequenceLevenshteinCaller(dummy),
				new KMerCaller(4),
				new KMerCaller(5),
				new KMerCaller(6),
				new KMerCaller(7),
				new TrimerCaller(dummy),
		};
		
		BarcodeDataset[] barcodeSets = BarcodeDataset.datasetsFromPaths(barcodesFiles);
		ArtificialDatasetCallingExperiment exp = new ArtificialDatasetCallingExperiment(
				barcodeSets, 
				candidates, 
				1000, 
				DEFAULT_MUTATION_RATE_STEPS, 
				DEFAULT_MAX_MUTATION_RATE
		);
		exp.setExperimentName("artificial-dataset-calling-comparison_baseline");
		
		return exp;
	}
	
	public static ArtificialDatasetCallingExperiment createNonBaselineComparison(
			Path[] barcodesFiles
	) {
		BarcodeDataset dummy = new BarcodeDataset();
		BarcodeCaller[] candidates = {
				new KMerCaller(0, 4),
				new KMerCaller(0, 5),
				new KMerCaller(0, 6),
				new KMerCaller(0, 7),
				new TrimerCaller(dummy),
		};
		
		BarcodeDataset[] barcodeSets = BarcodeDataset.datasetsFromPaths(barcodesFiles);
		ArtificialDatasetCallingExperiment exp = new ArtificialDatasetCallingExperiment(
				barcodeSets, 
				candidates, 
				10000, 
				DEFAULT_MUTATION_RATE_STEPS, 
				DEFAULT_MAX_MUTATION_RATE
		);
		exp.setExperimentName("artificial-dataset-calling-comparison_non-baseline");
		
		return exp;
	}
	
	public static ArtificialDatasetCallingExperiment createLargeScaleComparison(
			Path[] barcodesFiles
	) {
		
		BarcodeCaller[] candidates = {
				new KMerCaller(0, 4),
				new KMerCaller(0, 5),
				new KMerCaller(0, 6),
				new KMerCaller(0, 7),
		};
		
		BarcodeDataset[] barcodeSets = BarcodeDataset.datasetsFromPaths(barcodesFiles);
		double[] mutationRates = {0.1};
		
		ArtificialDatasetCallingExperiment exp = new ArtificialDatasetCallingExperiment(
				barcodeSets, 
				candidates, 
				1000000, 
				mutationRates,
				DEFAULT_ABUNDANCE_DISTRIBUTION
		);
		exp.setExperimentName("artificial-dataset-calling-comparison_large-scale");
		
		return exp;
	}
	
}
