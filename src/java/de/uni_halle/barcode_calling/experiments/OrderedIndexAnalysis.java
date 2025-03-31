package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;
import static de.uni_halle.barcode_calling.experiments.CallingExperiment.getRepresentativeBarcodeSet;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.indexes.KMerIndex;
import de.uni_halle.barcode_calling.callers.indexes.OrderedIndex;
import de.uni_halle.barcode_calling.callers.indexes.TrimerIndexOutsourced;
import de.uni_halle.barcode_calling.experiments.data.DNAErrorModel;
import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset.RandomDistribution;

/**
 * Implements Experiment 2 from the thesis.
 * @author Riko Uphoff
 *
 */
public class OrderedIndexAnalysis implements Experiment {

	/*
	 * Constants
	 */
	static final RandomDistribution DEFAULT_ABUNDANCE_DISTRIBUTION = RandomDistribution.CONSTANT;
	static final int DEFAULT_MUTATION_RATE_STEPS = 3;
	static final double DEFAULT_MUTATION_RATE = 0.2;
	static final double DEFAULT_MUTATION_RATE_STEP_SIZE = 0.1;
	static final int DEFAULT_READS_COUNT = 10000;
	static final int DEFAULT_MAX_ENTRIES = 2000;
	static final int DEFAULT_BARCODE_COUNT = 1000000;
	static final int DEFAULT_BARCODE_LENGTH = 34;
	static final int DEFAULT_BARCODE_DIST = 0;
	static final String[] CSV_HEADER_OUT = {
		"id", "index position", "recall"
	};
	static final String[] CSV_HEADER_EXTRA = {
		"id", "index name", "barcode count", "barcode length", "mutual edit distance", "error model", "mutation probability", "execution time", "initialization time", "memory usage"
	};
	static final String DEFAULT_EXPERIMENT_NAME = "artificial-dataset-index-analysis";
	
	/*
	 * Attribute declarations
	 */
	BarcodeDataset[] barcodeSets;
	RandomDistribution abundanceDistribution;
	OrderedIndex[] indexes;
	String name;
	Path outFilePath;
	Path extraFilePath;
	int entryCounter = 0;
	
	/*
	 * Constructors
	 */
	public OrderedIndexAnalysis(
			BarcodeDataset[] barcodeSets,
			OrderedIndex[] indexes,
			RandomDistribution abundanceDistribution
	) {
		this.barcodeSets = barcodeSets;
		this.indexes = indexes;
		this.abundanceDistribution = abundanceDistribution;
		this.setExperimentName(DEFAULT_EXPERIMENT_NAME);
	}
	
	public OrderedIndexAnalysis(
			BarcodeDataset[] barcodeSets,
			OrderedIndex[] indexes
	) {
		this(barcodeSets, indexes, DEFAULT_ABUNDANCE_DISTRIBUTION);
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
		
		this.outFilePath = this.resolveOutPath(outPath);
		this.extraFilePath = this.resolveExtraPath(outPath);
		
		List<Object[]> csvEntriesOut = new ArrayList<>();
		List<Object[]> csvEntriesExtra = new ArrayList<>();
		
		analyzeIndexes(csvEntriesOut, csvEntriesExtra);
		
		Main.writeToLog("Experiment ran successfully!");
	}
	
	public Path resolveExtraPath(Path outPath) {
		String resolvedOutPath = this.resolveOutPath(outPath).toString();
		resolvedOutPath = resolvedOutPath.substring(0, resolvedOutPath.length() - 4) + "_extra.csv";
		
		return Path.of(resolvedOutPath);
	}
	
	@Override
	public String getDefaultResultFileName() {
		StringBuilder sb = new StringBuilder(this.getExperimentName());
		sb.append("_r");
		sb.append(DEFAULT_READS_COUNT);
		sb.append("_");
		sb.append(getTimestampString());
		sb.append(".csv");
		
		return sb.toString();
	}
	
	private void analyzeIndexes(List<Object[]> csvEntriesOut, List<Object[]> csvEntriesExtra) {
		
		// Set default barcode set
		BarcodeDataset defaultBarcodeSet = getRepresentativeBarcodeSet(this.barcodeSets, DEFAULT_BARCODE_COUNT, DEFAULT_BARCODE_LENGTH, DEFAULT_BARCODE_DIST);
		DNAErrorModel errorModel = DNAErrorModel.UPHOFF;	
		
		for (OrderedIndex index : this.indexes) {
			for (BarcodeDataset barcodes : this.barcodeSets) {
				Main.writeToLog("Analyzing index" + index.getName());
				
				// Iterate over barcode count
//				if (barcodes == defaultBarcodeSet) {
//	
//					double log2 = Math.log(barcodes.getSize()) / Math.log(2);
//					
//					for (int i = 10; i < log2; i++) {
//						int barcodeCountNew = (int)Math.pow(2, i);
//						BarcodeDataset barcodesNew = barcodes.getSubset(barcodeCountNew);
//						
//						double[] initMeasures = initIndex(index, barcodesNew);
//						double initTime = initMeasures[0];
//						double initMemory = initMeasures[1];
//						
//						this.analyze(index, barcodesNew, errorModel, DEFAULT_MUTATION_RATE, DEFAULT_READS_COUNT, initTime, initMemory, csvEntriesOut, csvEntriesExtra);
//					}
//				}
				
				double[] initMeasures = initIndex(index, barcodes);
				double initTime = initMeasures[0];
				double initMemory = initMeasures[1];
				
				// Iterate over mutation probability
				double stepSize = DEFAULT_MUTATION_RATE_STEP_SIZE;
				double p;
				
				for (int i = 1; i <= DEFAULT_MUTATION_RATE_STEPS; i++) {
					p = i*stepSize;
					
					if (this.barcodeSets.length <= 3 || barcodes == defaultBarcodeSet || p == DEFAULT_MUTATION_RATE) {
						this.analyze(index, barcodes, errorModel, p, DEFAULT_READS_COUNT, initTime, initMemory, csvEntriesOut, csvEntriesExtra);
					}
				}
			}
			
			// Free some cache
			index.setBarcodes(new BarcodeDataset());
		}
	}
	
	private static double[] initIndex(OrderedIndex index, BarcodeDataset barcodes) {
		Main.writeToLog("Analyzing index for n = " + barcodes.getSize() + ", l = " + barcodes.getSequenceLength() + ", d = " + barcodes.getMinMutualEditDist() + "...");
		
		index.setBarcodesFractionProbed(1.0);
		
		// Measure init time
		long startTime = System.nanoTime();
		
		index.setBarcodes(barcodes);
		
		long finishTime = System.nanoTime();
		long finishMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		
		double timeElapsedInit = (double)(finishTime - startTime)/1000000;	// ms
		double memoryUsageInit = (double)(finishMemory)/1024;	// kB
		
		Main.writeToLog("> Index \"" + index.getName() + "\" initialized in " + timeElapsedInit + " ms with memory peak " + memoryUsageInit + " kB...");
		
		double[] out = {timeElapsedInit, memoryUsageInit};
		
		return out;
	}
	
	private void analyze(
			OrderedIndex index, BarcodeDataset barcodes, DNAErrorModel errorModel, double p, int readsCount, 
			double timeElapsedInit, double memoryUsageInit, List<Object[]> csvEntriesOut, List<Object[]> csvEntriesExtra) {
		
		int barcodeCount = barcodes.getSize();
		int barcodeLength = barcodes.getSequenceLength();
		int editDist = barcodes.getMinMutualEditDist();
		
		// Generate reads
		LabeledDataset reads = barcodes.corrupt(p/3, p/3, p/3, readsCount, errorModel);
		Main.writeToLog("Corrupted dataset with " + readsCount + " reads and p = " + p + " created with error model " + errorModel + "...");
		
		// Analyze
		String indexName = index.getName();
			
		// ---
		// Measure execution time
		long startTime = System.nanoTime();
		
		double[] recallAtPos = index.getRecallAtPos(reads);
		
		long finishTime = System.nanoTime();
		double timeElapsedCalling = (double)(finishTime - startTime)/1000000;  // in ms
		// ---
		
		int pos = 0;
		double pseudoPos;
		
		// Calculate recall
		for (int i = 0; i <= DEFAULT_MAX_ENTRIES; i++) {
			if (i < 10) {
				pos = i+1;
			}
			else {
				// Equal spacing on a log scale
				pseudoPos = (Math.pow(10, (double)(i-10) / (double)(DEFAULT_MAX_ENTRIES-10)) - 1) / 9 * barcodeCount;
				pos = Math.max(pos+1, (int)pseudoPos);
			}
			
			if (pos < barcodeCount) {
				double recall = recallAtPos[pos];
				
				Object[] entry = {
					this.entryCounter,
					pos,
					recall,
				};
				csvEntriesOut.add(entry);
			}
		}
		
		Object[] entry = {
			this.entryCounter,
			indexName,
			barcodeCount,
			barcodeLength,
			editDist,
			errorModel,
			p,
			timeElapsedCalling / (double)readsCount,	// Millisecond per read
			timeElapsedInit,
			memoryUsageInit,
		};
		csvEntriesExtra.add(entry);
		
		this.entryCounter++;
		
		Main.writeToLog("> Index \"" + indexName + "\" finished in " + timeElapsedCalling + " ms...");
		
		IOMethods.writeCSV(outFilePath, CSV_HEADER_OUT, csvEntriesOut);
		IOMethods.writeCSV(extraFilePath, CSV_HEADER_EXTRA, csvEntriesExtra);
		
		System.out.println();
	}
	
	
	/*
	 * Factory methods
	 */
	public static OrderedIndexAnalysis createIndexAnalysis(
			Path[] barcodesFiles
	) {
		OrderedIndex[] indexes = {
				new KMerIndex(4),
				new KMerIndex(5),
				new KMerIndex(6),
				new KMerIndex(7),
				new TrimerIndexOutsourced(),
		};
		
		BarcodeDataset[] barcodeSets = BarcodeDataset.datasetsFromPaths(barcodesFiles);
		OrderedIndexAnalysis exp = new OrderedIndexAnalysis(barcodeSets, indexes);
		
		return exp;
	}
	
}
