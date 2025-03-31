package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;
import static de.uni_halle.barcode_calling.experiments.CallingExperiment.getRepresentativeBarcodeSet;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.indexes.KMerIndexCounter;
import de.uni_halle.barcode_calling.callers.indexes.OrderedIndex;
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
public class OperationCountAnalysis implements Experiment {

	/*
	 * Constants
	 */
	static final RandomDistribution DEFAULT_ABUNDANCE_DISTRIBUTION = RandomDistribution.CONSTANT;
	static final int DEFAULT_MUTATION_RATE_STEPS = 3;
	static final double DEFAULT_MUTATION_RATE = 0.2;
	static final double DEFAULT_MUTATION_RATE_STEP_SIZE = 0.1;
	static final int DEFAULT_READS_COUNT = 10000;
	static final int DEFAULT_BARCODE_COUNT = 1000000;
	static final int DEFAULT_BARCODE_LENGTH = 34;
	static final int DEFAULT_BARCODE_DIST = 0;
	static final String[] CSV_HEADER = {
		"index name", "barcode count", "barcode length", "mutual edit distance", "error model", "mutation probability",
		"operation count", "candidates", "execution time measured", "execution time iterate k-mers total", 
		"execution time sorting", "execution time total"
	};
	static final String DEFAULT_EXPERIMENT_NAME = "artificial-dataset-index-operation-count-analysis";
	
	/*
	 * Attribute declarations
	 */
	BarcodeDataset[] barcodeSets;
	RandomDistribution abundanceDistribution;
	KMerIndexCounter[] indexes;
	String name;
	Path outFilePath;
	
	/*
	 * Constructors
	 */
	public OperationCountAnalysis(
			BarcodeDataset[] barcodeSets,
			KMerIndexCounter[] indexes,
			RandomDistribution abundanceDistribution
	) {
		this.barcodeSets = barcodeSets;
		this.indexes = indexes;
		this.abundanceDistribution = abundanceDistribution;
		this.setExperimentName(DEFAULT_EXPERIMENT_NAME);
	}
	
	public OperationCountAnalysis(
			BarcodeDataset[] barcodeSets,
			KMerIndexCounter[] indexes
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
		
		List<Object[]> csvEntries = new ArrayList<>();
		
		analyzeIndexes(csvEntries);
		
		Main.writeToLog("Experiment ran successfully!");
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
	
	private void analyzeIndexes(List<Object[]> csvEntries) {
		DNAErrorModel errorModel = DNAErrorModel.UPHOFF;
		
		// Set default barcode set
		BarcodeDataset defaultBarcodeSet = getRepresentativeBarcodeSet(this.barcodeSets, DEFAULT_BARCODE_COUNT, DEFAULT_BARCODE_LENGTH, DEFAULT_BARCODE_DIST);
		
		for (KMerIndexCounter index : this.indexes) {
			Main.writeToLog("Analyzing index" + index.getName());
			
			for (BarcodeDataset barcodes : this.barcodeSets) {
			
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
//						this.analyze(index, barcodesNew, errorModel, DEFAULT_MUTATION_RATE, DEFAULT_READS_COUNT, initTime, initMemory, csvEntries);
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
						this.analyze(index, barcodes, errorModel, p, DEFAULT_READS_COUNT, initTime, initMemory, csvEntries);
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
			KMerIndexCounter index, BarcodeDataset barcodes, DNAErrorModel errorModel, double p, int readsCount, 
			double timeElapsedInit, double memoryUsageInit, List<Object[]> csvEntries) {
		
		int barcodeCount = barcodes.getSize();
		int barcodeLength = barcodes.getSequenceLength();
		int editDist = barcodes.getMinMutualEditDist();
		
		// Generate reads
		LabeledDataset reads = barcodes.corrupt(p/3, p/3, p/3, readsCount, errorModel);
		Main.writeToLog("Corrupted dataset with " + readsCount + " reads and p = " + p + " created with error model " + errorModel + "...");
		
		// Analyze
		String indexName = index.getName();
		index.resetCounter();
			
		// ---
		// Measure execution time
		long startTime = System.nanoTime();
		
		index.getRecallAtPos(reads);
		
		long finishTime = System.nanoTime();
		double timeElapsedCalling = (double)(finishTime - startTime)/1000000;  // in ms
		// ---
		
		List<Integer> operationCounts = index.getOperationCounter();
		List<Integer> candidateCounts = index.getCandidateCounter();
		List<Double> timeInsertHitDistance = index.timeInsertHitDistance;
		List<Double> timeSorting = index.timeSorting;
		List<Double> timeTotal = index.timeTotal;
		
		for (int i = 0; i < readsCount; i++) {
			Object[] entry = {
				indexName,
				barcodeCount,
				barcodeLength,
				editDist,
				errorModel,
				p,
				operationCounts.get(i),
				candidateCounts.get(i),
				timeElapsedCalling,
				timeInsertHitDistance.get(i),
				timeSorting.get(i),
				timeTotal.get(i),
			};
			csvEntries.add(entry);
		}
		
		
		Main.writeToLog("> Index \"" + indexName + "\" finished in " + timeElapsedCalling + " ms...");
		
		IOMethods.writeCSV(outFilePath, CSV_HEADER, csvEntries);
		
		System.out.println();
	}
	
	
	/*
	 * Factory methods
	 */
	public static OperationCountAnalysis createIndexAnalysis(
			Path[] barcodesFiles
	) {
		KMerIndexCounter[] indexes = {
				new KMerIndexCounter(4),
				new KMerIndexCounter(5),
				new KMerIndexCounter(6),
				new KMerIndexCounter(7),
		};
		
		BarcodeDataset[] barcodeSets = BarcodeDataset.datasetsFromPaths(barcodesFiles);
		OperationCountAnalysis exp = new OperationCountAnalysis(barcodeSets, indexes);
		
		return exp;
	}
	
}
