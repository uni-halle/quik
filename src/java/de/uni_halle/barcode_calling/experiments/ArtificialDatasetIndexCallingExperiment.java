package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.DistanceMeasureBase;
import de.uni_halle.barcode_calling.callers.FilteredKMerCaller;
import de.uni_halle.barcode_calling.callers.IndexCaller;
import de.uni_halle.barcode_calling.callers.KMerCaller;
import de.uni_halle.barcode_calling.callers.LevenshteinBase;
import de.uni_halle.barcode_calling.callers.MaxDistanceCaller;
import de.uni_halle.barcode_calling.callers.SequenceLevenshteinBase;
import de.uni_halle.barcode_calling.callers.TrimerCaller;
import de.uni_halle.barcode_calling.experiments.data.DNAErrorModel;
import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

/**
 * Analyse the impact of the model's parameters
 * @author Riko Uphoff
 *
 */
public class ArtificialDatasetIndexCallingExperiment extends CallingExperiment {
	/*
	 * Constants
	 */
	static final double[] DEFAULT_MUTATION_RATES = {0.1, 0.2, 0.3};
	static final int DEFAULT_READS_COUNT = 10000; // TODO 100_000
	static final double DEFAULT_MIN_PRECISION = 0.99;
	static final double DEFAULT_MIN_RECALL = 0.8;
	static final String[] CSV_HEADER = {
		"caller", "barcode count", "barcode length", "mutual edit distance", "error model", "reads count", "mutation probability", "distance measure",
		"execution time index", "execution time distance", "execution time calling", "precision", "recall", "f1", "L", "pos"
	};
	static final String DEFAULT_EXPERIMENT_NAME = "artificial-dataset-parameter-analysis";
	
	protected BarcodeDataset[] barcodeSets;
	protected double[] mutationRates;
	protected int readsCount;
	protected IndexCaller[] candidates;
	private Path outFilePath;
	private DistanceMeasureBase distanceMeasure;
	
	/*
	 * Constructors
	 */
	public ArtificialDatasetIndexCallingExperiment(
			BarcodeDataset[] barcodeSets,
			IndexCaller[] candidates,
			DistanceMeasureBase distanceMeasure,
			int readsCount,
			double[] mutationRates
	) {
		this.barcodeSets = barcodeSets;
		this.candidates = candidates;
		this.mutationRates = mutationRates;
		this.readsCount = readsCount;
		this.distanceMeasure = distanceMeasure;
		this.setExperimentName(DEFAULT_EXPERIMENT_NAME);
	}
	
	public ArtificialDatasetIndexCallingExperiment(
			BarcodeDataset[] barcodeSets,
			IndexCaller[] candidates,
			DistanceMeasureBase distanceMeasure
	) {
		this(
				barcodeSets, 
				candidates, 
				distanceMeasure,
				DEFAULT_READS_COUNT,
				DEFAULT_MUTATION_RATES
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
		
		this.outFilePath = this.resolveOutPath(outPath);
		
		List<Object[]> csvEntries = new ArrayList<>();
		
		// Iterate over barcode-specific parameters
		for (BarcodeDataset barcodes : this.barcodeSets) {
			for (IndexCaller candidate : this.candidates) {
				candidate.setBarcodes(barcodes);
				
				this.doCalling(barcodes, candidate, this.readsCount, csvEntries);
				
				candidate.setBarcodes(new BarcodeDataset());
			}
		}
		
		IOMethods.writeCSV(outFilePath, CSV_HEADER, csvEntries);
		Main.writeToLog("Experiment ran successfully!");
	}
	
	protected void doCalling(BarcodeDataset barcodes, IndexCaller candidate, int readsCount, List<Object[]> csvEntries) {
		
		int barcodeCount = barcodes.getSize();
		int barcodeLength = barcodes.getSequenceLength();
		int barcodeMutualEditDist = barcodes.getMinMutualEditDist();
		
		// Call
		Main.writeToLog("Calling with n = " + barcodes.getSize() + ", l = " + barcodes.getSequenceLength() + ", d = " + barcodes.getMinMutualEditDist() + "...");
		
		for (double p : this.mutationRates) {
			// Generate reads
			DNAErrorModel errorModel = DNAErrorModel.UPHOFF;
			
			LabeledDataset reads = barcodes.corrupt(p/3, p/3, p/3, readsCount, errorModel);
			Main.writeToLog("Corrupted dataset with " + readsCount + " reads and p = " + p + " created...");
			
			int maxDistOpt = 8;  // candidate.getMaxDist();
			double posOpt = 100 / (double)barcodeCount;  // candidate.getBarcodesFractionProbed();
			
			// Set pos
			candidate.setBarcodesFractionProbed(posOpt);
			
			// Iterate over maxDist
			if(candidate instanceof MaxDistanceCaller) {
				MaxDistanceCaller candidateConv = (MaxDistanceCaller)candidate;
				
				for(int maxDist = 0; maxDist <= barcodeLength / 2; maxDist++) {
					candidateConv.setMaxDist(maxDist);
					
					// Call
					call(candidate, reads, barcodeCount, barcodeLength, barcodeMutualEditDist, p, errorModel, csvEntries);
				}
				
				// Reset
				candidateConv.setMaxDist(maxDistOpt);
			}
			
			// Iterate over pos
			if (candidate instanceof IndexCaller) {
				IndexCaller candidateConv = (IndexCaller)candidate;
				
				int iMax = 17;	// Up to pos ~= 120k
				
				for(int i = 0; i <= iMax; i++) {
					int pos = (int)Math.pow(2, i);
					
					if (pos < barcodeCount) {
						double barcodeFractionProbed = (double)pos / (double)barcodeCount;
						
						candidateConv.getIndex().setBarcodesFractionProbed(barcodeFractionProbed);
						
						// Call
						call(candidate, reads, barcodeCount, barcodeLength, barcodeMutualEditDist, p, errorModel, csvEntries);
					}
				}
				
				// Reset
				candidateConv.getIndex().setBarcodesFractionProbed(posOpt);
			}
		}
		
		System.out.println();
	}
	
	protected void call(
			IndexCaller candidate, LabeledDataset reads, 
			int barcodeCount, int barcodeLength, int barcodeMutualEditDist, double p, 
			DNAErrorModel errorModel, List<Object[]> csvEntries) {
		int maxDist = candidate.getMaxDist();
		double pos = candidate.getBarcodesFractionProbed();
		
		// Call
		long startTime = System.nanoTime();
		int[] calles = candidate.callParallel(reads.getSequences());
		long finishTime = System.nanoTime();

		double timeIndex = (double)(candidate.getExecTimeIndex())/1000000;	// Milliseconds
		double timeNN = (double)(candidate.getExecTime())/1000000;		// Milliseconds
		double timeTotal = timeIndex + timeNN;
		
		Main.writeToLog("> Caller \"" + candidate.getName() + "\" finished calling in " + timeTotal + " ms with L = " + maxDist + " and pos = " + pos + "...");
		Main.writeToLog("TEST: exec time = " + (double)(finishTime-startTime)/1000000);
		
		// Statistics
		double[] callingRates = getMatchingRates(reads.getLabels(), calles);
		
		Object[] entry = {
			candidate.getName(),
			barcodeCount,
			barcodeLength,
			barcodeMutualEditDist,
			errorModel,
			readsCount,
			p,
			distanceMeasure.getName(),
			timeIndex / (double)readsCount,	// Millisecond per read
			timeNN / (double)readsCount,	// Millisecond per read
			timeTotal / (double)readsCount,	// Millisecond per read
			callingRates[0],
			callingRates[1],
			callingRates[2],
			maxDist,
			pos,
		};
		
		csvEntries.add(entry);
		IOMethods.writeCSV(outFilePath, CSV_HEADER, csvEntries);  // Write intermediate results
		
		System.gc();	// Force garbage collection
	}
	
	
	/*
	 * Factory methods
	 */
	public static ArtificialDatasetIndexCallingExperiment createHeuristicAnalysis(
			Path[] barcodesFiles,
			DistanceMeasureBase distanceMeasure
	) {
		IndexCaller[] candidates = {
				new FilteredKMerCaller(4, 7),
				new FilteredKMerCaller(5, 7),
				new KMerCaller(4),
				new KMerCaller(5),
				new KMerCaller(6),
				new KMerCaller(7),
				new TrimerCaller(),
		};
		
		// TODO could be prettier
		for (IndexCaller candidate : candidates) {
			Main.writeToLog("Distance measure: " + distanceMeasure.getName());
			candidate.setDistanceModel(distanceMeasure);
		}
		
		
		BarcodeDataset[] barcodeSets = BarcodeDataset.datasetsFromPaths(barcodesFiles);
		ArtificialDatasetIndexCallingExperiment exp = new ArtificialDatasetIndexCallingExperiment(
				barcodeSets, 
				candidates,
				distanceMeasure
		);
		
		return exp;
	}
	
	public static ArtificialDatasetIndexCallingExperiment createHeuristicAnalysisL(
			Path[] barcodesFiles
	) {
		return createHeuristicAnalysis(barcodesFiles, new LevenshteinBase());
	}
	
	public static ArtificialDatasetIndexCallingExperiment createHeuristicAnalysisSL(
			Path[] barcodesFiles
	) {
		return createHeuristicAnalysis(barcodesFiles, new SequenceLevenshteinBase());
	}
}
