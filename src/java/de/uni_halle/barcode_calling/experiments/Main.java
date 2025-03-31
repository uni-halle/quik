package de.uni_halle.barcode_calling.experiments;

import java.io.File;
import java.nio.file.Path;

import static de.uni_halle.barcode_calling.Main.setMaxThreadCount;
import static de.uni_halle.barcode_calling.Main.setLogFile;
import static de.uni_halle.barcode_calling.Main.logHeader;

import de.uni_halle.barcode_calling.callers.indexes.PressCommunicator;
import de.uni_halle.barcode_calling.util.IOMethods;

/**
 * Executable class for running experiments.
 * @author Riko Uphoff
 *
 */
public class Main {
	
	public static void printHelp() {
		StringBuilder sb = new StringBuilder();
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -o --out         | [required] |  Set output file or directory for experiment results. In case of directory, a default file name is chosen.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -b --barcodes    | [required] |  Set source file or directory for barcodes. In case of directory, all files in the directory are interpreted as source files. Required even for experiments without barodes.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -t --threads     |            |  Set the number of threads that should be used for calculations. Defaults to the number of cores + 1.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -l --log         |            |  Set output file or directory for logging. In case of directory, a default file name is chosen. If not provided, logs are only printed on the command line.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -e --experiment  | [required] |  Set the type of experiment you want to run. Possibilities are");
		sb.append("\n                  |            |  > error-models: Compare error model used in Press (2022) to the error model by Uphoff et al. (2024).");
		sb.append("\n                  |            |  > index: Compare hit rate of k-mer indexes with k = 4, ..., 7 to the index proposed by Press (2022) on specified barcode sets. ATTENTION: Requires python listener for Press' code to already be running.");
		sb.append("\n                  |            |  > index-operation-count: Count inner loop iterations of k-mer indexes with k = 4, ..., 7 on specified barcode sets.");
		sb.append("\n                  |            |  > calling: Calculate precision and recall of k-mer calling algorithm with k = 4, ..., 7 and prefiltered variants using Levenshtein distance, aswell as by Press (2022) on specified barcode sets. ATTENTION: Requires python listener for Press' code to already be running.");
		sb.append("\n                  |            |  > calling-SL: Calculate precision and recall of k-mer calling algorithm with k = 4, ..., 7 and prefiltered variants using Sequence-Levenshtein distance, aswell as by Press (2022) on specified barcode sets. ATTENTION: Requires python listener for Press' code to already be running.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -h --help        |            |  Show possible arguments.");
		sb.append("\n-----------------------------------------------------");
		
		System.out.println(sb.toString());
	}
	
	private static Experiment resolveExperiment(String experimentName, Path[] barcodeFiles) {
		Experiment exp = null;
		
		switch (experimentName) {
//			case "calling":
//				exp = ArtificialDatasetIndexCallingExperiment.createHeuristicAnalysisL(barcodeFiles);
//				break;
//			case "calling-SL":
//				exp = ArtificialDatasetIndexCallingExperiment.createHeuristicAnalysisSL(barcodeFiles);
//				break;
			case "pseudo-distance":
				exp = OrderedIndexAnalysis.createIndexAnalysis(barcodeFiles);
				break;
			case "pseudo-distance-operation-count":
				exp = OperationCountAnalysis.createIndexAnalysis(barcodeFiles);
				break;
			case "error-models":
				exp = new ErrorModelComparison();
				break;
//			case "precision-prediction":
//				exp = new PrecisionPredictionExperiment(barcodeFiles[0]);
//				break;
			default:
				throw new RuntimeException("Experiment with tag '" + experimentName + "' does not exist!");
		}
		
		return exp;
	}
	
	public static void main(String[] args) {
		String outFileString = null, barcodesPathString = null, experimentName = null, logFileString = null;
		
		// Parse input arguments
		String arg;
		
		for (int i = 0; i < args.length; i++) {
			arg = args[i];
			
			if (arg.equals("-e") || arg.equals("--experiment")) {
				i++;
				experimentName = args[i];
			}
			else if (arg.equals("-o") || arg.equals("--out")) {
				i++;
				outFileString = args[i];
			}
			else if (arg.equals("-b") || arg.equals("--barcodes")) {
				i++;
				barcodesPathString = args[i];
			}
			else if (arg.equals("-t") || arg.equals("--threads")) {
				i++;
				setMaxThreadCount(Integer.parseInt(args[i]));
			}
			else if (arg.equals("-l") || arg.equals("--log")) {
				i++;
				logFileString = args[i];
			}
			else if (arg.equals("-h") || arg.equals("--help")) {
				printHelp();
				System.exit(0);
			}
			else {
				System.err.println("Warning: unknown parameter" + arg);
				System.err.println("Skipping parameter and hoping for the best...");
				System.err.println("Use -h or --help to display possible arguments.");
			}
		}
		
		// Stop if relevant information not given
		if (barcodesPathString == null) {
			throw new RuntimeException("Source file for barcodes is undefined! Use -h or --help to display possible arguments.");
		}
		else if (outFileString == null) {
			throw new RuntimeException("Output file or directory is undefined! Use -o to set output file or directory. Use -h or --help to display possible arguments.");
		}
		else if (experimentName == null) {
			throw new RuntimeException("Experiment name is undefined! Use -e to set the type of experiment. Use -h or --help to display possible arguments.");
		}
		
		// Get input paths
		Path inputPath = Path.of(barcodesPathString);
		Path[] barcodeFiles = IOMethods.resolveInputPath(inputPath);
		
		// Init experiment
		Experiment exp = resolveExperiment(experimentName, barcodeFiles);
		
		// Set log file
		setLogFile(logFileString, exp.getExperimentName());
		logHeader(args);
		
		// Do experiment
		Path outPath = Path.of(outFileString);
		
		exp.doExperiment(outPath);
		
		System.exit(0);
	}

}
