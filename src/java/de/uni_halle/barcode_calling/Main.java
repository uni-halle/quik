package de.uni_halle.barcode_calling;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.BiFunction;

import de.uni_halle.barcode_calling.callers.BarcodeCaller;
import de.uni_halle.barcode_calling.callers.FilteredKMerCaller;
import de.uni_halle.barcode_calling.callers.IndexCaller;
import de.uni_halle.barcode_calling.callers.KMerCaller;
import de.uni_halle.barcode_calling.callers.MaxDistLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.MaxDistSequenceLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.NaiveLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.NaiveSequenceLevenshteinCaller;
import de.uni_halle.barcode_calling.callers.TrimerCaller;
import de.uni_halle.barcode_calling.experiments.Experiment;
import de.uni_halle.barcode_calling.util.IOMethods;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.ReadDataset;

/**
 * Executable class for calling a set of reads against a barcode library.
 * @author Riko Uphoff
 *
 */
public class Main {
	
	/**
	 * Max. number of threads to be run concurrently in this application.
	 */
	public static int maxThreadCount = Runtime.getRuntime().availableProcessors();// + 1;
	/**
	 * Global singleton thread-pool for bound multi-threading.
	 */
	public static ExecutorService threadPool = Executors.newFixedThreadPool(maxThreadCount);
	
	private static File logFile = null;
	
	/**
	 * Sets the max. number of threads to be run concurrently in this application.
	 * @param newMaxThreadCount
	 */
	public static void setMaxThreadCount(int newMaxThreadCount) {
		maxThreadCount = newMaxThreadCount;
		threadPool = Executors.newFixedThreadPool(maxThreadCount);
	}
	
	public static void setLogFile(String logFileString, String baseName) {
		
		if (logFileString != null) {
			Path logFilePath = Path.of(logFileString);
			logFile = logFilePath.toFile();
			
			if (logFile.isDirectory()) {
				logFilePath = Path.of(logFileString, getDefaultLogFileName(baseName));
				logFile = logFilePath.toFile();
			}
		}
		else {
			logFile = null;
		}
		
	}
	
	public static void setLogFile(String logFileString) {
		setLogFile(logFileString, null);
	}
	
	public static void writeToLog(String logText) {
		// TODO Library for logging
		System.out.println("Logging:");
		System.out.println(logText);
		
		if (logFile != null) {
			try (FileWriter writer = new FileWriter(logFile, true)) {
				writer.write(String.format("\n[%s] ", getTimestampString()));
				writer.write(logText);
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			catch (RuntimeException e) {
				e.printStackTrace();
			}
		}
		
	}
	
	public static void writeToLog(Exception e) {
		// TODO Library for logging
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		e.printStackTrace(pw);
		
		writeToLog(sw.toString());
	}
	
	public static String getDefaultFileName(String baseName, String purpose, String outExtension) {
		
		StringBuilder sb = new StringBuilder();
		
		if (baseName != null) {
			int extensionSeperatorIndex = baseName.lastIndexOf('.');
			if (extensionSeperatorIndex >= 0) {
				String extension = baseName.substring(extensionSeperatorIndex);
				baseName = baseName.substring(0, extensionSeperatorIndex);
			}
			
			sb.append(baseName);
			sb.append("_");
		}
		
		sb.append(purpose);
		sb.append("_");
		sb.append(getTimestampString());
		sb.append(".");
		sb.append(outExtension);
		
		return sb.toString();
	}
	
	public static String getDefaultLogFileName(String baseName) {
		return getDefaultFileName(baseName, "log", "txt");
	}
	
	public static String getDefaultResultFileName(String baseName) {
		return getDefaultFileName(baseName, "result", "csv");
	}
	
	public static void printHelp() {
		StringBuilder sb = new StringBuilder();
		sb.append("\n----------------------------------------------------------");
		sb.append("\n -o --out           | [required] |  Set output file or directory for matching results. In case of directory, a default file name is chosen.");
		sb.append("\n----------------------------------------------------------");
		sb.append("\n -b --barcodes      | [required] |  Set source file for barcodes. Must be in .txt format, one sequence per line.");
		sb.append("\n----------------------------------------------------------");
		sb.append("\n -r --reads         | [required] |  Set source file for reads. Must be in .txt format, one sequence per line.");
		sb.append("\n----------------------------------------------------------");
		sb.append("\n -a --algorithm     | [required] |  Set calling algorithm. Must be one of:\n" + Main.getAvailableCallers().keySet());
		sb.append("\n----------------------------------------------------------");
		sb.append("\n -t --threads       |            |  Set the number of threads that should be used for calculations. Defaults to the number of cores + 1.");
		sb.append("\n----------------------------------------------------------");
		sb.append("\n -l --log           |            |  Set output file or directory for logging. In case of directory, a default file name is chosen. If not provided, logs are only printed on the command line.");
		sb.append("\n----------------------------------------------------------");
//		sb.append("\n -prec --precision  |            |  Set the minimum desired precision for calling results. Defaults to 0.");
//		sb.append("\n----------------------------------------------------------");
//		sb.append("\n -rec --recall      |            |  Set the minimum desired recall for calling results. Defaults to 0.");
//		sb.append("\n----------------------------------------------------------");
		sb.append("\n -h --help          |            |  Show possible arguments.");
		sb.append("\n----------------------------------------------------------");
		
		System.out.println(sb.toString());
	}
	
	public static void writeCallingResults(
			Path outFilePath, 
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond,
			ReadDataset reads, 
			BarcodeDataset barcodes
		) {
		
		List<Object[]> csvEntries = new ArrayList<>();
		
		String[] header = {"read_tag", "barcode_tag", "distance", "barcode_tag_second", "distance_second"};
		
		for (int i = 0; i < callingResults.length; i++) {
			Object[] entry = {
					reads.getTag(i), 
					barcodes.getTag(callingResults[i]), 
					callingDistances[i],
					barcodes.getTag(callingResultsSecond[i]), 
					callingDistancesSecond[i]
			};
			csvEntries.add(entry);
		}
		
		IOMethods.writeCSV(outFilePath, header, csvEntries);
		
		writeToLog("Results written to " + outFilePath);
	}
	
	public static void logHeader(String[] args) {
		StringBuilder sb = new StringBuilder();
		sb.append("Command line arguments used:\n");
		sb.append(String.join(" ", args));
		sb.append("\nMaximum number of threads: ");
		sb.append(maxThreadCount);
		
		writeToLog(sb.toString());
	}
	
	private static void logRecall(int[] callingResults) {
		// Log actual recall
		double hits = 0;
		for (int i = 0; i < callingResults.length; i++) {
			if (callingResults[i] >= 0)
				hits++;
		}
		
		double recall = hits / (double)callingResults.length;
		writeToLog("Recall recieved in calling: " + recall);
	}
	
	public static Map<String, BiFunction<Map<String, Object>, BarcodeDataset, BarcodeCaller>> getAvailableCallers() {
		Map<String, BiFunction<Map<String, Object>, BarcodeDataset, BarcodeCaller>> availableCallers = new HashMap<>();
		availableCallers.put("naive-levenshtein", NaiveLevenshteinCaller::fromParameters);
		availableCallers.put("naive-sequence-levenshtein", NaiveSequenceLevenshteinCaller::fromParameters);
		availableCallers.put("max-distance-levenshtein", MaxDistLevenshteinCaller::fromParameters);
		availableCallers.put("max-distance-sequence-levenshtein", MaxDistSequenceLevenshteinCaller::fromParameters);
		availableCallers.put("kmer", KMerCaller::fromParameters);
		availableCallers.put("filtered-kmer", FilteredKMerCaller::fromParameters);
//		availableCallers.put("press2022", new TrimerCaller());
		
		return availableCallers;
	}

	
	public static void main(String[] args) {
		String outFileString = null, logFileString = null, barcodesPathString = null, readsPathString = null, algorithmTag = null;
		double minPrec = 0;
		double minRec = 0;
		Map<String, BiFunction<Map<String, Object>, BarcodeDataset, BarcodeCaller>> availableCallers = Main.getAvailableCallers();
		Map<String, Object> algorithmParameters = new HashMap<>();
		
		// Parse input arguments
		// TODO Use library
		String arg;
		
		for (int i = 0; i < args.length; i++) {
			arg = args[i];
			
			if (arg.equals("-prec") || arg.equals("--precision")) {
				// TODO Not supported anymore
				i++;
				minPrec = Double.parseDouble(args[i]);
			}
			else if (arg.equals("-rec") || arg.equals("--recall")) {
				// TODO Not supported anymore
				i++;
				minRec = Double.parseDouble(args[i]);
			}
			else if (arg.equals("-t") || arg.equals("--threads")) {
				i++;
				setMaxThreadCount(Integer.parseInt(args[i]));
			}
			else if (arg.equals("-o") || arg.equals("--out")) {
				i++;
				outFileString = args[i];
			}
			else if (arg.equals("-l") || arg.equals("--log")) {
				i++;
				logFileString = args[i];
			}
			else if (arg.equals("-r") || arg.equals("--reads")) {
				i++;
				readsPathString = args[i];
			}
			else if (arg.equals("-b") || arg.equals("--barcodes")) {
				i++;
				barcodesPathString = args[i];
			}
			else if (arg.equals("-a") || arg.equals("--algorithm")) {
				i++;
				algorithmTag = args[i];
			}
			else if (arg.equals("-h") || arg.equals("--help")) {
				printHelp();
				System.exit(0);
			}
			else if (arg.startsWith("-p")) {
				// Parameter for algorithm
				// TODO Handle non-integer parameters
				arg = arg.substring(2);
				i++;
				algorithmParameters.put(arg, Integer.parseInt(args[i]));
			}
			else {
				System.err.println("Warning: unknown parameter" + arg);
				System.err.println("Skipping parameter and hoping for the best...");
			}
		}
		
		// Stop if relevant information not given
		if (barcodesPathString == null) {
			throw new RuntimeException("Source file for barcodes is undefined! Use -b to set input file or -h to view instructions.");
		}
		else if (Path.of(barcodesPathString).toFile().isDirectory()) {
			throw new RuntimeException("Source file for barcodes is a directory!");
		}
		if (readsPathString == null) {
			throw new RuntimeException("Source file for reads is undefined! Use -r to set input file or -h to view instructions.");
		}
		else if (Path.of(readsPathString).toFile().isDirectory()) {
			throw new RuntimeException("Source file for reads is a directory!");
		}
		else if (outFileString == null) {
			throw new RuntimeException("Output file or directory is undefined! Use -o to set output file or directory or -h to view instructions.");
		}
		else if (algorithmTag == null) {
			throw new RuntimeException("No calling algorithm specified! Use -a to set the algorithm or -h to view instructions.");
		}
		else if (!availableCallers.containsKey(algorithmTag)) {
			throw new RuntimeException("Unknown calling algorithm! Must be one of:\n" + availableCallers.keySet());
		}
		// TODO Warning when log file isn't set
		
		// Get input paths
		Path barcodeFilePath = Path.of(barcodesPathString);
		Path readsFilePath = Path.of(readsPathString);
		
		// Get output paths
		String barcodesFileName = barcodeFilePath.getFileName().toString();
		setLogFile(logFileString, barcodesFileName);
		//  - Result file
		Path outFilePath = Path.of(outFileString);
		
		if (outFilePath.toFile().isDirectory()) {
			outFilePath = Path.of(outFileString, getDefaultResultFileName(barcodesFileName));
		}
		
		// -- Calling --//
		long startTime, finishTime;
		
		logHeader(args);
		
		// Init calling
		try {
			
			BarcodeDataset barcodes = BarcodeDataset.fromFile(barcodeFilePath);
			writeToLog(barcodes.getSize() + " barcodes read from " + barcodesPathString + " with length " + barcodes.getSequenceLength() + "  and mutual distance of " + barcodes.getMinMutualEditDist());
			
			ReadDataset reads = ReadDataset.fromFile(readsFilePath, barcodes.getSequenceLength());
			writeToLog(reads.getSize() + " reads read from " + readsPathString);
			
			BarcodeCaller callingAlgorithm = availableCallers.get(algorithmTag).apply(
					algorithmParameters,
					barcodes
			);
			
			// Init calling algorithm
			startTime = System.nanoTime();
			
			writeToLog("Calling with " + callingAlgorithm);
			
			finishTime = System.nanoTime();
			double timeElapsedInit = (double)(finishTime - startTime)/1000000000;
			
			writeToLog("Initialization finished in " + timeElapsedInit + " s.");
			
			// Tune parameters of calling algorithm and measure time
//			startTime = System.nanoTime();
//			
//			callingAlgorithm.tuneParameters(reads, minPrec, minRec);
//			
//			finishTime = System.nanoTime();
//			double timeElapsedTuning = (double)(finishTime - startTime)/1000000000;
//			
//			writeToLog("Finished tuning in " + timeElapsedTuning + " s.");
			
			// Identify reads and measure time
			int[] callingResults = new int[reads.getSize()];
			int[] callingDistances = new int[reads.getSize()];
			int[] callingResultsSecond = new int[reads.getSize()];
			int[] callingDistancesSecond = new int[reads.getSize()];
			
			for (int i = 0; i < reads.getSize(); i++) {
				callingResults[i] = -1;
				callingResultsSecond[i] = -1;
			}
			
			startTime = System.nanoTime();
			
			callingAlgorithm.callParallel(
					reads.getSequences(), 
					callingResults, 
					callingDistances,
					callingResultsSecond, 
					callingDistancesSecond
				);
			
			finishTime = System.nanoTime();
			double timeElapsedMatching = (double)(finishTime - startTime)/1000000000;
			
			writeToLog("Finished calling in " + timeElapsedMatching + " s.");
			
			logRecall(callingResults);
			
			writeCallingResults(
					outFilePath, 
					callingResults, 
					callingDistances,
					callingResultsSecond, 
					callingDistancesSecond, 
					reads, 
					barcodes
				);
			
		} catch (Exception e) {
			
			Main.writeToLog(e);
			
		}
		
		System.exit(0);
	}

}
