package de.uni_halle.barcode_calling.experiments.data;

import static de.uni_halle.barcode_calling.Main.logHeader;
import static de.uni_halle.barcode_calling.Main.setLogFile;

import java.nio.file.Path;

import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

public class Main {
	
	public static void printHelp() {
		StringBuilder sb = new StringBuilder();
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -r --reads       | [required] |  Set output file for generated reads.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -l --labels      | [required] |  Set output file for the labels of the generated reads.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -b --barcodes    | [required] |  Set source file or directory for barcodes.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -n --number      | [required] |  Set the number of reads to by generated.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -s --pSub        | [required] |  Set the probability per base for substitution.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -i --pIns        | [required] |  Set the probability per base for insertion.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -d --pDel        | [required] |  Set the probability per base for deletion.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -l --log         |            |  Set output file or directory for logging. In case of directory, a default file name is chosen. If not provided, logs are only printed on the command line.");
		sb.append("\n-----------------------------------------------------");
		sb.append("\n -h --help        |            |  Show possible arguments.");
		sb.append("\n-----------------------------------------------------");
		
		System.out.println(sb.toString());
	}
	
	public static void main(String[] args) {
		String readsFileString = null, labelsFileString = null, barcodesPathString = null, logFileString = null;
		int n = -1;
		double pSub = -1, pIns = -1, pDel = -1;
		
		// Parse input arguments
		String arg;
		
		for (int i = 0; i < args.length; i++) {
			arg = args[i];
			
			if (arg.equals("-r") || arg.equals("--reads")) {
				i++;
				readsFileString = args[i];
			}
			else if (arg.equals("-l") || arg.equals("--labels")) {
				i++;
				labelsFileString = args[i];
			}
			else if (arg.equals("-b") || arg.equals("--barcodes")) {
				i++;
				barcodesPathString = args[i];
			}
			else if (arg.equals("-n") || arg.equals("--number")) {
				i++;
				n = Integer.parseInt(args[i]);
			}
			else if (arg.equals("-s") || arg.equals("--pSub")) {
				i++;
				pSub = Double.parseDouble(args[i]);
			}
			else if (arg.equals("-i") || arg.equals("--pIns")) {
				i++;
				pIns = Double.parseDouble(args[i]);
			}
			else if (arg.equals("-d") || arg.equals("--pDel")) {
				i++;
				pDel = Double.parseDouble(args[i]);
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
		else if (readsFileString == null) {
			throw new RuntimeException("Output file for reads is undefined! Use -r to set output file. Use -h or --help to display possible arguments.");
		}
		else if (labelsFileString == null) {
			throw new RuntimeException("Output file for labels is undefined! Use -l to set output file. Use -h or --help to display possible arguments.");
		}
		else if (n < 0) {
			throw new RuntimeException("Number of reads not specified or invalid. Use -n to set reads count. Use -h or --help to display possible arguments.");
		}
		else if (pSub < 0) {
			throw new RuntimeException("Substitution probability not specified or invalid. Use -s to set substitution probability. Use -h or --help to display possible arguments.");
		}
		else if (pIns < 0) {
			throw new RuntimeException("Insertion probability not specified or invalid. Use -i to set insertion probability. Use -h or --help to display possible arguments.");
		}
		else if (pDel < 0) {
			throw new RuntimeException("Deletion probability not specified or invalid. Use -d to set deletion probability. Use -h or --help to display possible arguments.");
		}
		
		// Get input paths
		Path barcodePath = Path.of(barcodesPathString);
		Path readsPath = Path.of(readsFileString);
		Path labelsPath = Path.of(labelsFileString);
		
		// Set log file
		setLogFile(logFileString);
		logHeader(args);
		
		// Generate reads
		BarcodeDataset barcodes = BarcodeDataset.fromFile(barcodePath);
		LabeledDataset reads = barcodes.corrupt(pSub, pIns, pDel, n);
		reads.toFile(readsPath, labelsPath);
		
		
		System.exit(0);
	}

}

