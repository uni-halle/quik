package de.uni_halle.barcode_calling.callers.indexes;

import py4j.GatewayServer;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

public interface PressCommunicator {
	
	/* 
	 * Global variables
	 */
	public static final String PYTHON_SRC_DIR = "src/press_2022";  // TODO Make dynamic?
	public static final String PYTHON_COMMAND = "venv/bin/python3.9";
	public static final String RAM_FILESYSTEM_PATH = "/media/ramdisk";

	
	/*
	 * Interface methods
	 */
	
	public void setBarcodes(BarcodeDataset barcodes);
	public void setBarcodesFile(String barcodesFile);
	public void setValidBarcodes(List<Integer> ids);

	public List<Integer> calcCosineIndex(DNASequence read, int k);
	public List<Integer> calcCosineIndex(DNASequence read, int k, List<Integer> validBarcodes);
	
	public List<Integer> calcHammingIndex(DNASequence read);
	public List<Integer> calcHammingIndex(DNASequence read, List<Integer> validBarcodes);
	
	public double[] getRecallAtPos(String readsFilename, String labelsFilename);
	public void writeFrequencyAtPos(String readsFilename, String labelsFilename, String outFilename);
	
	public void writeTriageCandidates(String readsFilename, String outFilename);
	public void writeTriageCandidates(String readsFilename, String outFilename, String valid_barcodesFilename);
	
	public void decodeReads(String readsFilename, String decodesFilename);
	public void decodeReads(String readsFilename, String decodesFilename, String execTimesFilename);
	public void decodeReads(String readsFilename, String decodesFilename, String execTimesFilename, String distsFilename);
	public void decodeReads(String readsFilename, String decodesFilename, String execTimesFilename, String distsFilename, int pos, int L);
	
	/*
	 * Static methods
	 */
	
	public static PressCommunicator createCommunicator() throws IOException {
		// Start python listener
		/*
		// FIXME Python seems to fail finding the JVM if it is already running...
		if (ConnectionManager.pythonProcess == null || !ConnectionManager.pythonProcess.isAlive()) {
			Path pythonInterpreterPath = Path.of(PressCommunicator.PYTHON_SRC_DIR, PressCommunicator.PYTHON_COMMAND);
			Path pythonScriptPath = Path.of(PressCommunicator.PYTHON_SRC_DIR, "java_communicator.py");
			
			ProcessBuilder processBuilder = new ProcessBuilder(
					pythonInterpreterPath.toString(), pythonScriptPath.toString()
			);
			
			processBuilder.inheritIO();  // Connect IO of subprocess to java program

			System.out.println("Starting python listener...");	// TODO
			ConnectionManager.pythonProcess = processBuilder.start();
		    System.out.println("Python listener started!");	// TODO
		}
		*/
		
		// Connect python listener to gateway server
		// GatewayServer.turnLoggingOff();  TODO
		System.out.println("Starting Java gateway...");	// TODO
		ConnectionManager.server = new GatewayServer();  // TODO Dynamic port
		ConnectionManager.server.start();
        System.out.println("Gateway server started on port " + ConnectionManager.server.getPort() + "!");	// TODO
        
        PressCommunicator communicator = (PressCommunicator) ConnectionManager.server.getPythonServerEntryPoint(
        		new Class[] { PressCommunicator.class }
		);
        System.out.println("Communicator class initialized!");	// TODO
        
        // server.shutdown();
        return communicator;
	}
	
	
	public abstract class ConnectionManager {
		public static Process pythonProcess = null;
		public static GatewayServer server = null;
		public static PressCommunicator communicator = null;
		public static BarcodeDataset barcodesInUse = null;
		public static List<Integer> validBarcodesInUse = null;
		
		public static File tempFilesystem = new File(RAM_FILESYSTEM_PATH);
		
		static {
			if (!tempFilesystem.isDirectory()) {
				tempFilesystem = new File("/tmp");
			}
		}
		
		static {
			try {
				communicator = createCommunicator();
			}
			catch (Exception e) {
				// TODO
				// e.printStackTrace();
				System.out.println("py4j listener not found.");
			}
		}
	}
	
}
