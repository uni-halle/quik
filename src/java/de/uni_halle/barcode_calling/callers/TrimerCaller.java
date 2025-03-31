package de.uni_halle.barcode_calling.callers;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import de.uni_halle.barcode_calling.callers.indexes.PressCommunicator;
import de.uni_halle.barcode_calling.callers.indexes.TrimerIndexOutsourced;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;


/**
 * CPU implementation of Press (2022).
 * @author Riko Uphoff
 *
 */
public class TrimerCaller extends IndexCaller {
	
	/*
	 * Instance variables
	 */
	PressCommunicator communicator;

	/*
	 * Constructors
	 */
	public TrimerCaller(BarcodeDataset barcodes, int maxDist, double indexFractionProbed, PressCommunicator communicator) {
		super(barcodes, new TrimerIndexOutsourced(indexFractionProbed), maxDist, indexFractionProbed);
		
		this.communicator = communicator;
	}
	
	public TrimerCaller(BarcodeDataset barcodes, int maxDist, PressCommunicator communicator) {
		this(barcodes, maxDist, 1.0, communicator);
	}
	
	public TrimerCaller(BarcodeDataset barcodes, double indexFractionProbed, PressCommunicator communicator) {
		this(barcodes, barcodes.getMinMutualEditDist() / 2, indexFractionProbed, communicator);
	}
	
	public TrimerCaller(BarcodeDataset barcodes, PressCommunicator communicator) {
		this(barcodes, barcodes.getMinMutualEditDist() / 2, communicator);
	}
	
	public TrimerCaller(BarcodeDataset barcodes, int maxDist, double indexFractionProbed) {
		this(barcodes, maxDist, indexFractionProbed, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerCaller(BarcodeDataset barcodes, int maxDist) {
		this(barcodes, maxDist, 1.0, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerCaller(BarcodeDataset barcodes, double indexFractionProbed) {
		this(barcodes, barcodes.getMinMutualEditDist() / 2, indexFractionProbed, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerCaller(BarcodeDataset barcodes) {
		this(barcodes, barcodes.getMinMutualEditDist() / 2, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerCaller() {
		this(new BarcodeDataset());
	}
	
	@Override
	public String getName() {
		return this.getClass().getSimpleName();
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(this.getName());
		
		sb.append("_L");
		sb.append(maxDist);
		
		sb.append("_pos");
		sb.append(this.index.getBarcodesFractionProbed());
		
		return sb.toString();
	}
	
	@Override
	public int call(DNASequence read) {
		int[] callDetails = this.callDetails(read);
		
		return callDetails[0];
	}
	
	@Override
	public int[] callDetails(DNASequence read) {
		System.err.println("Warning: calling single reads can be very inefficient with the current implementation of TrimerCaller. Consider calling in batches or switching to TrimerCallerCPU class.");
		DNASequence[] reads = {read};
		int[][] res = this.callDetails(reads);
		int[] out = {res[0][0], res[1][0]};
		
		return out;
	}
	
	/**
	 * Calls the reads and returns the presumed barcode ids and the Levenshtein distances
	 * @param read
	 * @return [callIds[], dists[]]
	 */
	public int[][] callDetails(DNASequence[] reads) {
		this.resetExecutionTimes();
		if (reads.length == 0)
			return new int[0][0];
		
		int barcodesCount = this.barcodes.getSize();
		int readsCount = reads.length;
		int[] labels = null;
		int[] dists = null;
		
		try {
			// Temp files for communication
			File readsFile = File.createTempFile("reads", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			File resultFile = File.createTempFile("decodes", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			File execTimesFile = File.createTempFile("exec_times", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			File distancesFile = File.createTempFile("distances", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			
			// Write reads to temp file
			Dataset readsSet = new Dataset(reads);
			readsSet.toFile(readsFile.getPath());
			
			// Run external Python script
			int pos = (int)(this.getBarcodesFractionProbed() * barcodesCount);
			this.communicator.decodeReads(
					readsFile.getAbsolutePath(), 
					resultFile.getAbsolutePath(), 
					execTimesFile.getAbsolutePath(), 
					distancesFile.getAbsolutePath(), 
					pos, this.maxDist
					);
			
			long[] execTimes = readExecTimesFile(execTimesFile);
			this.addExecTimeIndex(execTimes[0]);
			this.addExecTime(execTimes[1]);
			
		    // Load decoded reads
		    labels = readBinaryIntFile(resultFile, readsCount);
		    dists = readBinaryIntFile(distancesFile, readsCount);
		    
		    // Delete files
		    readsFile.delete();
		    resultFile.delete();
		    distancesFile.delete();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return new int[][] {labels, dists};
	}
	
	public int[] call(DNASequence[] reads) {
		return this.callDetails(reads)[0];
	}
	
	@Override
	public void call(DNASequence[] reads, int i, int[] callingResults, int[] callingDistances) {
		int[] callDetails = this.callDetails(reads[i]);
		callingResults[i] = callDetails[0];
		callingDistances[i] = callDetails[1];
	}
	
	public int[] callParallel(DNASequence[] reads) {
		// TODO Parallel GPU support
		return this.call(reads);
	}
	
	@Override
	public void callParallel(
			DNASequence[] reads, 
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond
		) {
		// TODO Add second best guess
		// TODO Add parallel support for multiple GPUs
		int[][] callDetails = this.callDetails(reads);
		
		for (int i = 0; i < callingResults.length; i++) {
			callingResults[i] = callDetails[0][i];
			callingDistances[i] = callDetails[1][i];
		}
	};
	
	/*
	 * Static methods
	 */
	public static int[] readBinaryIntFile(File intFile, int numints) throws FileNotFoundException, IOException {
		int[] out = new int[numints];
		byte[] intBuffer = new byte[4];
		
		try (
	            InputStream inputStream = new BufferedInputStream(new FileInputStream(intFile.getAbsolutePath()));
	    ) {
	    	DataInputStream dataStream = new DataInputStream(inputStream);
	    	
	    	for (int i = 0; i < numints; i++) {
	    		dataStream.read(intBuffer);
    			int newInt = (((intBuffer[3] & 0xff) << 24) | ((intBuffer[2] & 0xff) << 16) | ((intBuffer[1] & 0xff) <<  8) | (intBuffer[0] & 0xff));
    			out[i] = newInt;
	    	}
		}
		
		return out;
	}
	
	public static void writeBinaryIntFile(int[] data, File outFile) throws FileNotFoundException, IOException {
		byte[] intBuffer = new byte[4];
		
		try (
	            OutputStream outputStream = new BufferedOutputStream(new FileOutputStream(outFile.getAbsolutePath()));
	    ) {
	    	DataOutputStream dataStream = new DataOutputStream(outputStream);
	    	
	    	for (int x = 0; x < data.length; x++) {
	    		intBuffer[0] = (byte)(data[x]);
	    		intBuffer[1] = (byte)(data[x] >> 8);
	    		intBuffer[2] = (byte)(data[x] >> 16);
	    		intBuffer[3] = (byte)(data[x] >> 24);
	    		
	    		dataStream.write(intBuffer);
    		}
		}
	}
	
	public static long[] readExecTimesFile(File file) {
		long[] out = new long[2];
		
		try (
				BufferedReader reader = new BufferedReader(new FileReader(file))
		) {
			for (int i = 0; i < out.length; i++)
				out[i] = (long)(Double.parseDouble(reader.readLine()) * 1000000000);	// To nano seconds
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		
		return out;
	}

}