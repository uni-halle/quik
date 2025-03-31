package de.uni_halle.barcode_calling.callers.indexes;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.TrimerCaller;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

public class TrimerIndexOutsourced implements OrderedIndex {
	
	protected static final double DEFAULT_BARCODES_FRACTION_PROBED = 0.01;
	
	protected double barcodesFractionProbed;
	protected BarcodeDataset barcodes;
	private PressCommunicator communicator;
	
	public TrimerIndexOutsourced(BarcodeDataset barcodes, double barcodesFractionProbed, PressCommunicator communicator) {
		this.communicator = communicator;
		this.barcodesFractionProbed = barcodesFractionProbed;
		this.setBarcodes(barcodes);
	}
	
	public TrimerIndexOutsourced(BarcodeDataset barcodes, double barcodesFractionProbed) {
		this(barcodes, barcodesFractionProbed, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerIndexOutsourced(BarcodeDataset barcodes, PressCommunicator communicator) {
		this(barcodes, DEFAULT_BARCODES_FRACTION_PROBED, communicator);
	}
	
	public TrimerIndexOutsourced(BarcodeDataset barcodes) {
		this(barcodes, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerIndexOutsourced(double barcodesFractionProbed, PressCommunicator communicator) {
		this(new BarcodeDataset(), barcodesFractionProbed, communicator);
	}
	
	public TrimerIndexOutsourced(double barcodesFractionProbed) {
		this(barcodesFractionProbed, PressCommunicator.ConnectionManager.communicator);
	}
	
	public TrimerIndexOutsourced(PressCommunicator communicator) {
		this(new BarcodeDataset(), communicator);
	}
	
	public TrimerIndexOutsourced() {
		this(PressCommunicator.ConnectionManager.communicator);
	}
	
	@Override
	public void setBarcodes(BarcodeDataset barcodes) {
		this.barcodes = barcodes;
		
		// Write barcodes to temp file
		try {
			File barcodesFile = File.createTempFile("barcodes", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			barcodes.toFile(barcodesFile.getAbsolutePath());
			this.communicator.setBarcodesFile(barcodesFile.getAbsolutePath());
			barcodesFile.delete();
		}
		catch (IOException e) {
			Main.writeToLog(e);
		}
	}
	
	@Override
	public BarcodeDataset getBarcodes() {
		return this.barcodes;
	}
	
	@Override
	public double getBarcodesFractionProbed() {
		return this.barcodesFractionProbed;
	}

	@Override
	public void setBarcodesFractionProbed(double barcodesFractionProbed) {
		this.barcodesFractionProbed = barcodesFractionProbed;
	}

	@Override
	public double[] getRecallAtPos(LabeledDataset reads) {
		int barcodesCount = this.barcodes.getSize();
		
		if (reads.getSize() == 0)
			return new double[barcodesCount];
		
		int[] posCount = null;
		
		try {
			// Temp files for communication
			File readsFile = File.createTempFile("reads", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			File labelsFile = File.createTempFile("labels", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			File outFile = File.createTempFile("out", ".txt", PressCommunicator.ConnectionManager.tempFilesystem);
			
			// Write reads to temp file
			reads.toFile(readsFile.getPath());
			TrimerCaller.writeBinaryIntFile(reads.getLabels(), labelsFile);
			
			// Run external Python script to get count of correct barcodes at positions in list
			this.communicator.writeFrequencyAtPos(readsFile.getAbsolutePath(), labelsFile.getAbsolutePath(), outFile.getAbsolutePath());
			
			posCount = TrimerCaller.readBinaryIntFile(outFile, barcodesCount);
			
			readsFile.delete();
			labelsFile.delete();
			outFile.delete();
			
		} catch (Exception e) {
			Main.writeToLog(e);
		}
		
		double recall = 0;
		double[] recallAtPos = new double[barcodesCount];
		
		for (int pos = 0; pos < barcodesCount; pos++) {
			recall += (double)posCount[pos] / (double)reads.getSize();
			recallAtPos[pos] = recall;
		}
		
		return recallAtPos;
	}
	
	@Override
	public String getName() {
		return "TrimerIndex";
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(this.getName());
		
		sb.append("_pos");
		sb.append(this.barcodesFractionProbed);
		
		return sb.toString();
	}

	@Override
	public int[] getOrderedCandidates(DNASequence seq) {
		// TODO Auto-generated method stub
		return null;
	}

	
}
