package de.uni_halle.barcode_calling.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class IOMethods {
	
	// Paths
	public static Path[] resolveInputPath(Path inputPath) {
		
		List<Path> barcodeFiles = new ArrayList<Path>();
		Path barcodeFile;
		
		if (inputPath.toFile().isDirectory()) {
			File[] files = inputPath.toFile().listFiles();
			
			for (int i = 0; i < files.length; i++) {
				
				barcodeFile = files[i].toPath();
				if (barcodeFile.toFile().isFile()) {
					barcodeFiles.add(barcodeFile);
				}
				
			}
		}
		else {
			barcodeFiles.add(inputPath);
		}
		
		return barcodeFiles.toArray(new Path[0]);
	}
	
	public static String getFileExtension(Path filePath) {
		String filePathString = filePath.getFileName().toString();
		String extension = null;
		
		int extensionSeperatorIndex = filePathString.lastIndexOf('.');
		if (extensionSeperatorIndex >= 0) {
			extension = filePathString.substring(extensionSeperatorIndex);
		}
		else if (filePath.toFile().isDirectory()) {
			throw new IllegalArgumentException("Path " + filePath + " is a directory!");
		}
		
		return extension;
	}
	
	// CSV
	// TODO Make separate class
	public static boolean writeCSV(Path csvFilePath, String[] header, Iterable<Object[]> entries) {
		return IOMethods.writeCSV(csvFilePath, header, entries, ',');
	}

	public static boolean writeCSV(Path csvFilePath, String[] header, Iterable<Object[]> entries, char separator) {
		int columnCount = header.length;
		
		try (FileWriter writer = new FileWriter(csvFilePath.toFile())) {
			writer.write(toCSVRow(header, separator));
			
			for (Object[] entry : entries) {
				if (entry.length != columnCount) {
					throw new RuntimeException("Attempted to write a row with different column count to csv!");
				}
				writer.write('\n');
				writer.write(toCSVRow(entry, separator));
			}
		}
		catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		catch (RuntimeException e) {
			e.printStackTrace();
			return false;
		}
		
		return true;
	}

	private static String toCSVRow(Object[] entry, char separator) {
		StringBuilder row = new StringBuilder();
		
		for (Object o : entry) {
			row.append(toCSVField(o));
			row.append(separator);
		}
		
		return row.substring(0, row.length() - 1);
	}

	private static String toCSVField(Object o) {
		String out;
		
		if (o instanceof Number) {
			out = o.toString();
		}
		else {
			String rawString = o.toString();
			out = "\"" + rawString.replace("\"", "\"\"") + "\"";
		}
		
		return out;
	}
	
}
