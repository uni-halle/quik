package de.uni_halle.barcode_calling.experiments;

import static de.uni_halle.barcode_calling.util.UtilMethods.getTimestampString;

import java.nio.file.Path;

/**
 * Base for experiments.
 * @author Riko Uphoff
 *
 */
public interface Experiment {
	
	public void doExperiment(Path outFilePath);
	
	public default void doExperiment(String outFilePath) {
		Path outFile = Path.of(outFilePath);
		this.doExperiment(outFile);
	}
	
	public String getExperimentName();
	
	public void setExperimentName(String newName);
	
	public default String getDefaultResultFileName() {
		StringBuilder sb = new StringBuilder(this.getExperimentName());
		sb.append("_");
		sb.append(getTimestampString());
		sb.append(".csv");
		
		return sb.toString();
	}
	
	/**
	 * Resolves outPath to a proper file path by appending the default result file name if it is a directory.
	 * Else it is left unchanged. 
	 */
	public default Path resolveOutPath(Path outPath) {
		
		Path resolvedOutPath = null;
		
		if (outPath.toFile().isFile()) {
			// File
			if (outPath.toString().endsWith(".csv")) {
				// CSV
				resolvedOutPath = outPath;
			}
			else {
				throw new RuntimeException("Cannot handle file extension of path " + outPath.toString());
			}
		}
		else {
			// Directory
			resolvedOutPath = Path.of(outPath.toString(), this.getDefaultResultFileName());
		}
		
		return resolvedOutPath;
	}
	
	
	/*
	 * Static implemented methods
	 */
	
}
