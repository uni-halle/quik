package de.uni_halle.barcode_calling.experiments;

import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;

/**
 * Abstract base for matching-/calling- experiments providing basic functionalities.
 * @author rikso
 *
 */
public abstract class CallingExperiment implements Experiment {
	
	private String experimentName = "matching-experiment";
	
	
	public String getExperimentName() {
		return this.experimentName;
	}
	
	public void setExperimentName(String newName) {
		this.experimentName = newName;
	}
	
	/**
	 * Calculates precision [0], recall [1] and f1 [2] for given data.
	 */
	public static double[] getMatchingRates(int[] labels, int[] predictions) {
		
		if (labels.length != predictions.length) {
			throw new IllegalArgumentException("labels and matches are of different size!");
		}
		
		double size = labels.length;
		double matches = 0;
		double misses = 0;
		
		for (int i = 0; i < labels.length; i++) {
			if (labels[i] < 0) {
				throw new RuntimeException("label contains a negative value (miss)...");
			}
			else if (predictions[i] < 0) {
				misses++;
			}
			else if (labels[i] == predictions[i]) {
				matches++;
			}
		}
		
		double precision = size > misses ? matches / (size - misses) : 1;
		double recall = (size - misses) / size;
		double f1 = precision+recall != 0 ? 2*precision*recall / (precision+recall) : 0;
		
		double[] out = {
			precision,
			recall,
			f1
		};
		
		return out;
	}
	
	/**
	 * Picks a default dataset among the available barcodes by choosing the one "closest" to the average of all variables.
	 * More specifically, the set with the minimum summed difference to the average of all variables is selected.
	 */
	public static BarcodeDataset getRepresentativeBarcodeSet(BarcodeDataset[] barcodeSets) {
		
		int avgN = 0, avgL = 0, avgD = 0;
		int count = barcodeSets.length;
		
		for (BarcodeDataset barcodeSet : barcodeSets) {
			avgN += barcodeSet.getSize();
			avgL += barcodeSet.getSequenceLength();
			avgD += barcodeSet.getMinMutualEditDist();
		}
		
		avgN = avgN / count;
		avgL = avgL / count;
		avgD = avgD / count;
		
		return getRepresentativeBarcodeSet(barcodeSets, avgN, avgL, avgD);
	}
	
	/**
	 * Picks a default dataset among the available barcodes by choosing the one "closest" to the average of all variables.
	 * More specifically, the set with the minimum summed difference to the average of all variables is selected.
	 */
	public static BarcodeDataset getRepresentativeBarcodeSet(
			BarcodeDataset[] barcodeSets, int defaultN, int defaultL, int defaultD
		) {
		
		BarcodeDataset defaultSet = null;
		double diffN, diffL, diffD;
		double diffTotal;
		double minDif = Double.MAX_VALUE;
		
		for (BarcodeDataset barcodeSet : barcodeSets) {
			diffN = (double)barcodeSet.getSize() / (double)defaultN;
			diffL = (double)barcodeSet.getSequenceLength() / (double)defaultL;
			diffD = (double)(barcodeSet.getMinMutualEditDist() + 1) / (double)(defaultD + 1);
			
			diffTotal = Math.abs(diffN - 1) + Math.abs(diffL - 1) + Math.abs(diffD - 1);
			
			if (diffTotal < minDif) {
				defaultSet = barcodeSet;
				minDif = diffTotal;
			}
		}
		
		return defaultSet;
	}
	
}
