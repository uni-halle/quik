package de.uni_halle.barcode_calling.util.datasets;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.callers.BarcodeCaller;
import de.uni_halle.barcode_calling.callers.DistanceMeasureBase;
import de.uni_halle.barcode_calling.callers.IndexCaller;
import de.uni_halle.barcode_calling.callers.KMerCaller;
import de.uni_halle.barcode_calling.callers.SequenceLevenshteinBase;
import de.uni_halle.barcode_calling.callers.indexes.KMerIndex;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.UtilMethods;

/**
 * Abstract class that provides the analytic functionalities from Section 3.4.5
 * @author Riko Uphoff
 *
 */
public abstract class DatasetAnalyzer {
	
	public static final int DEFAULT_K = 10000;
	public static final double DEFAULT_MAX_ERROR_P = 0.001;
	public static final DistanceMeasureBase DEFAULT_DISTANCE_MEASURE = new SequenceLevenshteinBase();
	
	/**
	 * Returns an array with observed error probabilities in the tau-calls amongst the reads and stores the tau-calls in the provided array tauCalls.
	 * 0: Fraction of reads which are tau-calls
	 * 1: Probability of substitution
	 * 2: Probability of insertion
	 * 3: Probability of deletion
	 */
	public static double[] analyzeTauCalls(Dataset reads, BarcodeDataset barcodes, int[] tauCalls) {
		
		int tau = UtilMethods.max(barcodes.getMinMutualEditDist() / 2, 3);	// TODO Mention
		int k = Math.min(DEFAULT_K, reads.getSize());
		int barcodesCount = barcodes.getSize();
		Dataset readsSample = reads.getSubset(k);
		
		// Find tau-calls
		IndexCaller caller = new KMerCaller(
				barcodes,
				7,
				tau,
				2/(double)barcodesCount
		);

		
		// TODO
		System.out.println("tau-Caller: " + caller.getName());
		
		// Inits
		double pSub = 0, pIns = 0, pDel = 0, count = 0;
		double[] distanceDetails;
		DNASequence read, barcode;
		int[] dummy = new int[k];	// Detailed information not needed
		
		caller.callParallel(readsSample, tauCalls, dummy, dummy, dummy);
		
		// Get error probs. for tau-calls
		for (int i = 0; i < k; i++) {
			
			if (tauCalls[i] != -1) {
				read = reads.getSequence(i);
				barcode = barcodes.getSequence(tauCalls[i]);
				
				distanceDetails = SequenceLevenshteinBase.inDistanceDetails(barcode, read, tau);
				
				if (distanceDetails[0] > 0) {
					pSub += distanceDetails[1] / distanceDetails[0];
					pIns += distanceDetails[2] / distanceDetails[0];
					pDel += distanceDetails[3] / distanceDetails[0];
				}
				count++;
			}
		}
		
		// Normalize
		pSub /= count;
		pIns /= count;
		pDel /= count;
		double fracTauCalls = count / (double)k;
		
		double[] analysis = {fracTauCalls, pSub, pIns, pDel};
		
		return analysis;
	}
	
	/**
	 * Returns an array with observed error probabilities in the tau-calls amongst the reads.
	 * 0: Fraction of reads which are tau-calls
	 * 1: Probability of substitution
	 * 2: Probability of insertion
	 * 3: Probability of deletion
	 */
	public static double[] analyzeTauCalls(Dataset reads, BarcodeDataset barcodes) {
		return analyzeTauCalls(reads, barcodes, new int[reads.getSize()]);
	}
	
	public static double[] estimateErrorProbabilities(Dataset reads, BarcodeDataset barcodes, double maximumErrorP) {
		/**
		 * Returns an array with the estimated error probabilities in the reads:
		 * 0: Probability of substitution
		 * 1: Probability of insertion
		 * 2: Probability of deletion
		 */
		int tau = barcodes.getMinMutualEditDist() / 2;
		double[] tauCallAnalysis = analyzeTauCalls(reads, barcodes);
		double fracTauCalls = tauCallAnalysis[0];
		double pSub = tauCallAnalysis[1];
		double pIns = tauCallAnalysis[2];
		double pDel = tauCallAnalysis[3];
		
		// Get actual error prob.
		double p = 0.5;
		double k = DEFAULT_K;
		int maxI = -(int)(Math.log(maximumErrorP) / Math.log(2.0)) - 1;
		DNASequence seq = DNASequence.random(barcodes.getSequenceLength());
		DNASequence seqMut;
		double countArtif = 0;
		double fracTauCallsArtif;
		
		for (int i = 1; i <= maxI; i++) {
			
			countArtif = 0;
			
			for (int j = 0; j < k; j++) {
				seqMut = seq.corrupt(p*pSub, p*pIns, p*pDel);
				
				if (DEFAULT_DISTANCE_MEASURE.inDistance(seq, seqMut, tau)) {
					countArtif++;
				}
			}
			
			fracTauCallsArtif = countArtif / k;
			
			if (fracTauCallsArtif > fracTauCalls) {
				p += Math.pow(2, -(i+1));
			}
			else {
				p -= Math.pow(2, -(i+1));
			}
		}
		
		pSub *= p;
		pIns *= p;
		pDel *= p;
		
		double[] errorProbabilities = {pSub, pIns, pDel};
		
		// Log estimates
		StringBuilder logText = new StringBuilder();
		logText.append("Estimated error probabilities:");
		logText.append("\n - Substitution: " + pSub);
		logText.append("\n - Insertion: " + pIns);
		logText.append("\n - Deletion: " + pDel);
		logText.append("\n -> Total: " + p);
		Main.writeToLog(logText.toString());
		
		return errorProbabilities;
	}
	
	public static double[] estimateErrorProbabilities(Dataset reads, BarcodeDataset barcodes) {
		/**
		 * Returns an array with the estimated error probabilities in the reads:
		 * 0: Probability of substitution
		 * 1: Probability of insertion
		 * 2: Probability of deletion
		 */
		return estimateErrorProbabilities(reads, barcodes, DEFAULT_MAX_ERROR_P);
	}
	
	public static LabeledDataset generateSimilarReads(Dataset reads, BarcodeDataset barcodes, double maximumErrorP) {
		
		double[] errorProbs = estimateErrorProbabilities(reads, barcodes, maximumErrorP);
		double pSub = errorProbs[0];
		double pIns = errorProbs[1];
		double pDel = errorProbs[2];
		
		int k = DEFAULT_K;
		
		LabeledDataset readsArtif = barcodes.corrupt(pSub, pIns, pDel, k);
		
		return readsArtif;
	}
	
	public static LabeledDataset generateSimilarReads(Dataset reads, BarcodeDataset barcodes) {
		return generateSimilarReads(reads, barcodes, DEFAULT_MAX_ERROR_P);
	}
	
}
