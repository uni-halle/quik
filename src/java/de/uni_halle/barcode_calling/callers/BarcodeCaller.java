package de.uni_halle.barcode_calling.callers;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.experiments.CallingExperiment;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.DatasetAnalyzer;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;
import de.uni_halle.barcode_calling.util.datasets.ReadDataset;

/**
 * Base class for barcode calling algorithms.
 * @author Riko Uphoff
 *
 */
public interface BarcodeCaller {
	
	/**
	 * Calls the read and returns the presumed barcode id and the distance for best and second best guess.
	 * This is the minimum a barcode caller has to implement.
	 * @param read
	 * @return [callIdBest, distBest, callIdSecond, distSecond]
	 */
	public int[] callDetails(DNASequence read);
	
	public default int call(DNASequence read) {
		int[] callDetails = this.callDetails(read);
		
		return callDetails[0];
	}
	
	/**
	 * Result of calling read i is written into callingResults[i] (-1 if miss) and in case of call, the distance is written into callingDistances[i].
	 * The same goes for the second best guess.
	 * @param reads
	 * @param i
	 * @param callingResults
	 * @param callingDistances
	 * @param callingResultsSecond
	 * @param callingDistancesSecond
	 */
	public default void call(
			DNASequence[] reads, 
			int i, 
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond
		) {
		
		DNASequence read = reads[i];
		int[] callingDetails = this.callDetails(read);
		
		callingResults[i] = callingDetails[0];
		callingDistances[i] = callingDetails[1];
		callingResultsSecond[i] = callingDetails[2];
		callingDistancesSecond[i] = callingDetails[3];
	}
	
	/**
	 * Result of calling read i is written into callingResults[i] and in case of call, the distance is written into callingDistances[i]
	 * @param reads
	 * @param i
	 * @param callingResults
	 * @param callingDistances
	 */
	public default void call(
			DNASequence[] reads, 
			int i, 
			int[] callingResults, 
			int[] callingDistances
		) {
		int[] dummy = new int[reads.length];
		
		this.call(
				reads, 
				i, 
				callingResults, 
				callingDistances, 
				dummy, 
				dummy
			);
	}
	
	/**
	 * Result of calling read i is written into callingResults[i] and in case of call, the distance is written into callingDistances[i]
	 * @param reads
	 * @param callingResults
	 * @param callingDistances
	 */
	public default void callParallel(
			DNASequence[] reads, 
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond
		) {
		this.resetExecutionTimes();	// Track performance
		
		int maxThreadCount = Main.maxThreadCount;
		ExecutorService threadPool = Main.threadPool;
		int stepSize = reads.length / maxThreadCount;
		Future<?>[] threadHandles = new Future<?>[maxThreadCount + 1];
		CallerThread[] threads = new CallerThread[maxThreadCount + 1];
		
		for (int t = 0; t < maxThreadCount; t++) {
			threads[t] = new CallerThread(
					this, 
					reads, 
					callingResults, 
					callingDistances,
					callingResultsSecond, 
					callingDistancesSecond, 
					t*stepSize, 
					(t+1)*stepSize
				);
			threadHandles[t] = threadPool.submit(threads[t]);
		}
		
		threads[maxThreadCount] = new CallerThread(
				this, 
				reads, 
				callingResults, 
				callingDistances,
				callingResultsSecond, 
				callingDistancesSecond, 
				maxThreadCount*stepSize, 
				reads.length
			);
		threadHandles[maxThreadCount] = threadPool.submit(threads[maxThreadCount]);
		
		try {
			
			for (Future<?> handle : threadHandles) {
				handle.get();
			}
			
			// Set exec time to longest running thread
			long maxExecTime = 0;
			for (CallerThread thread : threads) {
				maxExecTime = thread.execTime > maxExecTime ? thread.execTime : maxExecTime;
			}

			this.setExecTime(maxExecTime);
			
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
	};
	
	public default int[] callParallel(DNASequence[] reads) {
		int[] callingResults = new int[reads.length];
		for (int i = 0; i < reads.length; i++)
			callingResults[i] = -1;
		
		int[] dummy = new int[reads.length];
		
		this.callParallel(
				reads, 
				callingResults, 
				dummy,
				dummy, 
				dummy);
		
		return callingResults;
	}
	
	public default int[] call(DNASequence[] reads) {
		return this.callParallel(reads);
	}
	
	public default int[] call(Dataset reads) {
		return this.call(reads.getSequences());
	}
	
	/**
	 * Result of calling read i is written into callingResults[i] and in case of call, the distance is written into callingDistances[i]
	 * @param reads
	 * @param callingResults
	 * @param callingDistances
	 */
	public default void callParallel(
			Dataset reads,
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond
			) {
		this.callParallel(
				reads.getSequences(), 
				callingResults, 
				callingDistances,
				callingResultsSecond, 
				callingDistancesSecond
			);
	}
	
	public void resetExecutionTimes();
	
	public void addExecTime(long t);
	
	public void setExecTime(long t);
	
	public long getExecTime();
	
	public void setBarcodes(BarcodeDataset barcodes);
	
	public BarcodeDataset getBarcodes();
	
	public DNASequence getBarcode(int id);
	
	public default int getMaxDist() {
		return this.getBarcodes().getSequenceLength();
	}
	
	public default double getBarcodesFractionProbed() {
		return 1.0;
	}
	
	public default String getName() {
		return this.getClass().getSimpleName();
	}
	
	/**
	 * EXPERIMENTAL
	 * Sets the parameters of the model so that maximum F1-score can be expected with (if possible) precision and recall of at least minPres and minRec.
	 * @param reads Reads to tune parameters to.
	 * @param minPres Minimum desired precision.
	 * @param minRec Minimum desired recall.
	 * @return Array with [0] predicted precision and [1] predicted recall, and [2] predicted F1-score or NULL if minimum criteria can not be satisfied.
	 */
	public default double[] tuneParameters(DNASequence[] reads, double minPres, double minRec) {
		return this.tuneParameters(new ReadDataset(reads), minPres, minRec);
	}
	
	/**
	 * EXPERIMENTAL
	 * Sets the parameters of the model so that maximum F1-score can be expected with (if possible) precision and recall of at least minPres and minRec.
	 * @param reads Reads to tune parameters to.
	 * @param minPres Minimum desired precision.
	 * @param minRec Minimum desired recall.
	 * @return Array with [0] predicted precision, [1] predicted recall, and [2] predicted F1-score or NULL if minimum criteria can not be satisfied.
	 */
	public default double[] tuneParameters(ReadDataset reads, double minPres, double minRec) {
		
		LabeledDataset readsArtif = DatasetAnalyzer.generateSimilarReads(reads, this.getBarcodes());
		
		return this.tuneParameters(readsArtif, minPres, minRec);
	};
	
	/**
	 * EXPERIMENTAL
	 * Sets the parameters of the model so that maximum F1-score can be expected with (if possible) precision and recall of at least minPres and minRec.
	 * @param readsArtif Artificial reads to tune parameters on.
	 * @param minPres Minimum desired precision.
	 * @param minRec Minimum desired recall.
	 * @return Array with [0] predicted precision, [1] predicted recall, and [2] predicted F1-score or NULL if minimum criteria can not be satisfied.
	 */
	public default double[] tuneParameters(LabeledDataset readsArtif, double minPres, double minRec) {
		
		int[] calles = this.call(readsArtif);
		double[] callingRates = CallingExperiment.getMatchingRates(readsArtif.getLabels(), calles);
		double pres = callingRates[0];
		double rec = callingRates[1];
		double f1 = callingRates[2];
		
		double[] out = {pres, rec, f1};
		
		return out;
	};
	
	
	
	/**
	 * Thread class for parallel calling into input arrays that also outputs the distance to the closest barcode.
	 */
	public static class CallerThread extends Thread {
		private BarcodeCaller caller;
		DNASequence[] reads;
		private int[] callingResults;
		private int[] callingDistances;
		private int[] callingResultsSecond;
		private int[] callingDistancesSecond;
		private int start, end;
		public long execTime;
		
		public CallerThread(
				BarcodeCaller caller, 
				DNASequence[] reads, 
				int[] callingResults, 
				int[] callingDistances,
				int[] callingResultsSecond, 
				int[] callingDistancesSecond, 
				int start, 
				int end) {
			this.caller = caller;
			this.reads = reads;
			this.callingResults = callingResults;
			this.callingDistances = callingDistances;
			this.callingResultsSecond = callingResultsSecond;
			this.callingDistancesSecond = callingDistancesSecond;
			this.start = start;
			this.end = end;
		}
		
		@Override
		public void run() {
			long startTime = System.nanoTime();
			for (int i = start; i < end; i++) {
				try {
					if (callingResults[i] == -1)
						caller.call(reads, i, callingResults, callingDistances, callingResultsSecond, callingDistancesSecond);
				}
				catch (NullPointerException e) {
					System.out.println("Read " + i + " contains null!");
					for (int j = 0; j < this.reads[i].getLength(); j++) {
						if (this.reads[i].at(j) == null) {
							System.out.print("(NULL)");
						}
						else {
							System.out.print(this.reads[i].at(j));
						}
					}
				}
			}
			long finishTime = System.nanoTime();
			this.execTime = finishTime - startTime;
		}
	}
}
