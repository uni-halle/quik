package de.uni_halle.barcode_calling.callers.indexes;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.UtilMethods;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;
import de.uni_halle.barcode_calling.util.datasets.DatasetAnalyzer;
import de.uni_halle.barcode_calling.util.datasets.LabeledDataset;

/**
 * Index which orders barcodes after their probability of being a match.
 * @author Riko Uphoff
 *
 */
public interface OrderedIndex {
	
	public void setBarcodes(BarcodeDataset barcodes);
	
	public BarcodeDataset getBarcodes();
	
	public int[] getOrderedCandidates(DNASequence seq);
	
	public default String getName() {
		return this.getClass().getSimpleName();
	}
	
	public double getBarcodesFractionProbed();
	
	public void setBarcodesFractionProbed(double barcodesFractionProbed);
	
	/**
	 * Sets the parameters of the index so that the position of the correct candidate is as low as possible with recall of at least minRec.
	 * @param readsArtif Artificial reads to tune parameters to.
	 * @param minRec Minimum desired recall.
	 */
	public default void tuneParameters(LabeledDataset readsArtif, double minRec) {
		
		System.out.println("Tuning index " + this.getName());
		
		int barcodesCount = this.getBarcodes().getSize();
		double[] recallAtPos = getRecallAtPos(readsArtif);
		
		for (int pos = 0; pos < barcodesCount; pos++) {
			
			if (pos >= 9 && recallAtPos[pos] >= minRec) {
				this.setBarcodesFractionProbed((double)(pos+1) / (double)barcodesCount);
				break;
			}
		}
		
		// TODO
		System.out.println("Tuning of index yielded pos = " + this.getBarcodesFractionProbed());
	}
	
	/**
	 * Sets the parameters of the index so that the position of the correct candidate is as low as possible with recall of at least minRec.
	 * @param reads Reads to tune parameters to.
	 * @param minRec Minimum desired recall.
	 */
	public default void tuneParameters(Dataset reads, double minRec) {
		
		BarcodeDataset barcodes = this.getBarcodes();
		LabeledDataset readsArtif = DatasetAnalyzer.generateSimilarReads(reads, barcodes);
		this.tuneParameters(readsArtif, minRec);
		
	}
	
	/**
	 * Returns an array in which index i corresponds to the recall at position i.
	 * @param reads
	 * @return
	 */
	public default double[] getRecallAtPos(LabeledDataset reads) {
		
		// For parallelization
		int maxThreadCount = Main.maxThreadCount;
		ExecutorService threadPool = Main.threadPool;
		int stepSize = reads.getSize() / maxThreadCount;
		Thread thread;
		Future<?>[] threadHandles = new Future<?>[maxThreadCount + 1];
		
		// Collect positions of labels in index
		int barcodesCount = this.getBarcodes().getSize();
		int readsCount = reads.getSize();
		int[][] posCounts = new int[maxThreadCount+1][barcodesCount];
		int[] candidatesCounts = new int[maxThreadCount+1];
		
		for (int t = 0; t < maxThreadCount; t++) {
			thread = new OrderedIndexThread(this, reads, posCounts, candidatesCounts, t*stepSize, (t+1)*stepSize, t);
			threadHandles[t] = threadPool.submit(thread);
		}
		
		thread = new OrderedIndexThread(this, reads, posCounts, candidatesCounts, maxThreadCount*stepSize, readsCount, maxThreadCount);
		threadHandles[maxThreadCount] = threadPool.submit(thread);
		
		try {
			
			for (Future<?> handle : threadHandles) {
				handle.get();
			}
			
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		
		// Merge results
		int[] posCount = new int[barcodesCount];
		long candidatesCount = 0;
		
		for (int t = 0; t <= maxThreadCount; t++) {
			for (int pos = 0; pos < barcodesCount; pos++) {
				posCount[pos] += posCounts[t][pos];
			}
			candidatesCount += candidatesCounts[t];
		}
		
		// Calculate recall at positions and stop once we attain minRec
		double recall = 0;
		
		double[] recallAtPos = new double[barcodesCount];
		
		for (int pos = 0; pos < barcodesCount; pos++) {
			recall += (double)posCount[pos] / (double)readsCount;
			recallAtPos[pos] = recall;
		}
		
		// Log average candidates count
		Main.writeToLog("Average number of candidates for " + this.getName() + ": " + (candidatesCount / (double)readsCount));
		
		return recallAtPos;
	}
	
	
	/**
	 * Thread class for parallel matching into input arrays.
	 */
	public static class OrderedIndexThread extends Thread {
		private OrderedIndex index;
		private LabeledDataset reads;
		private int[][] posCounts;
		private int[] candidatesCounts;
		private int start, end;
		private int threadId;
		
		public OrderedIndexThread(OrderedIndex index, LabeledDataset reads, int[][] posCounts, int[] candidatesCounts, int start, int end, int threadId) {
			this.index = index;
			this.reads = reads;
			this.posCounts = posCounts;
			this.candidatesCounts = candidatesCounts;
			this.start = start;
			this.end = end;
			this.threadId = threadId;
		}
		
		@Override
		public void run() {
			int posLabel;
			
			for (int i = start; i < end; i++) {
				DNASequence read = reads.getSequence(i);
				int[] orderedCandidates = index.getOrderedCandidates(read);
				this.candidatesCounts[this.threadId] += orderedCandidates.length;
				posLabel = UtilMethods.indexOf(orderedCandidates, reads.getLabel(i));
				
				if (posLabel != -1) {
					this.posCounts[this.threadId][posLabel]++;
				}
			}
		}
	}
	
}
