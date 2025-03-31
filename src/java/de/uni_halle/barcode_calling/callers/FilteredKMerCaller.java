package de.uni_halle.barcode_calling.callers;

import java.util.Map;

import de.uni_halle.barcode_calling.callers.indexes.KMerIndex;
import de.uni_halle.barcode_calling.exceptions.MissingParameterException;
import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.datasets.BarcodeDataset;
import de.uni_halle.barcode_calling.util.datasets.Dataset;

public class FilteredKMerCaller extends IndexCaller {
	
	static final int DEFAULT_K_FILTER = 7;
	static final int DEFAULT_POS_FILTER = 100;
	static final int DEFAULT_K = 5;
	static final int DEFAULT_POS = 100;
	
	protected int pos, posFilter;
	protected KMerCaller filter;
	protected KMerCaller caller;
	
	public FilteredKMerCaller (int k, int maxDist, int pos, int kFilter, int maxDistFilter, int posFilter) {
		super(new KMerIndex(k));	// Dummy initialization
		this.filter = new KMerCaller(kFilter, maxDistFilter);
		this.caller = new KMerCaller(k, maxDist);
		// TODO set pos
		this.pos = pos;
		this.posFilter = posFilter;
		this.setMaxDist(maxDist);
		
		super.setIndex(this.caller.getIndex());
	}
	
	public FilteredKMerCaller (int k, int maxDist, int pos) {
		this(
				k, 
				maxDist, 
				pos, 
				DEFAULT_K_FILTER, 
				maxDist - 1,
				DEFAULT_POS_FILTER
		);
	}
	
	public FilteredKMerCaller (BarcodeDataset barcodes, int k, int kFilter) {
		this(
				k,
				IndexCaller.DEFAULT_DISTANCE_MEASURE.getDefaultMaxDistance(barcodes.getSequenceLength()),
				DEFAULT_POS,
				kFilter,
				IndexCaller.DEFAULT_DISTANCE_MEASURE.getDefaultMaxDistance(barcodes.getSequenceLength()) - 1,
				DEFAULT_POS_FILTER
		);
		this.setBarcodes(barcodes);
	}
	
	public FilteredKMerCaller (int k, int kFilter) {
		this(k, 0, DEFAULT_POS, kFilter, 0, DEFAULT_POS_FILTER);
	}
	
	public FilteredKMerCaller (int k) {
		this(k, 0, DEFAULT_POS);
	}
	
	public FilteredKMerCaller () {
		this(DEFAULT_K);
	}

	@Override
	public void setBarcodes(BarcodeDataset barcodes) {
		super.setBarcodes(barcodes);
		this.filter.setBarcodes(barcodes);
		this.caller.setBarcodes(barcodes);
		
		int barcodesCount = this.barcodes.getSize();
		double barcodesFractionProbed = this.pos / (double)barcodesCount;
		double barcodesFractionProbedFilter = this.posFilter / (double)barcodesCount;
		
		this.setBarcodesFractionProbed(barcodesFractionProbed);
		this.filter.setBarcodesFractionProbed(barcodesFractionProbedFilter);
	}
	
	@Override
	public void setBarcodesFractionProbed(double barcodesFractionProbed) {
		super.setBarcodesFractionProbed(barcodesFractionProbed);
		this.caller.setBarcodesFractionProbed(barcodesFractionProbed);
	}

	@Override
	public void setMaxDist(int maxDist) {
		super.setMaxDist(maxDist);
		this.caller.setMaxDist(maxDist);
	}
	
	public void setBarcodesFractionProbedFilter(double barcodesFractionProbed) {
		this.filter.setBarcodesFractionProbed(barcodesFractionProbed);
	}
	
	public void setMaxDistFilter(int maxDist) {
		this.filter.setMaxDist(maxDist);
	}
	
	@Override
	public void setDistanceModel(DistanceMeasureBase distanceMeasure) {
		super.setDistanceModel(distanceMeasure);
		this.filter.setDistanceModel(distanceMeasure);
		this.filter.setMaxDist(this.filter.getMaxDist() - 1);  // Increase precision
		this.caller.setDistanceModel(distanceMeasure);
	}
	
	/**
	 * Calls the read and returns the presumed barcode id and the SL distance
	 * @param read
	 * @return [callId, dist]
	 */
	@Override
	public int[] callDetails(DNASequence read) {
		long startTime = System.nanoTime();	// Track performance
		
		int[] res = this.filter.callDetails(read);
		
		long finishTime = System.nanoTime();	// Track performance
		this.addExecTimeIndex(finishTime - startTime);	// Track performance
		
		startTime = System.nanoTime();	// Track performance
		
		if (res[0] == -1) {
			res = this.caller.callDetails(read);
		}
		
		finishTime = System.nanoTime();	// Track performance
		this.addExecTime(finishTime - startTime);	// Track performance
		
		return res;
	}
	
	/**
	 * Result of calling read i is written into callingResults[i] and in case of call, the distance is written into callingDistances[i]
	 * @param reads
	 * @param callingResults
	 * @param callingDistances
	 */
	@Override
	public void callParallel(
			DNASequence[] reads, 
			int[] callingResults, 
			int[] callingDistances,
			int[] callingResultsSecond, 
			int[] callingDistancesSecond
		) {
		// FIXME Seems to be a lot slower somehow...
		this.resetExecutionTimes();	// Track performance
		this.filter.callParallel(reads, callingResults, callingDistances, callingResultsSecond, callingDistancesSecond);
		this.caller.callParallel(reads, callingResults, callingDistances, callingResultsSecond, callingDistancesSecond);
		
		this.setExecTimeIndex(this.filter.getExecTimeIndex() + this.caller.getExecTimeIndex());	// Track performance
		this.setExecTime(this.filter.getExecTime() + this.caller.getExecTime());	// Track performance
	}
	
	
	@Override
	public String getName() {
		StringBuilder sb = new StringBuilder(this.getClass().getSimpleName());
		
		sb.append("_k");
		sb.append(this.caller.getK());
		sb.append("_kFilter");
		sb.append(this.filter.getK());
		
		return sb.toString();
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(this.getName());
		
		sb.append("_L");
		sb.append(maxDist);
		sb.append("_pos");
		sb.append(this.getBarcodesFractionProbed());
		
		return sb.toString();
	}
	
	public static FilteredKMerCaller fromParameters(Map<String, Object> params, BarcodeDataset barcodes) {
		int k, kf;
		
		if (params.containsKey("k")) {
			k = (Integer)params.get("k");
		}
		else {
			throw new MissingParameterException("Parameter 'k' missing!");
		}
		
		if (params.containsKey("kf")) {
			kf = (Integer)params.get("kf");
		}
		else {
			throw new MissingParameterException("Parameter 'kf' missing!");
		}
		
		FilteredKMerCaller alg = new FilteredKMerCaller(barcodes, k, kf);
		
		if (params.containsKey("L")) {
			int maxDist = (Integer)params.get("L");
			alg.setMaxDist(maxDist);
		}
		
		if (params.containsKey("p")) {
			int pos = (Integer)params.get("p");
			alg.setBarcodesFractionProbed((double)pos / (double)barcodes.getSize());
		}
		
		if (params.containsKey("Lf")) {
			int maxDist = (Integer)params.get("Lf");
			alg.setMaxDistFilter(maxDist);
		}
		
		if (params.containsKey("pf")) {
			int pos = (Integer)params.get("pf");
			alg.setBarcodesFractionProbedFilter((double)pos / (double)barcodes.getSize());
		}
		
		return alg;
	}
	
}
