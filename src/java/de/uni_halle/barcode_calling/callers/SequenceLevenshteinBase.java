package de.uni_halle.barcode_calling.callers;

import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;

import de.uni_halle.barcode_calling.util.DNASequence;
import static de.uni_halle.barcode_calling.util.UtilMethods.max;
import static de.uni_halle.barcode_calling.util.UtilMethods.min;

/**
 * Base for Sequence-Levenshtein baselines (Buschmann et al. 2013)
 * @author Riko Uphoff
 *
 */
public class SequenceLevenshteinBase implements DistanceMeasureBase {
	
	public static final String NAME = "Sequence-Levenshtein";
	
	public int getDefaultMaxDistance(int length) {
		// TODO This has not been tested thoroughly
		return (int)Math.round(length / 4.5);
	}
	
	/**
	 * Classical dynamic program for computing Levensthein distance extended to return the SL distance
	 * @param a Sequence A
	 * @param b Sequence B
	 * @return Sequence-Levenshtein distance between A and B.
	 */
	public int distance(DNASequence a, DNASequence b) {
		int lengthA = a.getLength();
		int lengthB = b.getLength();
		
		// Distance matrix
		int[][] dm = LevenshteinBase.computeDistMatrix(a, b);
		
		// Sequence part
		// Elongating
		int minDist = min(dm[lengthA]);
		
		// Truncating
		for (int i = 0; i <= lengthA; i++) {
			minDist = min(minDist, dm[i][lengthB]);
		}
	
		return minDist;
	}
	
	/**
	 * Implements the thresholded Levenshtein calculation after Zorita et al. (2015) extended to return the SL distance
	 * @param a Sequence A
	 * @param b Sequence B
	 * @param maxDist Maximum computable distance
	 * @return Sequence-Levenshtein distance between a and b up to maxDist
	 */
	public int inDistanceResult(DNASequence a, DNASequence b, int maxDist) {
		// Account for negative maxDist
		if (maxDist < 0) return 0;
		
		// Use shorter string as reference
		if (a.getLength() > b.getLength()) {
			DNASequence temp = a;
			a = b;
			b = temp;
		}
		
		int minLength = a.getLength();
		int lengthDif = b.getLength() - a.getLength();
		int posA, posB, editCost;
		int slDistance = 0;	// Initial value for sequences of length 0
		int columnSize = maxDist+1 > lengthDif ? maxDist+1 - lengthDif : 0;
		int rowSize = maxDist+1 + lengthDif;
		int[][] columns = new int[minLength+1][columnSize];	// Column is indexed top to bottom
		int[][] rows = new int[minLength+1][rowSize];		// Row is indexed left to right
		int[] centers = new int[minLength+1];
		
		centers[0] = lengthDif;
		
		
		for (int depth = 0; depth <= minLength; depth++) {
			slDistance = Integer.MAX_VALUE;
			
			// Column
			for (int i = 0; i < columnSize; i++) {
				if (i < max(1, columnSize+1 - depth) || depth == 0) {
					// Border
					columns[depth][i] = Math.abs(maxDist+1 - i);
				}
				else {
					// Inside
					posA = depth-1 - columnSize + i;
					posB = depth-1 + lengthDif;
					
					editCost = a.at(posA) == b.at(posB) ? 0 : 1;
					
					columns[depth][i] = min(
							columns[depth-1][i] + editCost,	// replace
							columns[depth][i-1] + 1, 				// delete
							i == columnSize-1 ? centers[depth-1] + 1 : columns[depth-1][i+1] + 1	// insert
					);
				}
				
				// Minimum for Sequence Levensthein distance
				slDistance = columns[depth][i] < slDistance ? columns[depth][i] : slDistance;
			}
			
			// Row
			for (int i = 0; i < rowSize; i++) {
				if (i < max(1, rowSize+1 - depth - lengthDif) || depth == 0) {
					// Border
					rows[depth][i] = Math.abs(maxDist+1 - i);
				}
				else {
					// Inside
					posA = depth-1;
					posB = depth-1 - rowSize + i + lengthDif;
					
					editCost = a.at(posA) == b.at(posB) ? 0 : 1;
					
					rows[depth][i] = min(
							rows[depth-1][i] + editCost,	// replace
							rows[depth][i-1] + 1, 			// insert
							i == rowSize-1 ? centers[depth-1] + 1 : rows[depth-1][i+1] + 1	// delete
					);
				}
				
				// Minimum for Sequence Levensthein distance
				slDistance = rows[depth][i] < slDistance ? rows[depth][i] : slDistance;
			}
			
			// Center
			if (maxDist+1 <= lengthDif || depth == 0) {
				// Border
				centers[depth] = lengthDif;
			}
			else {
				// Inside
				posA = depth-1;
				posB = depth-1 + lengthDif;
				editCost = a.at(posA) == b.at(posB) ? 0 : 1;
				
				centers[depth] = min(
						centers[depth-1] + editCost,		// replace
						rows[depth][rowSize-1] + 1,			// insert
						columns[depth][columnSize-1] + 1	// delete
				);
			}
			
			// Minimum for Sequence Levensthein distance
			slDistance = centers[depth] < slDistance ? centers[depth] : slDistance;
			
			// Early stopping
			if (slDistance > maxDist) {
				return slDistance;
			}
		}
		
		return slDistance;
	}
	
	
	/**
	 * Implements the thresholded Levenshtein calculation after Zorita et al. (2015) extended to return the SL distance
	 * @param a Sequence A
	 * @param b Sequence B
	 * @param maxDist Maximum computable distance
	 * @return Are A and B within a Sequence-Levenshtein distance of maxDist?
	 */
	public boolean inDistance(DNASequence a, DNASequence b, int maxDist) {
		return inDistanceResult(a, b, maxDist) <= maxDist;
	}
	
	/**
	 * Returns an array with detailed information about the Sequence-Levenshtein distance: 
	 * 0: Sequence-Levenshtein distance
	 * 1: Average number of substitutions
	 * 2: Average number of insertions
	 * 3: Average number of deletions
	 */
	public static double[] distanceDetails(double[][][] detailedDistMatrix) {
		
		int lengthA = detailedDistMatrix[0].length;
		int lengthB = detailedDistMatrix[0][0].length;
		int posA = lengthA-1;
		int posB = lengthB-1;
		double minDist = detailedDistMatrix[0][posA][posB];
		
		for (int a = 0; a < lengthA-1; a++) {
			if (detailedDistMatrix[0][a][lengthB-1] < minDist) {
				posA = a;
				posB = lengthB-1;
				minDist = detailedDistMatrix[0][posA][posB];
			}
		}
		for (int b = 0; b < lengthA-1; b++) {
			if (detailedDistMatrix[0][lengthA-1][b] < minDist) {
				posA = lengthA-1;
				posB = b;
				minDist = detailedDistMatrix[0][posA][posB];
			}
		}
		
		double[] out = {
				detailedDistMatrix[0][posA][posB],
				detailedDistMatrix[1][posA][posB],
				detailedDistMatrix[2][posA][posB],
				detailedDistMatrix[3][posA][posB]
		};
		
		return out;
	}
	
	/**
	 * Returns an array with detailed information about the Sequence-Levenshtein distance: 
	 * 0: Sequence-Levenshtein distance
	 * 1: Average number of substitutions
	 * 2: Average number of insertions
	 * 3: Average number of deletions
	 */
	public static double[] distanceDetails(DNASequence a, DNASequence b) {
		
		double[][][] distanceMatrix = LevenshteinBase.computeDetailedDistMatrix(a, b);
		return distanceDetails(distanceMatrix);
	}
	
	/**
	 * Returns an array with detailed information about the Sequence-Levenshtein distance up to a distance of maxDist: 
	 * 0: Sequence-Levenshtein distance
	 * 1: Average number of substitutions
	 * 2: Average number of insertions
	 * 3: Average number of deletions
	 */
	public static double[] inDistanceDetails(double[][] distanceRowCol) {
		int lengthRow = distanceRowCol[0].length / 2;
		int pos = lengthRow;
		double minDist = distanceRowCol[0][pos];
		
		for (int i = 0; i < distanceRowCol[0].length; i++) {
			if (distanceRowCol[0][i] < minDist) {
				pos = i;
				minDist = distanceRowCol[0][pos];
			}
		}
		
		double[] out = {
				distanceRowCol[0][pos],
				distanceRowCol[1][pos],
				distanceRowCol[2][pos],
				distanceRowCol[3][pos]
		};
		
		return out;
	}
	
	/**
	 * Returns an array with detailed information about the Sequence-Levenshtein distance up to a distance of maxDist: 
	 * 0: Sequence-Levenshtein distance
	 * 1: Average number of substitutions
	 * 2: Average number of insertions
	 * 3: Average number of deletions
	 */
	public static double[] inDistanceDetails(DNASequence a, DNASequence b, int maxDist) {
		
		// TODO Use inDistance algorithm to calculate matrix
		//double[][] distanceRowCol = LevenshteinBase.computeDetailedInDistanceResult(a, b);
		//return inDistanceDetails(distanceRowCol);
		return distanceDetails(a, b);
	}

	@Override
	public String getName() {
		return NAME;
	}
	
	@Override
	public String toString() {
		return this.getName();
	}
}
