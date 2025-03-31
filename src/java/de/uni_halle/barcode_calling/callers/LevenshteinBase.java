package de.uni_halle.barcode_calling.callers;

import static de.uni_halle.barcode_calling.util.UtilMethods.max;
import static de.uni_halle.barcode_calling.util.UtilMethods.min;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;

import de.uni_halle.barcode_calling.util.DNASequence;
import de.uni_halle.barcode_calling.util.UtilMethods;

/**
 * Base for Levenshtein baselines
 * @author Riko Uphoff
 *
 */
public class LevenshteinBase implements DistanceMeasureBase {
	
	public static final String NAME = "Levenshtein";
	
	public int getDefaultMaxDistance(int length) {
		// TODO This has not been tested thoroughly
		return (int)Math.round(length / 4.0);
	}
	
	/**
	 * Classical dynamic program for computing Levensthein distance.
	 * @param a Sequence A
	 * @param b Sequence B
	 * @return Levenshtein distance between A and B.
	 */
	public int distance(DNASequence a, DNASequence b) {
		return computeDistMatrix(a, b)[a.getLength()][b.getLength()];
	}
	
	/**
	 * Classical dynamic program for computing Levensthein distance.
	 * @param a Sequence A
	 * @param b Sequence B
	 * @return Distance matrix of the dynamic program.
	 */
	public static int[][] computeDistMatrix(DNASequence a, DNASequence b) {
		int lengthA = a.getLength();
		int lengthB = b.getLength();
		int editCost;
		
		// Distance matrix
		int[][] dm = new int[lengthA + 1][lengthB + 1];
		
		
		for (int i = 0; i <= lengthA; i++) {
			for (int j = 0; j <= lengthB; j++) {
		
				if (i == 0) {
					dm[i][j] = j;
				}
				else if (j == 0) {
					dm[i][j] = i;
				}
				else {
					editCost = a.at(i-1) == b.at(j-1) ? 0 : 1;
					dm[i][j] = UtilMethods.min(
							dm[i - 1][j - 1] + editCost, // replace
							dm[i - 1][j] + 1, // delete
							dm[i][j - 1] + 1 // insert
					);
				}
			}
		}
	
		return dm;
	}
	
	/**
	 * Implements the thresholded Levenshtein calculation after Zorita et al. (2015)
	 * @param a Sequence A
	 * @param b Sequence B
	 * @param maxDist Maximum computable distance
	 * @return Are A and B within a Levenshtein distance of maxDist?
	 */
	public boolean inDistance(DNASequence a, DNASequence b, int maxDist) {
		return inDistanceResult(a, b, maxDist) <= maxDist;
	}
	
	/**
	 * Implements the thresholded Levenshtein calculation after Zorita et al. (2015)
	 * @param a Sequence A
	 * @param b Sequence B
	 * @param maxDist Maximum computable distance
	 * @return Levenshtein distance between a and b up to maxDist
	 */
	public int inDistanceResult(DNASequence a, DNASequence b, int maxDist) {
		// Account for negative maxDist
		if (maxDist < 0) return maxDist+1;
		
		// Use shorter string as reference
		if (a.getLength() > b.getLength()) {
			DNASequence temp = a;
			a = b;
			b = temp;
		}
		
		int minLength = a.getLength();
		int lengthDif = b.getLength() - a.getLength();
		int posA, posB, editCost;
		int columnSize = maxDist+1 > lengthDif ? maxDist+1 - lengthDif : 0;
		int rowSize = maxDist+1 + lengthDif;
		int[][] columns = new int[minLength+1][columnSize];	// Column is indexed top to bottom
		int[][] rows = new int[minLength+1][rowSize];		// Row is indexed left to right
		int[] centers = new int[minLength+1];
		
		centers[0] = lengthDif;
		
		
		for (int depth = 0; depth <= minLength; depth++) {
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
			
			// Early stopping
			if (centers[depth] > maxDist) {
				return maxDist+1;
			}
		}
		
		return centers[minLength];
	}
	
	/**
	 * Computes the average frequencies of the error types according to the procedure in Section 3.2, Figure 3
	 * @param a Sequence A
	 * @param b Sequence B
	 * @return Four matrices containing the average number of errors for each error type:
	 * [0]: Matrix containing total distance
	 * [1]: Matrix containing average number of subs.
	 * [2]: Matrix containing average number of ins.
	 * [3]: Matrix containing average number of del.
	 */
	public static double[][][] computeDetailedDistMatrix(DNASequence a, DNASequence b) {
		int lengthA = a.getLength();
		int lengthB = b.getLength();
		int editCost;
		
		// Deatiled distance matrix
		/*
		 * [0]: Matrix containing total distance
		 * [1]: Matrix containing average number of subs.
		 * [2]: Matrix containing average number of ins.
		 * [3]: Matrix containing average number of del.
		 */
		double[][][] dm = new double[4][lengthA + 1][lengthB + 1];
		int[][] numPaths = new int[lengthA + 1][lengthB + 1];
		
		for (int i = 0; i <= lengthA; i++) {
			for (int j = 0; j <= lengthB; j++) {
		
				if (i == 0 && j == 0) {
					dm[0][i][j] = 0;
					numPaths[i][j] = 1;
				}
				else if (i == 0) {
					dm[0][i][j] = j;
					dm[1][i][j] = dm[1][i][j-1];
					dm[2][i][j] = dm[2][i][j-1] + 1;
					dm[3][i][j] = dm[3][i][j-1];
					numPaths[i][j] = numPaths[i][j-1];
				}
				else if (j == 0) {
					dm[0][i][j] = i;
					dm[1][i][j] = dm[1][i-1][j];
					dm[2][i][j] = dm[2][i-1][j];
					dm[3][i][j] = dm[3][i-1][j] + 1;
					numPaths[i][j] = numPaths[i-1][j];
				}
				else {
					editCost = a.at(i-1) == b.at(j-1) ? 0 : 1;
					dm[0][i][j] = UtilMethods.min(
							(int)dm[0][i - 1][j - 1] + editCost, // substitute
							(int)dm[0][i][j - 1] + 1, // insert
							(int)dm[0][i - 1][j] + 1 // delete
					);
					
					if (dm[0][i][j] == dm[0][i - 1][j - 1] + editCost) {
						dm[1][i][j] += numPaths[i-1][j-1] * (dm[1][i-1][j-1]);
						dm[2][i][j] += numPaths[i-1][j-1] * (dm[2][i-1][j-1]);
						dm[3][i][j] += numPaths[i-1][j-1] * (dm[3][i-1][j-1]);
						
						if (editCost == 1) {
							// Substitution
							dm[1][i][j] += numPaths[i-1][j-1];
						}
						
						numPaths[i][j] += numPaths[i-1][j-1];
					}
					if (dm[0][i][j] == dm[0][i][j - 1] + 1) {
						// Insertion
						dm[1][i][j] += numPaths[i][j-1] * (dm[1][i][j-1]);
						dm[2][i][j] += numPaths[i][j-1] * (dm[2][i][j-1] + 1);
						dm[3][i][j] += numPaths[i][j-1] * (dm[3][i][j-1]);
						numPaths[i][j] += numPaths[i][j-1];
					}
					if (dm[0][i][j] == dm[0][i - 1][j] + 1) {
						// Deletion
						dm[1][i][j] += numPaths[i-1][j] * (dm[1][i-1][j]);
						dm[2][i][j] += numPaths[i-1][j] * (dm[2][i-1][j]);
						dm[3][i][j] += numPaths[i-1][j] * (dm[3][i-1][j] + 1);
						numPaths[i][j] += numPaths[i-1][j];
					}
					
					for (int e = 1; e < 4; e++) {
						dm[e][i][j] /= (double)numPaths[i][j];
					}
				}
			}
		}
	
		return dm;
	}
	
	/**
	 * Returns an array with detailed information about the Levenshtein distance: 
	 * 0: Levenshtein distance
	 * 1: Average number of substitutions
	 * 2: Average number of insertions
	 * 3: Average number of deletions
	 */
	public static double[] distanceDetails(double[][][] detailedDistMatrix) {
		
		int posA = detailedDistMatrix[0].length - 1;
		int posB = detailedDistMatrix[0][0].length - 1;
		double[] out = {
				detailedDistMatrix[0][posA][posB],
				detailedDistMatrix[1][posA][posB],
				detailedDistMatrix[2][posA][posB],
				detailedDistMatrix[3][posA][posB]
		};
		
		return out;
	}
	
	public static double[][] computeDetailedInDistanceResult(DNASequence a, DNASequence b, int maxDist) {
		/*TODO Implement with thresholded algorithm from Zorita et al. (2015)
		// Account for negative maxDist
		if (maxDist < 0) return null;
		
		// Use shorter string as reference
		boolean sequencesSwitched = false;
		
		if (a.getLength() > b.getLength()) {
			DNASequence temp = a;
			a = b;
			b = temp;
			sequencesSwitched = true;
		}
		
		int minLength = a.getLength();
		int lengthDif = b.getLength() - a.getLength();
		int posA, posB, editCost;
		int columnSize = maxDist+1 > lengthDif ? maxDist+1 - lengthDif : 0;
		int rowSize = maxDist+1 + lengthDif;
		/* 
		 * 0: Levenshtein distance
		 * 1: Average number of substitutions
		 * 2: Average number of insertions
		 * 3: Average number of deletions
		 * 4: Number of paths
		*/
		/*TODO
		double[][][] columns = new double[5][minLength+1][columnSize];	// Column is indexed top to bottom
		double[][][] rows = new double[5][minLength+1][rowSize];		// Row is indexed left to right
		double[][] centers = new double[5][minLength+1];
		
		centers[0][0] = lengthDif;
		
		
		for (int depth = 0; depth <= minLength; depth++) {
			// Column
			for (int i = 0; i < columnSize; i++) {
				if (i < max(1, columnSize+1 - depth) || depth == 0) {
					// Border
					columns[0][depth][i] = Math.abs(maxDist+1 - i);		// Distance
					columns[1][depth][i] = 0;							// Substitutions
					columns[2][depth][i] = Math.abs(maxDist+1 - i);		// Insertions
					columns[3][depth][i] = 0;							// Deletions
					columns[4][depth][i] = 1;							// Number of paths
				}
				else {
					// Inside
					posA = depth-1 - columnSize + i;
					posB = depth-1 + lengthDif;
					
					editCost = a.at(posA) == b.at(posB) ? 0 : 1;
					double leftEntry = i == columnSize-1 ? centers[0][depth-1] + 1 : columns[0][depth-1][i+1] + 1;
					
					columns[0][depth][i] = min(
							(int)columns[0][depth-1][i] + editCost,			// Substitution
							(int)columns[0][depth][i-1] + 1, 				// Insertion
							(int)leftEntry									// Deletion
					);
					
					if (columns[0][depth][i] == columns[0][depth-1][i] + editCost) {
						columns[1][depth][i] += columns[5][depth-1][i] * columns[1][depth-1][i];
						columns[2][depth][i] += columns[5][depth-1][i] * columns[2][depth-1][i];
						columns[3][depth][i] += columns[5][depth-1][i] * columns[3][depth-1][i];
						
						if (editCost == 1) {
							// Substitution
							columns[1][depth][i] += columns[5][depth-1][i];
						}
						
						columns[5][depth][i] += columns[5][depth-1][i];
					}
					if (columns[0][depth][i] == columns[0][depth][i-1] + 1) {
						// Insertion
						dm[1][i][j] += numPaths[i][j-1] * (dm[1][i][j-1]);
						dm[2][i][j] += numPaths[i][j-1] * (dm[2][i][j-1] + 1);
						dm[3][i][j] += numPaths[i][j-1] * (dm[3][i][j-1]);
						numPaths[i][j] += numPaths[i][j-1];
					}
					if (dm[0][i][j] == dm[0][i - 1][j] + 1) {
						// Deletion
						dm[1][i][j] += numPaths[i-1][j] * (dm[1][i-1][j]);
						dm[2][i][j] += numPaths[i-1][j] * (dm[2][i-1][j]);
						dm[3][i][j] += numPaths[i-1][j] * (dm[3][i-1][j] + 1);
						numPaths[i][j] += numPaths[i-1][j];
					}
					
					for (int e = 1; e < 4; e++) {
						dm[e][i][j] /= (double)numPaths[i][j];
					}
				}
			}
			
			// Row
			for (int i = 0; i < rowSize; i++) {
				if (i < max(1, rowSize+1 - depth - lengthDif) || depth == 0) {
					// Border
					rows[0][depth][i] = Math.abs(maxDist+1 - i);		// Distance
					rows[1][depth][i] = 0;								// Substitutions
					rows[2][depth][i] = 0;								// Insertions
					rows[3][depth][i] = Math.abs(maxDist+1 - i);		// Deletions
					rows[4][depth][i] = 1;								// Number of paths
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
			
			// Early stopping
			if (centers[depth] > maxDist) {
				return maxDist+1;
			}
		}
		
		return centers[minLength];
		TODO*/
		return null;
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
