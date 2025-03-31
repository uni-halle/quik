package de.uni_halle.barcode_calling.experiments.data;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import de.uni_halle.barcode_calling.util.DNABase;
import de.uni_halle.barcode_calling.util.DNASequence;

public enum DNAErrorModel {
	UPHOFF, PRESS;
	
	public static DNASequence corrupt(DNASequence seq, double substitutionProb, double insertionProb, double deletionProb, DNAErrorModel model) {
		switch (model) {
			case UPHOFF:
				return corruptUphoff(seq, substitutionProb, insertionProb, deletionProb);
			case PRESS:
				return corruptPress(seq, substitutionProb, insertionProb, deletionProb);
			default:
				return null;	// TODO Throw exception
		}
	}
	
	public static DNASequence corrupt(DNASequence seq, double p, DNAErrorModel model) {
		return corrupt(seq, p/3, p/3, p/3, model);
	}
	
	private static DNASequence corruptUphoff(DNASequence seq, double substitutionProb, double insertionProb, double deletionProb) {
		/*
		 * Returns the result of every base being corrupted with corresponding probabilities according to Uphoff et. al.
		 */
		// TODO Memorize when insertion and deletion happened and prevent insertion from being deleted again (and vice versa)
		int length = seq.getLength();
		DNABase[] newSequence = new DNABase[length];
		DNABase[] oldSequence = seq.getSequence();
		int readingPosition = 0;
		int writingPosition = 0;
		int justInserted = 0;
		int justDeleted = 0;
		double rand;
		
		while (readingPosition != length && writingPosition != length) {
			rand = Math.random();
			
			if (rand <= substitutionProb) {
				// Mutation
				newSequence[writingPosition] = DNABase.randomBaseExcept(oldSequence[readingPosition]);
				readingPosition++;
				writingPosition++;
				
				justInserted = 0;
				justDeleted = 0;
			}
			else if (rand <= substitutionProb + (1-justDeleted) * insertionProb + (justInserted) * deletionProb) {
				// Insertion
				newSequence[writingPosition] = DNABase.randomBase();
				writingPosition++;
				
				justInserted = 1;
				justDeleted = 0;
			}
			else if (rand <= substitutionProb + insertionProb + deletionProb) {
				// Deletion
				readingPosition++;
				
				justInserted = 0;
				justDeleted = 1;
			}
			else {
				// No corruption
				newSequence[writingPosition] = oldSequence[readingPosition];
				readingPosition++;
				writingPosition++;
				
				justInserted = 0;
				justDeleted = 0;
			}
		}
		
		// Fill up with random bases
		for (; writingPosition < length; writingPosition++) {
			newSequence[writingPosition] = DNABase.randomBase();
		}
		
		return new DNASequence(newSequence);
	}
	
	private static DNASequence corruptPress(DNASequence seq, double substitutionProb, double insertionProb, double deletionProb) {
		/*
		 * Returns the result of every base being corrupted with corresponding probabilities according to Press
		 */
		Random generator = new Random();
		int length = seq.getLength();
		DNABase[] out = seq.getSequence();
		int randPos;
		
		// Substitutions
		int numSub = getBinomial(length, substitutionProb * 4/3);	// TODO What if substitutionProb > 3/4?
		Set<Integer> posSub = new HashSet<Integer>();
		// Generate positions
		for (int i = 0; i < numSub; i++) {
			randPos = generator.nextInt(length);
			posSub.add(randPos);
		}
		// Generate Bases
		for (int pos : posSub) {
			out[pos] = DNABase.randomBase();
		}
		
		// Deletions
		int numDel = getBinomial(length, deletionProb);
		Set<Integer> posDel = new HashSet<Integer>();
		// Generate positions
		for (int i = 0; i < numDel; i++) {
			randPos = generator.nextInt(length);
			posDel.add(randPos);
		}
		// Delete Bases
		DNABase[] newOut = new DNABase[length - posDel.size()];
		int posNewOut = 0;
				
		for (int i = 0; i < length; i++) {
			if (!posDel.contains(i)) {
				newOut[posNewOut] = out[i];
				posNewOut++;
			}
		}
		
		out = newOut;
		length = newOut.length;
		
		// Insertions
		int numIns = getBinomial(length, insertionProb);
		Map<Integer, Integer> posIns = new HashMap<Integer, Integer>();
		// Generate positions
		for (int i = 0; i < numIns; i++) {
			randPos = generator.nextInt(length+1);
			if (posIns.containsKey(randPos)) {
				posIns.replace(randPos, posIns.get(randPos)+1);
			}
			else {
				posIns.put(randPos, 1);
			}
		}
		// Insert Bases
		newOut = new DNABase[length + numIns];
		posNewOut = 0;
				
		for (int i = 0; i <= length; i++) {
			if (posIns.containsKey(i)) {
				for (int j = 0; j < posIns.get(i); j++) {
					newOut[posNewOut] = DNABase.randomBase();
					posNewOut++;
				}
			}
			if (i < length) {
				newOut[posNewOut] = out[i];
				posNewOut++;
			}
		}
		
		out = newOut;
		length = newOut.length;
		
		// Pad/Truncate
		newOut = new DNABase[seq.getLength()];
		for (int i = 0; i < newOut.length; i++) {
			if (i < out.length) {
				newOut[i] = out[i];
			}
			else {
				newOut[i] = DNABase.randomBase();
			}
		}
		
		return new DNASequence(newOut);
	}
	
	public static int getBinomial(int n, double p) {
		double log_q = Math.log(1.0 - p);
	    int x = 0;
	    double sum = 0;
	    for(;;) {
	       sum += Math.log(Math.random()) / (n - x);
	       if(sum < log_q) {
	          return x;
	       }
	       x++;
	    }
	}
		
}
