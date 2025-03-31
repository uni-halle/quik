package de.uni_halle.barcode_calling.util;

import java.lang.StringBuilder;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;

import de.uni_halle.barcode_calling.Main;
import de.uni_halle.barcode_calling.exceptions.UnknownBaseException;
import de.uni_halle.barcode_calling.experiments.data.DNAErrorModel;

public class DNASequence implements Iterable<DNABase>, Comparable<DNASequence> {
	
	/*
	 * Attribute declarations
	 */
	private DNABase[] sequence;
	
	
	/*
	 * Constructors
	 */
	public DNASequence(DNABase[] sequence) {
		this.sequence = sequence;
	}
	
	public DNASequence(String sequence) {
		this.sequence = DNASequence.stringToSequence(sequence);
	}
	
	
	/*
	 * Methods
	 */
	public DNABase[] getSequence() {
		return this.sequence;
	}
	
	public int getLength() {
		return this.sequence.length;
	}
	
	public void setSequence(DNABase[] sequence) {
		this.sequence = sequence;
	}
	
	public DNABase at(int idx) {
		return this.sequence[idx];
	}
	
	public DNASequence corrupt(double p) {
		/*
		 * Returns the result of every base being corrupted with probability p. Mutation, insertion and deletion are equally likely.
		 */
		return this.corrupt(p/3, p/3, p/3);
	}
	
	public DNASequence corrupt(double substitutionProb, double insertionProb, double deletionProb, DNAErrorModel errorModel) {
		/*
		 * Returns the result of every base being corrupted with corresponding probabilities.
		 */
		return DNAErrorModel.corrupt(this, substitutionProb, insertionProb, deletionProb, errorModel);
	}
	
	public DNASequence corrupt(double substitutionProb, double insertionProb, double deletionProb) {
		/*
		 * Returns the result of every base being corrupted with corresponding probabilities. Defaults to the UPHOFF error model.
		 */
		return this.corrupt(substitutionProb, insertionProb, deletionProb, DNAErrorModel.UPHOFF);
	}
	
	public DNASequence corrupt(int mutations, Random generator) {
		/*
		 * Generates DNA sequence with exactly "mutations" corrupted bases determined by the seed of the generator.
		 * No distance to the original sequence can be guaranteed as the corruptions can lead back to the original sequence.
		 */
		// TODO Resolve duplicates with other corrupt method
		double rand;
		int mutationsDone = 0;
		
		int length = this.getLength();
		DNABase[] newSequence = new DNABase[length];
		DNABase[] oldSequence = this.getSequence();
		int readingPosition = 0;
		int writingPosition = 0;
		
		while (readingPosition != length && writingPosition != length) {
			rand = generator.nextDouble()%1;
			
			if (rand <= (double)(mutations-mutationsDone) / (double)(length-writingPosition)) {
				// No corruption
				newSequence[writingPosition] = oldSequence[readingPosition];
				readingPosition++;
				writingPosition++;
			}
			else {
				rand = generator.nextDouble()%1;
				// Corruption
				switch ((int)(3*rand)) {
				case 0:
					// Mutation
					newSequence[writingPosition] = DNABase.randomBaseExcept(oldSequence[readingPosition]);
					readingPosition++;
					writingPosition++;
					break;
				case 1:
					// Insertion
					newSequence[writingPosition] = DNABase.randomBase();
					writingPosition++;
					break;
				case 2:
					// Deletion
					readingPosition++;
					break;
				default:
					throw new RuntimeException("None of the cases fired while corrupting DNA sequence.");
				}
			}
		}
		
		// Fill up with random bases
		for (; writingPosition < length; writingPosition++) {
			newSequence[writingPosition] = DNABase.randomBase();
		}
		
		return new DNASequence(newSequence);
	}
	
	public int sharedPrefixLength(DNASequence other) {
		for (int i = 0; i < this.getLength() && i < other.getLength(); i++) {
			if (this.at(i) != other.at(i)) {
				return i;
			}
		}
		
		return UtilMethods.min(this.getLength(), other.getLength());
	}
	
	public DNASequence subSequence(int start, int end) {
		end = Math.min(this.getLength(), end);  // If seq is too short, return from start to last position
		DNABase[] newSequence = Arrays.copyOfRange(this.sequence, start, end);
		return new DNASequence(newSequence);
	}

	@Override
	public Iterator<DNABase> iterator() {
		return new Iterator<DNABase>() {
			int i = 0;

			@Override
			public boolean hasNext() {
				return i < sequence.length;
			}

			@Override
			public DNABase next() {
				i++;
				return sequence[i-1];
			}
			
		};
	}

	@Override
	public int compareTo(DNASequence other) {
		for (int i = 0; i < this.getLength() && i < other.getLength(); i++) {
			if (this.at(i) != other.at(i)) {
				return this.at(i).ordinal() - other.at(i).ordinal();
			}
		}
		
		return this.getLength() - other.getLength();
	}
	
	@Override
	public String toString() {
		return DNASequence.sequenceToString(this.getSequence());
	}
	
	@Override
	public boolean equals(Object other) {
		if(other instanceof DNASequence) {
			DNASequence otherSequence = (DNASequence) other;
			
            return Arrays.equals(this.sequence, otherSequence.getSequence());
        }
        return false; 
	}
	
	@Override
	public int hashCode() {
		int hash = 0;
		
		for (DNABase b : this) {
			hash = (hash << 2) + b.ordinal();
		}
		
		return hash;
	}
	
	
	/*
	 * Static methods
	 */
	public static DNABase[] stringToSequence(String s) {
		if (s == null) {
			throw new IllegalArgumentException("Emtpy string connot be converted to DNASequence!");
		}
		
		int length = s.length();
		DNABase[] sequence = new DNABase[length];
		int countN = 0;
		char c;
		
		for (int i = 0; i < length; i++) {
			try {
				c = s.charAt(i);
				sequence[i - countN] = DNABase.fromChar(c);
			}
			catch (UnknownBaseException e) {
				countN++;
			}
		}
		
		if (countN == length) {
			throw new IllegalArgumentException("Read contains only Ns!");
		}
		if (countN > 0) {
			sequence = Arrays.copyOfRange(sequence, 0, length - countN);
			Main.writeToLog("Notice: " + countN + " bases of type 'N' deleted.");
		}
		
		return sequence;
	}
	
	public static String sequenceToString(DNABase[] sequence) {
		int length = sequence.length;
		StringBuilder s = new StringBuilder(length);
		char c;
		
		for (int i = 0; i < length; i++) {
			c = sequence[i].toChar();
			s.append(c);
		}
		
		return s.toString();
	}
	
	public static DNASequence random(int length, Random generator) {
		DNABase[] seq = new DNABase[length];
		
		for (int i = 0; i < length; i++) {
			seq[i] = DNABase.randomBase(generator);
		}
		
		return new DNASequence(seq);
	}
	
	public static DNASequence random(int length) {
		return random(length, new Random());
	}
	
}
