package de.uni_halle.barcode_calling.util;

import java.lang.Math;
import java.util.Random;

import de.uni_halle.barcode_calling.exceptions.UnknownBaseException;

public enum DNABase {
	A, T, G, C;
	
	public char toChar() {
		return this.toString().charAt(0);
	}
	
	public static DNABase fromChar(char c) {
		// TODO Handle Ns
		switch (c) {
		case 'A':
			return DNABase.A;
		case 'T':
			return DNABase.T;
		case 'G':
			return DNABase.G;
		case 'C':
			return DNABase.C;
		case 'N':
			throw new UnknownBaseException("N-base detected!");
		default:
			throw new IllegalArgumentException("Invalid base: " + c);
		}
	}
	
	public static DNABase fromOrdinal(int ordinal) {
		switch (ordinal) {
		case 0:
			return DNABase.A;
		case 1:
			return DNABase.T;
		case 2:
			return DNABase.G;
		case 3:
			return DNABase.C;
		default:
			throw new RuntimeException("None of the cases fired while generating random base " + ordinal + ".");
		}
	}
	
	public static DNABase randomBase() {
		double rand = Math.random();
		
		return fromOrdinal((int)(rand / 0.25));
	}
	
	public static DNABase randomBase(Random generator) {
		int rand = Math.abs(generator.nextInt());
		
		return fromOrdinal(rand % 4);
	}
	
	public static DNABase randomBaseExcept(DNABase exception) {
		DNABase out;
		do {
			out = DNABase.randomBase();
		}
		while (out == exception);
		
		return out;
	}
}
