package de.uni_halle.barcode_calling.util;

public abstract class Tuples {
	
	public static class IdDoubleTuple implements Comparable<IdDoubleTuple> {
		public int id;
		public double d;
		
		public IdDoubleTuple(int i, double d) {
			this.id = i;
			this.d = d;
		}
		
		@Override
		public String toString() {
			return "(" + this.id + ", " + this.d + ")";
		}

		@Override
		public int compareTo(IdDoubleTuple other) {
			return (int)Math.signum(this.d - other.d);
		}
	}
	
	public static class IdIntegerTuple implements Comparable<IdIntegerTuple> {
		public int id;
		public long i;
		
		public IdIntegerTuple(int i1, long i2) {
			this.id = i1;
			this.i = i2;
		}
		
		@Override
		public String toString() {
			return "(" + this.id + ", " + this.i + ")";
		}
		
		@Override
		public int compareTo(IdIntegerTuple other) {
			return (int)Math.signum(this.i - other.i);	// TODO Test without signum
		}
	}
	
	public static class SequenceIndexTuple implements Comparable<SequenceIndexTuple> {
		public DNASequence seq;
		public Integer idx;
		
		public SequenceIndexTuple(DNASequence seq, Integer idx) {
			this.seq = seq;
			this.idx = idx;
		}
		
		@Override
		public int compareTo(SequenceIndexTuple other) {
			return this.seq.compareTo(other.seq);
		}
	}

}

