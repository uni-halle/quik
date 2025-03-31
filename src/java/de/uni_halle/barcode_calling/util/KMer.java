package de.uni_halle.barcode_calling.util;

import java.util.Iterator;

public class KMer implements Iterable<DNABase> {
	/**
	 * Class to efficiently iterate the k-mers of a sequence.
	 */
	
	int k;
	DNABase[] mer;
	int pos;
	int hash;
	int maxHash;
	
	public KMer(DNABase[] kmer) {
		this.k = kmer.length;
		this.mer = kmer;
		this.pos = 0;
		this.hash = this.calcHashCode();
		this.maxHash = 1 << (k << 1);	// 2^(2k)
	}
	
	public KMer(int k) {
		this.k = k;
		this.mer = new DNABase[k];
		this.pos = -k;
		this.hash = 0;
		this.maxHash = 1 << (k << 1);	// 2^(2k)
	}
	
	public KMer(DNASequence kmer) {
		this(kmer.getSequence());
	}
	
	public int getPosition() {
		return this.pos;
	}
	
	public void append(DNABase b) {
		this.mer[(this.pos+k) % k] = b;
		pos++;
		this.hash = ((this.hash << 2) + b.ordinal()) % this.maxHash;
	}
	
	public DNASequence toSequence() {
		if (this.pos >= 0) {
			DNABase[] seq = new DNABase[this.k];
			int i = 0;
			
			for (DNABase b : this) {
				seq[i] = b;
				i++;
			}
			
			return new DNASequence(seq);
		}
		else {
			return null;
		}
	}
	
	public static Iterable<KMer> iterateKMers(DNASequence seq, int k) {
		return new Iterable<KMer>() {

			@Override
			public Iterator<KMer> iterator() {
				return new Iterator<KMer>() {
					
					KMer kmer;
					int i = 0;

					@Override
					public boolean hasNext() {
						return i < seq.getLength() - (k-1);
					}

					@Override
					public KMer next() {
						if (i == 0) {
							kmer = new KMer(seq.subSequence(0, k));
						}
						else {
							kmer.append(seq.at(i+k-1));
						}
						
						i++;
						return kmer;
					}
					
				};
			}
			
		};
	}
	
	@Override
	public int hashCode() {
		return this.hash;
	}
	
	public int calcHashCode() {
		int hash = 0;
		
		for (DNABase b : this) {
			hash = (hash << 2) + b.ordinal();
		}
		
		return hash;
	}

	@Override
	public Iterator<DNABase> iterator() {
		return new Iterator<DNABase>() {
			
			int i = 0;
			int curPos = pos;

			@Override
			public boolean hasNext() {
				return i < k && curPos >= 0;
			}

			@Override
			public DNABase next() {
				DNABase out = mer[curPos % k];
				i++;
				curPos++;
				return out;
			}
			
		};
	}
	
	@Override
	public String toString() {
		StringBuilder s = new StringBuilder(k);
		
		for (DNABase b : this) {
			s.append(b);
		}
		
		return s.toString();
	}
	
}
