package Peppy;

import Math.HasValue;


/**
 * Is a data class that stores:
 * 1) a sequence of amino acids.
 * 2) the theoretical mass
 * 3) the begin index (end index can be calculated)
 * 4) boolean forward (false = reverse)
 * 5) the reading frame
 * @author Brian Risk
 *
 */
public class Peptide implements Comparable<Peptide>, HasValue {
	
//	private String acidSequence;
	private byte [] acidSequence;
	private double mass;
	private int startIndex;
	private int stopIndex;
	private int intronStartIndex;
	private int intronStopIndex;
	private boolean forward;
	private Sequence_DNA parentSequence;
	private Protein protein;
	private boolean isSpliced;
	private boolean isMatched = false;
	private int lengthMinusOne;
	private int cleavageAcidCount;
	
	
	public boolean isMatched() {
		return isMatched;
	}

	public void setMatched(boolean isMatched) {
		this.isMatched = isMatched;
	}

	/**
	 * just gets an amino acid sequence.
	 * @param sequence
	 */
	public Peptide(String sequence) {
		this(sequence, 0, sequence.length(), -1, -1, true, null, null, false);
	}
	
	public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence_DNA parentSequence, boolean isSpliced) {
		this(acidSequence, startIndex, stopIndex, intronStartIndex, intronStopIndex, forward, parentSequence, null, isSpliced);
	}
	
	public Peptide(String acidSequence, int startIndex, int stopIndex, int intronStartIndex, int intronStopIndex, boolean forward, Sequence_DNA parentSequence, Protein protein, boolean isSpliced) {
		this.acidSequence = AminoAcids.getByteArrayForString(acidSequence);
		this.mass = calculateMass();
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
		this.intronStartIndex = intronStartIndex;
		this.intronStopIndex = intronStopIndex;
		this.forward = forward;
		this.parentSequence = parentSequence;
		this.protein = protein;
		this.isSpliced = isSpliced;
		this.lengthMinusOne = this.acidSequence.length - 1;
		
		cleavageAcidCount = 0;
		for (int i = 0; i < this.acidSequence.length; i++) {
			if (this.acidSequence[i] == AminoAcids.K || this.acidSequence[i] == AminoAcids.R) {
				cleavageAcidCount++;
			}
		}
	}


	public int getLengthMinusOne() {
		return lengthMinusOne;
	}

	@Override
	public String toString() {
//		return mass + "\t" + getAcidSequenceString() + "\t" + startIndex + "\t" + proteinName;
		return  getAcidSequenceString() + "\t" + getMass() + "\t" + getStartIndex() + "\t" +  getStopIndex()  + "\t" + forward;
//		return  getAcidSequenceString();
	}


	public int compareTo(Peptide other) {
		if (mass > other.getMass()) return  1;
		if (mass < other.getMass()) return -1;
		return 0;
	}
	
	/**
	 * Okay, this equals is not in line with the way things work for compareTo.
	 * this compares acid sequences for equality.  compareTo compares masses.
	 * 
	 * the real trick for equality is ignoring any trailing stop (".") codon
	 */
	public boolean equals(byte [] otherAcidSequence) {
		int ourLength = acidSequence.length;
		int theirLength = otherAcidSequence.length;
		
		//ignoring terminating STOPs
		if (acidSequence[acidSequence.length - 1] == AminoAcids.STOP) ourLength--;
		if (otherAcidSequence[otherAcidSequence.length - 1] == AminoAcids.STOP) theirLength--;
		
		//if peptids not same length, they are not equal
		if (ourLength != theirLength) return false;
		
		//compare each acid
		for (int i = 0; i < ourLength; i++) {
			if (acidSequence[i] != otherAcidSequence[i]) return false;
		}
		
		//if reached this point, they are equal
		return true;
	}
	
	/**
	 * Equals if every sequential acid weighs the same as that of the other sequence
	 */
	public boolean equalsByAcidMasses(byte [] otherAcidSequence) {
		if (!equals(otherAcidSequence)) return false;
		
		for (int i = 0; i < acidSequence.length; i++) {
			if (AminoAcids.getWeightMono(acidSequence[i]) != AminoAcids.getWeightMono(otherAcidSequence[i])) return false;
		}
		
		return true;
			
	}
	
	public boolean equals(Peptide peptide) {
		if (mass == peptide.getMass()) {
			return equals(peptide.getAcidSequence());
		} else {
			return false;
		}
	}
	
	public boolean equals(String acidSequenceString) {
		return equals(AminoAcids.getByteArrayForString(acidSequenceString));
	}


	/**
	 * @return the sequence
	 */
	public byte [] getAcidSequence() {
		return acidSequence;
	}
	
	public String getAcidSequenceString() {
		return AminoAcids.getStringForByteArray(acidSequence);
	}
	
	public int getLength() {
		return acidSequence.length;
	}


	/**
	 * @return the mass
	 */
	public double getMass() {
		return mass;
	}


	/**
	 * @return the index
	 */
	public int getStartIndex() {
		if (forward) {
			return startIndex;
		} else {
			return stopIndex + 1;
		}
	}
	
	public int getStopIndex() {
		if (forward) {
			return stopIndex;
		} else {
			return startIndex + 1;
		}
	}

	public int getIntronStartIndex() {
		if (forward) {
			return intronStartIndex;
		} else {
			return intronStopIndex + 1;
		}
	}

	public int getIntronStopIndex() {
		if (forward) {
			return intronStopIndex;
		} else {
			return intronStartIndex + 1;
		}
	}

	public Protein getProtein() {
		return protein;
	}


	public int getCleavageAcidCount() {
		return cleavageAcidCount;
	}


	/**
	 * @return the forward
	 */
	public boolean isForward() {
		return forward;
	}


	public boolean isSpliced() {
		return isSpliced;
	}

	/**
	 * @return the parentSequence
	 */
	public Sequence_DNA getParentSequence() {
		return parentSequence;
	}
	
	/**
	 * This will calculate either the mono mass or the average mass depending on the setting
	 * in your Properties object.
	 * @return
	 */
	private double calculateMass() {
		double mass = 0.0;
		if (Properties.useMonoMass) {
			for (int i = 0; i < acidSequence.length; i++) {
				if (AminoAcids.isValid(acidSequence[i])) {
					mass += AminoAcids.getWeightMono(acidSequence[i]);
				} else {
					mass = -1;
					return mass;
				}
			}
			mass += Definitions.WATER_MONO;
		} else {
			for (int i = 0; i < acidSequence.length; i++) {
				if (AminoAcids.isValid(acidSequence[i])) {
					mass += AminoAcids.getWeightAverage(acidSequence[i]);
				} else {
					mass = -1;
					return mass;
				}
			}
			mass += Definitions.WATER_AVERAGE;
		}
		return mass;
	}

	public double getValue() {
		return getMass();
	}
	
	public double getHydrophobicProportion() {
		double out = 0;
		for (int i = 0; i < acidSequence.length; i++) {
			/* (leucine, valine, isoleucine, phenylalanine, methionine, cysteine and tryptophan */
			if (acidSequence[i] == AminoAcids.G) out++;
			if (acidSequence[i] == AminoAcids.A) out++;
			if (acidSequence[i] == AminoAcids.V) out++;
			if (acidSequence[i] == AminoAcids.L) out++;
			if (acidSequence[i] == AminoAcids.I) out++;
			if (acidSequence[i] == AminoAcids.M) out++;
			if (acidSequence[i] == AminoAcids.F) out++;
			if (acidSequence[i] == AminoAcids.W) out++;
			if (acidSequence[i] == AminoAcids.P) out++;
		}
		out /= acidSequence.length;
		return out;
	}
	
	public double getHydrophillicProportion() {
		double out = 0;
		for (int i = 0; i < acidSequence.length; i++) {
			/* (leucine, valine, isoleucine, phenylalanine, methionine, cysteine and tryptophan */
			if (acidSequence[i] == AminoAcids.S) out++;
			if (acidSequence[i] == AminoAcids.T) out++;
			if (acidSequence[i] == AminoAcids.C) out++;
			if (acidSequence[i] == AminoAcids.Y) out++;
			if (acidSequence[i] == AminoAcids.N) out++;
			if (acidSequence[i] == AminoAcids.Q) out++;
		}
		out /= acidSequence.length;
		return out;
	}
	

}
