package BedFileValidation.Objects;

import PersonalProteome.Definitions;

public class BedFileLine {

	
	int id = 0;
	int chromosomeName = -1;
	int startLocation = -1;
	int stopLocation = -1;
	String sequence = "";
	String score = "";
	int strand = -1;
	String color = "";
	String blockCount = "";
	String blockSize = "";
	String blockStart = "";
	

	String dnaSeq = "";
	/**
	 * @param id
	 * @param chromosomeName
	 * @param startLocation
	 * @param stopLocation
	 * @param sequence
	 * @param score
	 * @param strand
	 * @param blockCount
	 * @param blockSize
	 * @param blockStart
	 */
	public BedFileLine(int id, int chromosomeName, int startLocation,
			int stopLocation, String sequence, String score, int strand, String color,
			String blockCount, String blockSize, String blockStart) {
		super();
		this.id = id;
		this.chromosomeName = chromosomeName;
		this.startLocation = startLocation;
		this.stopLocation = stopLocation;
		this.sequence = sequence;
		this.score = score;
		this.strand = strand;
		this.color = color;
		this.blockCount = blockCount;
		this.blockSize = blockSize;
		this.blockStart = blockStart;
	}
	
	
	/** 
	 * Returns a valid bed file line representing this object.
	 */
	public String toString() {
		
		String strand = "+";
		if(this.strand == Definitions.genomicStrandNEGATIVE){
			strand = "-";
		}
		
		return  Definitions.convertChrmNumToString(chromosomeName) + "\t"
				+ startLocation + "\t"
				+ stopLocation + "\t"
				+ sequence + "\t"
				+ score + "\t"
				+ strand + "\t"
				+ startLocation + "\t"
				+ stopLocation + "\t"
				+ color + "\t"
				+ blockCount + "\t"
			    + blockSize + "\t"
				+ blockStart + "\t";
	}


	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}
	/**
	 * @return the chromosomeName
	 */
	public int getChromosomeName() {
		return chromosomeName;
	}
	/**
	 * @return the startLocation
	 */
	public int getStartLocation() {
		return startLocation;
	}
	/**
	 * @return the stopLocation
	 */
	public int getStopLocation() {
		return stopLocation;
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}
	/**
	 * @return the score
	 */
	public String getScore() {
		return score;
	}
	/**
	 * @return the strand
	 */
	public int getStrand() {
		return strand;
	}
	/**
	 * @return the blockCount
	 */
	public String getBlockCount() {
		return blockCount;
	}
	/**
	 * @return the blockSize
	 */
	public String getBlockSize() {
		return blockSize;
	}
	/**
	 * @return the blockStart
	 */
	public String getBlockStart() {
		return blockStart;
	}
	
	/**
	 * @return the color
	 */
	public String getColor() {
		return color;
	}


	/**
	 * @return the dnaSeq
	 */
	public String getDnaSeq() {
		return dnaSeq;
	}


	/**
	 * @param dnaSeq the dnaSeq to set
	 */
	public void setDnaSeq(String dnaSeq) {
		this.dnaSeq = dnaSeq;
	}
	
}
