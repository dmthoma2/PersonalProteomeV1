package PersonalProteome.PeptideAnalysisTool.SubClasses;

import PersonalProteome.Definitions;


/**
 * PeptideLocationLine represents a single line from a peptide file.  This is used when labeling peptides based on their location.  It preserves informaton from the input file, and its toString()
 * method returns this peptide just as it was read it.  It additionally knows information about its location type from the definitions file.  This is used to determine which type of peptide this is
 * and what file to output it to.
 * @author David "Corvette" Thomas
 *
 */
public class PeptideLocationLine implements Comparable<PeptideLocationLine> {
	//Variables to store information as read in from the input file
	private int chromosomeName = -1;
	private int startLocation = -1;
	private int stopLocation = -1;
	private String sequence = "";
	private String score = "";
	private int strand = -1;
	private String restOfLine = "";
//	private String score3 = "";
//	private String score4 = "";
//	private String score5 = "";
	
	//Variables to store information about this peptides location.
	private int startLocType = Definitions.LOCATION_TYPE_UNKNOWN;
	private int stopLocType = Definitions.LOCATION_TYPE_UNKNOWN;
	private boolean isSameSubUnit = true;

	/**
	 * Peptide Location Line represents a single line from a peptide location input file.  It preserves all of the data from the peptides input.
	 * @param chromosomeId
	 * @param startLocation
	 * @param stopLocation
	 * @param sequence
	 * @param score
	 * @param strand
	 * @param restOfLine
	 */
	public PeptideLocationLine(int chromosomeName, int startLocation,
			int stopLocation, String sequence, String score, int strand,
			String restOfLine) {
		super();
		this.chromosomeName = chromosomeName;
		this.startLocation = startLocation;
		this.stopLocation = stopLocation;
		this.sequence = sequence;
		this.score = score;
		this.strand = strand;
		this.restOfLine = restOfLine;
//		this.score3 = score3;
//		this.score4 = score4;
//		this.score5 = score5;
	}

	/**
	 * toString() returns this PeptideLocationLine as a string representation as it would appear in the input file.  
	 */
	public String toString() {
		String chrmStr = "";
		String strandStr = "";
		
		if(chromosomeName == Definitions.chromosomeM){
			chrmStr = "chrM";
		}else if(chromosomeName == Definitions.chromosomeX){
			chrmStr = "chrX";
		}else if(chromosomeName == Definitions.chromosomeY){
			chrmStr = "chrY";
		}else{
			chrmStr = "chr" + (chromosomeName);
		}
		
		if(strand == Definitions.genomicStrandPOSITIVE){
			strandStr = "+";
		}else {
			strandStr = "-";
		}
		
		
		return chrmStr
				+ "\t" + startLocation + "\t"
				+ stopLocation + "\t" + sequence + "\t" + score
				+ "\t" + strandStr + "\t" + restOfLine;
	}
	/**
	 * @return the chromosomeId
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
	 * @return
	 */
	public String getRestOfLine(){
		return restOfLine;
	}
	/**
	 * @return the startLocType
	 */
	public int getStartLocType() {
		return startLocType;
	}

	/**
	 * @return the stopLocType
	 */
	public int getStopLocType() {
		return stopLocType;
	}

	/**
	 * @return the isSameSubUnit
	 */
	public boolean isSameSubUnit() {
		return isSameSubUnit;
	}

	/**
	 * The location type as defined by the Definitions file.
	 * @param startLocType the startLocType to set
	 * @return true if this type was changed, false otherwise
	 */
	public void setStartLocType(int startLocType) {
		this.startLocType = startLocType;
	}
	

	/**
	 *	The location type as defined by the Definitions file.
	 * @param stopLocType the stopLocType to set
	 * @return true if this type was changed, false otherwise
	 */
	public void setStopLocType(int stopLocType) {
		this.stopLocType = stopLocType;
	}

	/**
	 * Determines whether the start and stop location fall within the same subunit.
	 * @param isSameSubUnit the isSameSubUnit to set
	 */
	public void setSameSubUnit(boolean isSameSubUnit) {

		this.isSameSubUnit = isSameSubUnit;
	}
	

	/**
	 * Compares this LocationLine to another.  This method sorts on chromosomeName, then on start location.
	 */
	public int compareTo(PeptideLocationLine p) {
		if(this.chromosomeName < p.getChromosomeName()){
			return -1;
		}
		if(this.chromosomeName > p.getChromosomeName()){
			return 1;
		}
		if(this.startLocation < p.getStartLocation()){
			return -1;
		}
		if(this.startLocation > p.getStartLocation()){
			return 1;
		}

		return 0;
	}
	
	
	/**
	 * convertLocTypeToString allows for a locType from the definitions file to be passed to it, and it returns the string representation of that type.
	 * ex.  Definitions.LOCATION_TYPE_EXON returns "EXON"
	 * @param locType  The locType to convert.
	 * @return The string conversion of the passed in locType.
	 */
	public static String convertLocTypeToString(int locType){
		String out = null;
		switch(locType){
			case Definitions.LOCATION_TYPE_EXON:
				out = "EXON";
				break;
			case Definitions.LOCATION_TYPE_TRANSCRIPT:
				out = "TRANSCRIPT";
				break;
			case Definitions.LOCATION_TYPE_INTRON:
				out = "INTRON";
				break;
			case Definitions.LOCATION_TYPE_INTERGENIC:
				out = "INTERGENIC";
				break;
			case Definitions.LOCATION_TYPE_UNKNOWN:
				out = "UNKNOWN";
				break;
			default:
				out = "ERROR OF LOCATION TYPE";
				
		}
		
		
		return out;
	}
	
	
}
