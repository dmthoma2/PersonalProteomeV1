package MiscTools.GENCODEPPCompariosonTool;

import PersonalProteome.Definitions;

/**
 * Translation stores a single protein with it associated ID.
 * @author David "Corvette" Thomas
 *
 */
public class Translation {

	private String geneID;
	private String transcriptID;
	private String sequence;
	private int matchType = Definitions.MATCH_TYPE_UNMATCHED;
	private int variants = 0;

	/**
	 * @param geneID
	 * @param transcriptID
	 * @param sequence
	 */
	public Translation(String geneID, String transcriptID, String sequence) {
		super();
		this.geneID = geneID;
		this.transcriptID = transcriptID;
		this.sequence = sequence;
	}
	
	/**
	 * @return a string representation of this object.  TranscirptID|GeneID|Sequence
	 */
	public String toString(){
		return transcriptID + "|" + geneID + "|" + sequence;
		
		
	}
	/**
	 * @return the geneID
	 */
	public String getGeneID() {
		return geneID;
	}
	/**
	 * @return the transcriptID
	 */
	public String getTranscriptID() {
		return transcriptID;
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}
	/**
	 * 
	 */
	public void setSequence(String newSeq){
		this.sequence = newSeq;
	}
	/**
	 * Based on matched type in definitions.
	 * @return the matchType
	 */
	public int getMatchType() {
		return matchType;
	}
	/**
	 * Set based on matched type in definitions.
	 * @param matchType the matchType to set
	 */
	public void setMatchType(int matchType) {
		this.matchType = matchType;
	}
	
	/**
	 * @return
	 */
	public int getVariants(){
		return variants;
	}
	/**
	 * Increments the variant count by 1.
	 */
	public void incrementVariants(){
		variants++;
	}
	
	
}
