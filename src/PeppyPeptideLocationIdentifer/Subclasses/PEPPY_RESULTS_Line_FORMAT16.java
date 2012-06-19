package PeppyPeptideLocationIdentifer.Subclasses;

import java.util.ArrayList;


/**
 * This object represents a single line a peppy results file.
 * @author David tab Thomas
 *
 */
public class PEPPY_RESULTS_Line_FORMAT16 {

	private int ID;
	private int SpectrumID;
	private String MD5;
	private String FileName;
	private double score;
	private double PrecursorMZ;
	private double PrecursorNeutralMass;
	private double eValue;
	private String peptideSequence;
	private int start;
	private int stop;
	private String proteinName;
	private int matchRank;
	private int RankCount;
	private int IonCount;
	private boolean labeled;
	private int charge; 
	private int cleavageAcidCount;
	private boolean inORF;
	private double hydrophobic;
	private double hydrophilic;
	boolean isModified;
	double modMass;
	int modIndex;
	
	/*Appended Information*/
	private int chrm;
	private int strand;
	ArrayList<Integer> peptideStarts = new ArrayList<Integer>();
	ArrayList<Integer> peptideStops = new ArrayList<Integer>();
	
	
	/**
	 * @param iD
	 * @param spectrumID
	 * @param mD5
	 * @param fileName
	 * @param score
	 * @param precursorMZ
	 * @param precursorNeutralMass
	 * @param eValue
	 * @param peptideSequence
	 * @param start
	 * @param stop
	 * @param proteinName
	 * @param matchRank
	 * @param rankCount
	 * @param ionCount
	 * @param labeled
	 * @param charge
	 * @param cleavageAcidCount
	 * @param inORF
	 * @param hydrophobic
	 * @param hydrophilic
	 */
	public PEPPY_RESULTS_Line_FORMAT16(int iD, int spectrumID, String mD5,
			String fileName, double score, double precursorMZ,
			double precursorNeutralMass, double eValue, String peptideSequence,
			int start, int stop, String proteinName, int matchRank,
			int rankCount, int ionCount, boolean labeled, int charge,
			int cleavageAcidCount, boolean inORF, double hydrophobic,
			double hydrophilic, boolean isModified, double modMass, int modIndex) {
		super();
		ID = iD;
		SpectrumID = spectrumID;
		MD5 = mD5;
		FileName = fileName;
		this.score = score;
		PrecursorMZ = precursorMZ;
		PrecursorNeutralMass = precursorNeutralMass;
		this.eValue = eValue;
		this.peptideSequence = peptideSequence;
		this.start = start;
		this.stop = stop;
		this.proteinName = proteinName;
		this.matchRank = matchRank;
		RankCount = rankCount;
		IonCount = ionCount;
		this.labeled = labeled;
		this.charge = charge;
		this.cleavageAcidCount = cleavageAcidCount;
		this.inORF = inORF;
		this.hydrophobic = hydrophobic;
		this.hydrophilic = hydrophilic;
		this.isModified = isModified;
		this.modMass = modMass;
		this.modIndex = modIndex;
	}


	
	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String tab = "\t";
		
		return  + ID + tab
				+ SpectrumID + tab
				+ MD5 + tab
				+ FileName+ tab
				+ score + tab
				+ PrecursorMZ+ tab
				+ PrecursorNeutralMass + tab 
				+ eValue + tab 
				+ peptideSequence + tab
				+ start + tab 
				+ stop + tab 
				+ proteinName + tab
				+ matchRank + tab 
				+ RankCount + tab 
				+ IonCount + tab
				+ labeled + tab
				+ charge + tab
				+ cleavageAcidCount + tab
				+ inORF + tab
				+ hydrophobic + tab
				+ hydrophilic + tab
				+ isModified + tab
				+ modMass + tab
				+ modIndex + tab;
	}





	/**
	 * @return the iD
	 */
	public int getID() {
		return ID;
	}


	/**
	 * @return the spectrumID
	 */
	public int getSpectrumID() {
		return SpectrumID;
	}


	/**
	 * @return the mD5
	 */
	public String getMD5() {
		return MD5;
	}


	/**
	 * @return the fileName
	 */
	public String getFileName() {
		return FileName;
	}


	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}


	/**
	 * @return the precursorMZ
	 */
	public double getPrecursorMZ() {
		return PrecursorMZ;
	}


	/**
	 * @return the precursorNeutralMass
	 */
	public double getPrecursorNeutralMass() {
		return PrecursorNeutralMass;
	}


	/**
	 * @return the eValue
	 */
	public double geteValue() {
		return eValue;
	}


	/**
	 * @return the peptideSequence
	 */
	public String getPeptideSequence() {
		return peptideSequence;
	}


	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}


	/**
	 * @return the stop
	 */
	public int getStop() {
		return stop;
	}


	/**
	 * @return the proteinName
	 */
	public String getProteinName() {
		return proteinName;
	}


	/**
	 * @return the matchRank
	 */
	public int getMatchRank() {
		return matchRank;
	}


	/**
	 * @return the rankCount
	 */
	public int getRankCount() {
		return RankCount;
	}


	/**
	 * @return the ionCount
	 */
	public int getIonCount() {
		return IonCount;
	}


	/**
	 * @return the labeled
	 */
	public boolean isLabeled() {
		return labeled;
	}


	/**
	 * @return the charge
	 */
	public double getCharge() {
		return charge;
	}


	/**
	 * @return the cleavageAcidCount
	 */
	public int getCleavageAcidCount() {
		return cleavageAcidCount;
	}


	/**
	 * @return the inORF
	 */
	public boolean isInORF() {
		return inORF;
	}


	/**
	 * @return the hydrophobic
	 */
	public double getHydrophobic() {
		return hydrophobic;
	}


	/**
	 * @return the hydrophilic
	 */
	public double getHydrophilic() {
		return hydrophilic;
	}



	/*Appended variables*/

	/**
	 * @return the chrm
	 */
	public int getChrm() {
		return chrm;
	}
	/**
	 * @return the strand
	 */
	public int getStrand() {
		return strand;
	}
	/**
	 * @return the peptideStarts
	 */
	public ArrayList<Integer> getPeptideStarts() {
		return peptideStarts;
	}
	/**
	 * @return the peptideStops
	 */
	public ArrayList<Integer> getPeptideStops() {
		return peptideStops;
	}
	/**
	 * @param chrm the chrm to set
	 */
	public void setChrm(int chrm) {
		this.chrm = chrm;
	}
	/**
	 * @param strand the strand to set
	 */
	public void setStrand(int strand) {
		this.strand = strand;
	}
	/**
	 * @param peptideStarts the peptideStarts to set
	 */
	public void setPeptideStarts(ArrayList<Integer> peptideStarts) {
		this.peptideStarts = peptideStarts;
	}
	/**
	 * @param peptideStops the peptideStops to set
	 */
	public void setPeptideStops(ArrayList<Integer> peptideStops) {
		this.peptideStops = peptideStops;
	}

	

}
