package PeppyPeptideLocationIdentifer.Subclasses;

import java.util.ArrayList;




/**
 * This object represents a single line a peppy results file.
 * @author David tab Thomas
 *
 */
public class PEPPY_RESULTS_Line_FORMAT16 {

	//Columns for variois 
	private final static String START = "start";
	private final static String STOP = "stop";
	private final static String SEQUENCE_NAME = "sequencename";
	private final static String STRAND = "strand";
	private final static String PEPTIDE_SEQUENCE = "peptidesequence";
	
	
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
	private boolean isModified;
	private double modMass;
	private int modIndex;
	
	private ArrayList<String> lineValue;
	private ArrayList<String> headerValues;
	private String header;
	
	
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
		//Call default super
		super();
		
		//loadData
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
		
		//Null out the forced format
		lineValue = null;
		headerValues = null;
		header = null;
	}//PEPPY_RESULTS_Line_FORMAT16
	
	
	public PEPPY_RESULTS_Line_FORMAT16(ArrayList<String> lineValueHash, String header){
		//Call default super
		super();
		
		String[] chunks = header.split("\\s");
		
		headerValues = new ArrayList<String>();
		
		for(int i = 0; i < chunks.length; i++){
			this.headerValues.add(chunks[i]);
		}
		
		//loadData
		this.lineValue = lineValueHash;
		this.header = header;
		
		//null out the forced format
		ID = -1;
		SpectrumID = -1;
		MD5 = null;
		FileName = null;
		this.score = -1;
		PrecursorMZ = -1;
		PrecursorNeutralMass = -1;
		this.eValue = -1;
		this.peptideSequence = null;
		this.start = -1;
		this.stop = -1;
		this.proteinName = null;
		this.matchRank = -1;
		RankCount = -1;
		IonCount = -1;
		this.labeled = false;
		this.charge = -1;
		this.cleavageAcidCount = -1;
		this.inORF = false;
		this.hydrophobic = -1;
		this.hydrophilic = -1;
		this.isModified = false;
		this.modMass = -1;
		this.modIndex = -1;
		
	}//PEPPY_RESULTS_Line_FORMAT16
	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String tab = "\t";
		
		//The header is null when the data fields are populated
		if(header == null){
		
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
		
		//Use the hashmap
		}else{
			StringBuffer sb = new StringBuffer();
			
			for(int i = 0; i < lineValue.size(); i++){
				sb.append(lineValue.get(i) + tab);
			}//for
			
			return sb.toString();
		}
	}//toString()





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
		if(header == null){
			return peptideSequence;
		}else{
			for(int i = 0; i < headerValues.size(); i++){
				if(headerValues.get(i).equalsIgnoreCase(PEPTIDE_SEQUENCE)){
					return lineValue.get(i);
				}
				
			}//for
			return null;
		}
	}


	/**
	 * @return the start -1 if none found
	 */
	public int getStart() {
		
		if(header == null){
			return start;
		}else{
			for(int i = 0; i < headerValues.size(); i++){
				if(headerValues.get(i).equalsIgnoreCase(START)){
					return Integer.parseInt(lineValue.get(i));
				}
				
			}//for
			return -1;
		}//else
	}//getStart


	/**
	 * @return the stop -1 if none found
	 */
	public int getStop() {
		if(header == null){
			return stop;
		}else{
			for(int i = 0; i < headerValues.size(); i++){
				if(headerValues.get(i).equalsIgnoreCase(STOP)){
					return Integer.parseInt(lineValue.get(i));
				}
				
			}//for
			return -1;
		}//else
	}


	/**
	 * @return the proteinName  null if none found
	 */
	public String getProteinName() {
		if(header == null){
			return proteinName;
		}else{
			for(int i = 0; i < headerValues.size(); i++){
				if(headerValues.get(i).equalsIgnoreCase(SEQUENCE_NAME)){
					return lineValue.get(i);
				}
				
			}//for
			return null;
		}//else
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
	 * @return the strand  -1 if there is no strand value
	 */
	public int getStrand() {
		if(header == null){
			return strand;
		
		}else{
			for(int i = 0; i < headerValues.size(); i++){
				if(headerValues.get(i).equalsIgnoreCase(STRAND)){
					return Integer.parseInt(lineValue.get(i));
				}
				
			}//for
			return -1;
		}//else
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
