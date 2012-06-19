package GENCODESubmissionTool.Subclasses;

import PersonalProteome.Definitions;

public class PeppyLine implements Comparable<PeppyLine> {
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
	private String sequenceFile;
	private String sequenceDescription;
	private int intronStart;
	private int intronStop;
	private int strand;
	private boolean isSpliced;
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
	private int ChromosomeName = -1;
//	private int Strand;
	private int StartLocation;
	private int StopLocation;
	private int BlockCount; 
	private String BlockSize;
	private String BlockStart;

	private int locationCount = 0;
	private int index = -1;
	private String uniqueKey = "";
	private boolean dna = true;
	private boolean both = false;
	private int bedScore = -1;
	
//	public boolean delete = false;
	
	/**
	 * 
	 * UNMODIFIED
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
	 * @param sequenceFile
	 * @param sequenceDescription
	 * @param intronStart
	 * @param intronStop
	 * @param matchRank
	 * @param rankCount
	 * @param ionCount
	 * @param labeled
	 * @param charge
	 * @param cleavageAcidCount
	 * @param inORF
	 * @param hydrophobic
	 * @param hydrophilic
	 * @param isModified
	 * @param modMass
	 * @param modIndex
	 */
	public PeppyLine(int iD, int spectrumID, String mD5, String fileName,
			double score, double precursorMZ, double precursorNeutralMass,
			double eValue, String peptideSequence, int start, int stop,
			String sequenceFile, String sequenceDescription, int intronStart,
			int intronStop, int strand, boolean isSpliced, int matchRank, int rankCount, int ionCount,
			boolean labeled, int charge, int cleavageAcidCount, boolean inORF,
			double hydrophobic, double hydrophilic, boolean isModified,
			double modMass, int modIndex) {
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
		this.sequenceFile = sequenceFile;
		this.sequenceDescription = sequenceDescription;
		this.intronStart = intronStart;
		this.intronStop = intronStop;
		this.strand = strand;
		this.isSpliced = isSpliced;
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

	/**
	 * 
	 *MODIFIED
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
	 * @param sequenceFile
	 * @param sequenceDescription
	 * @param intronStart
	 * @param intronStop
	 * @param matchRank
	 * @param rankCount
	 * @param ionCount
	 * @param labeled
	 * @param charge
	 * @param cleavageAcidCount
	 * @param inORF
	 * @param hydrophobic
	 * @param hydrophilic
	 * @param isModified
	 * @param modMass
	 * @param modIndex
	 * @param chromosomeName
	 * @param strand
	 * @param startLocation
	 * @param stopLocation
	 * @param blockCount
	 * @param blockSize
	 * @param blockStart
	 */
	public PeppyLine(int iD, int spectrumID, String mD5, String fileName,
			double score, double precursorMZ, double precursorNeutralMass,
			double eValue, String peptideSequence, int start, int stop,
			String proteinName, int matchRank, int rankCount, int ionCount,
			boolean labeled, int charge, int cleavageAcidCount, boolean inORF,
			double hydrophobic, double hydrophilic, boolean isModified,
			double modMass, int modIndex, int chromosomeName, int Strand,
			int startLocation, int stopLocation, int blockCount,
			String blockSize, String blockStart) {
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
		this.ChromosomeName = chromosomeName;
		this.strand = Strand;
		this.StartLocation = startLocation;
		this.StopLocation = stopLocation;
		this.BlockCount = blockCount;
		this.BlockSize = blockSize;
		this.BlockStart = blockStart;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String tab = "\t";
		String restOfLine = "";
//		String pepStrand = "+";
//		if(strand == Definitions.genomicStrandNEGATIVE){
//			pepStrand = "-";
//		}
		String strand = "+";
		if(this.strand == Definitions.genomicStrandNEGATIVE){
			strand = "-";
		}
		String midChunk = sequenceFile + tab
		+ sequenceDescription + tab 
		+ intronStart + tab 
		+ intronStop + tab 
		+ strand + tab 
		+ isSpliced;
		if(ChromosomeName != -1){
			
			midChunk = proteinName;
			restOfLine = Definitions.convertChrmNumToString(ChromosomeName) + tab + strand + tab + StartLocation + tab + StopLocation + tab + BlockCount + tab + BlockSize + tab + BlockStart;
		}

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
				+ midChunk + tab
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
				+ modIndex + tab
				+ restOfLine;
	}//toString

	
	
	/**
	 * @return the bedScore
	 */
	public int getBedScore() {
		return bedScore;
	}

	/**
	 * @param bedScore the bedScore to set
	 */
	public void setBedScore(int bedScore) {
		this.bedScore = bedScore;
	}

	/**
	 * @return the dna
	 */
	public boolean isDna() {
		return dna;
	}

	/**
	 * @return the both
	 */
	public boolean isBoth() {
		return both;
	}

	/**
	 * @param both the both to set
	 */
	public void setBoth(boolean both) {
		this.both = both;
	}

	/**
	 * @param dna the dna to set
	 */
	public void setDna(boolean dna) {
		this.dna = dna;
	}

	/**
	 * @return the uniqueKey
	 */
	public String getUniqueKey() {
		return uniqueKey;
	}

	/**
	 * @param uniqueKey the uniqueKey to set
	 */
	public void setUniqueKey(String uniqueKey) {
		this.uniqueKey = uniqueKey;
	}

	public void setIndex(int index){
		this.index = index;
	}
	
	public int getIndex(){
		return this.index;
	}
	
	public void setLocCount(int LocCount){
		this.locationCount = LocCount;
	}
	
	public int getlocationCount(){
		return this.locationCount;
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
	 * The genome based location in genome searches, the protein based location in protein searches.
	 * @return the start
	 */
	public int getStart() {
		return start;
	}

	/**
	 * The genome based location in genome searches, the protein based location in protein searches.
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
	 * @return the sequenceFile
	 */
	public String getSequenceFile() {
		return sequenceFile;
	}

	/**
	 * @return the sequenceDescription
	 */
	public String getSequenceDescription() {
		return sequenceDescription;
	}

	/**
	 * @return the intronStart
	 */
	public int getIntronStart() {
		return intronStart;
	}

	/**
	 * @return the intronStop
	 */
	public int getIntronStop() {
		return intronStop;
	}

	/**
	 * @return the strand
	 */
	public int getFirstStrand() {
		return strand;
	}

	/**
	 * @return the isSpliced
	 */
	public boolean isSpliced() {
		return isSpliced;
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
	public int getCharge() {
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

	/**
	 * @return the isModified
	 */
	public boolean isModified() {
		return isModified;
	}

	/**
	 * @return the modMass
	 */
	public double getModMass() {
		return modMass;
	}

	/**
	 * @return the modIndex
	 */
	public int getModIndex() {
		return modIndex;
	}

	/**
	 * @return the chromosomeName
	 */
	public int getChromosomeName() {
		return ChromosomeName;
	}
//
//	/**
//	 * @return the strand
//	 */
//	public int getSecondStrand() {
//		return Strand;
//	}

	/**
	 * The Genome based location in proteins.
	 * @return the startLocation
	 */
	public int getStartLocation() {
		return StartLocation;
	}

	/**
	 * The Genome based location in proteins.
	 * @return the stopLocation
	 */
	public int getStopLocation() {
		return StopLocation;
	}

	/**
	 * @return the blockCount
	 */
	public int getBlockCount() {
		return BlockCount;
	}

	/**
	 * @return the blockSize
	 */
	public String getBlockSize() {
		return BlockSize;
	}

	/**
	 * @return the blockStart
	 */
	public String getBlockStart() {
		return BlockStart;
	}

	/**
	 * Sorts based on MD5, Score, Peptide Sequence, Gen Location.  returns -1 if this object is less then the specified object, 0 if equal, and 1 if greater than.
	 */
	public int compareTo(PeppyLine pl) {
		//MD5, Score Descending, Peptide Sequence, Gen Location
		if(this.MD5.compareTo(pl.getMD5()) < 0){
			return -1;
		}else if(this.MD5.compareTo(pl.getMD5()) > 0){
			return 1;
		}
		
		if(this.getScore() > pl.getScore()){
			return -1;
		}else if(this.getScore() < pl.getScore()){
			return 1;
		}
		
		if(this.getPeptideSequence().compareTo(pl.getPeptideSequence()) < 0){
			return -1;
		}else if(this.getPeptideSequence().compareTo(pl.getPeptideSequence()) > 0){
			return 1;
		}
	
		int thisStart = 0;
		int plStart = 0;
		
		if(this.isDna()){
			thisStart = this.getStart();
		}else{
			thisStart = this.getStartLocation();
		}
		if(pl.isDna()){
			plStart = pl.getStart();
		}else{
			plStart = pl.getStartLocation();
		}
		
		if(thisStart < plStart){
			return -1;
		}else if(thisStart > plStart){
			return 1;
		}
		
		
		
		return 0;
	}


}//PeppyLine
