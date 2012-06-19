package PersonalProteome;

/**
 * GENCODE_GTF_Line is a object representation of a single line in a GENCODE GTF format file.  It contains all of the information from the annotation file it came from.
 * Information about the GENCODE gtf format can be found: http://www.gencodegenes.org/gencodeformat.html
 * @author David "Corvette" Thomas
 *
 */

public class GENCODE_GTF_Line{

	/*Variables*/
	//mandatory
	//iD represents a GENCODES GTF Lines place in the GTFlist
	private final int iD;
	private int chromosomeName;
	private int annotationSource;
	private int featureType;
	private int startLocation;
	private int stopLocation;
	private int score;
	private int genomicStrand;
	private int genomicPhase;
	private String geneID;
	private String transcriptID;
	private String gene_Type;
	private int gene_Status;
	private String gene_Name;
	private String transcript_Type;
	private int transcript_Status;
	private String transcript_Name;
	private int level;
	
	//optional
	private String tag;
	private String ccdsid;
	private String havana_Gene;
	private String havana_Transcript;
	private String ont;
	
	
	/**
	 * Use constants in Peppy.Definitions to populate each variable.
	 * Constructor for mandatory fields only.
	 */
	public GENCODE_GTF_Line(int iD, int chromosomeName, int annotationSource, int featureType, int startLocation, int stopLocation, int score, int genomicStrand, int genomicPhase,
								String geneID, String transcriptID, String gene_Type, int gene_Status, String gene_Name, String transcript_Type, int transcript_Status,
								String transcript_Name, int level){
		//Set all of the mandatory fields
		this.iD = iD;
		this.chromosomeName = chromosomeName;
		this.annotationSource = annotationSource;
		this.featureType = featureType;
		this.startLocation = startLocation;
		this.stopLocation = stopLocation;
		this.score = score;
		this.genomicStrand = genomicStrand;
		this.genomicPhase = genomicPhase;
		this.geneID = geneID;
		this.transcriptID = transcriptID;
		this.gene_Type = gene_Type;
		this.gene_Status = gene_Status;
		this.gene_Name = gene_Name;
		this.gene_Status = gene_Status;
		this.gene_Name = gene_Name;
		this.transcript_Type = transcript_Type;
		this.transcript_Status = transcript_Status;
		this.transcript_Name = transcript_Name;
		this.level = level;
	
		
	}
	
	/**
	 * Use constants in 'Definitions' to populate each variable.
	 * Constructor for mandatory and optional fields.
	 */
	public GENCODE_GTF_Line(int iD, int chromosomeName, int annotationSource, int featureType, int startLocation, int stopLocation, int score, int genomicStrand, int genomicPhase,
								String geneID, String transcriptID, String gene_Type, int gene_Status, String gene_Name, String transcript_Type, int transcript_Status,
								String transcript_Name, int level, String tag, String ccdsid, String havana_Gene, String havana_Transcript, String ont){
		//Set all of the mandatory fields
		this(iD, chromosomeName,  annotationSource,  featureType,  startLocation,  stopLocation,  score,  genomicStrand,  genomicPhase,
				 geneID,  transcriptID,  gene_Type,  gene_Status,  gene_Name,  transcript_Type,  transcript_Status,
				 transcript_Name,  level);
		//Set the optional fields
		this.tag = tag;
		this.ccdsid = ccdsid;
		this.havana_Gene = havana_Gene;
		this.havana_Transcript = havana_Transcript;
		this.ont = ont;
	}

	public String toString(){
	
		return  iD + " " + chromosomeName + " " +   annotationSource + " " +   featureType + " " +   startLocation + " " +   stopLocation + " " +   score + " " +   genomicStrand + " " +   genomicPhase + " " + 
				 geneID + " " +   transcriptID + " " +   gene_Type + " " +   gene_Status + " " +   gene_Name + " " +   transcript_Type + " " +  transcript_Status + " " +  transcript_Name + " " + 
				 level + " " +   tag + " " +   ccdsid + " " +   havana_Gene + " " +   havana_Transcript + " " +   ont;
	}

	/**
	 * @return the chromosomeName
	 */
	public int getiD() {
		return iD;
	}
	
	/**
	 * @return the chromosomeName
	 */
	public int getChromosomeName() {
		return chromosomeName;
	}


	/**
	 * @return the annotationSource
	 */
	public int getAnnotationSource() {
		return annotationSource;
	}


	/**
	 * @return the featureType
	 */
	public int getFeatureType() {
		return featureType;
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
	 * @return the score
	 */
	public int getScore() {
		return score;
	}


	/**
	 * @return the genomicStrand
	 */
	public int getGenomicStrand() {
		return genomicStrand;
	}


	/**
	 * @return the genomicPhase
	 */
	public int getGenomicPhase() {
		return genomicPhase;
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
	 * @return the gene_Type
	 */
	public String getGene_Type() {
		return gene_Type;
	}


	/**
	 * @return the gene_Status
	 */
	public int getGene_Status() {
		return gene_Status;
	}


	/**
	 * @return the gene_Name
	 */
	public String getGene_Name() {
		return gene_Name;
	}


	/**
	 * @return the transcript_Type
	 */
	public String getTranscript_Type() {
		return transcript_Type;
	}


	/**
	 * @return the transcirpt_Status
	 */
	public int getTranscirpt_Status() {
		return transcript_Status;
	}


	/**
	 * @return the transcript_Name
	 */
	public String getTranscript_Name() {
		return transcript_Name;
	}


	/**
	 * @return the level
	 */
	public int getLevel() {
		return level;
	}


	/**
	 * @return the tag
	 */
	public String getTag() {
		return tag;
	}


	/**
	 * @return the ccdsid
	 */
	public String getCcdsid() {
		return ccdsid;
	}


	/**
	 * @return the havana_gene
	 */
	public String getHavana_gene() {
		return havana_Gene;
	}


	/**
	 * @return the havana_transcript
	 */
	public String getHavana_transcript() {
		return havana_Transcript;
	}


	/**
	 * @return the ont
	 */
	public String getOnt() {
		return ont;
	}
	

	
}
