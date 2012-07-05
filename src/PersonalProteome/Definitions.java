package PersonalProteome;

/**
 * A class to store constants.  All variables should be declared final and static.
 * @author Brian Risk,  David "Corvette" Thomas
 *
 */
public class Definitions {	

	public final static char[] aminoAcidList =     {'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','.','Y','.','Y','S','S','S','S','.','C','W','C','L','F','L','F'};
	public final static char[] mitoAminoAcidList = {'K','N','K','N','T','T','T','T','.','S','.','S','M','I','M','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','.','Y','.','Y','S','S','S','S','W','C','W','C','L','F','L','F'};
	
	
	public final static char NULL_CHAR = '%';
	
	/*Generic date format to use accross prersonal proteome*/
	public final static String DATE_FORMAT = "MM-dd-yyyy HH:mm:ss";
	/*GENCODE GTF values*/
	//chromosomeNames
	public final static int chromosomeM = 23;
	public final static int chromosomeX = 24;
	public final static int chromosomeY = 25;
	//Annotation Source
	public final static int ENSEMBL = 0;
	public final static int HAVANA = 1;
	//feature-type
	public final static int featureTypeGENE = 0;
	public final static int featureTypeTRANSCRIPT = 1;
	public final static int featureTypeEXON = 2;
	public final static int featureTypeCDS = 3;
	public final static int featureTypeSTART_CODON = 4;
	public final static int featureTypeSTOP_CODON = 5;
	public final static int featureTypeUTR = 6;
	public final static int featureTypeSELENOCYSTEINE = 7;
	//Strand
	public final static int genomicStrandPOSITIVE = 0;
	public final static int genomicStrandNEGATIVE = 1;
	//Phase
	public final static int genomicPhaseZERO = 0;
	public final static int genomicPhaseONE = 1;
	public final static int genomicPhaseTWO = 2;
	public final static int genomicPhasePERIOD = 3;
	public final static char genomicPhaseCHAR = 'X';
	//Gene Status
	public final static int geneStatusKNOWN = 0;
	public final static int geneStatusNULL = 1;
	public final static int geneStatusNOVEL = 2;
	public final static int geneStatusUNKNOWN = 3;
	//Transcript Status
	public final static int transcriptStatusKNOWN = 0;
	public final static int transcriptStatusNULL = 1;
	public final static int transcriptStatusNOVEL = 2;
	public final static int transcriptStatusUNKNOWN = 3;
	//GENCODE Level
	public final static int GENCODELevelONE = 0;
	public final static int GENCODELevelTWO = 1;
	public final static int GENCODELevelTHREE = 2;
	
	//SpecialCharacters used to delimit sequences
	public final static char END_NOT_FOUND_CHAR = '+';
	public final static String END_NOT_FOUND_SEQUENCE = "+++";
	public final static char START_NOT_FOUND_CHAR = '&';
	public final static String START_NOT_FOUND_SEQUENCE = "&&&";
	public final static char SELENO_CHAR = '#';
	public final static String SELENO_SEQUENCE = "###";
	public final static char SELENO_AMINO_ACID = 'U';
	
	//Statistical error codes for Personal Proteome
	public final static int ERROR_DELAYED_START = 0;
	public final static int ERROR_EXTENDED_STOP = 1;
	public final static int ERROR_SELENOCYSTEINE = 2;
	
	//Location Type
	public final static int LOCATION_TYPE_UNKNOWN = -1;
	public final static int LOCATION_TYPE_INTERGENIC = 0;
	public final static int LOCATION_TYPE_TRANSCRIPT = 1;
	public final static int LOCATION_TYPE_INTRON = 2;
	public final static int LOCATION_TYPE_EXON = 3;
	
	//Match Type
	public final static int MATCH_TYPE_MATCHED = 1;
	public final static int MATCH_TYPE_UNMATCHED = 0;
	
	//Protein Types
	public final static String PROTEIN_TYPE_TOO_SHORT = "PROTEIN TOO SHORT";
	public final static String STOP_CODON_STRING = ".";
	

	/**
	 * Assumes you are using the Definitions File chromosomes
	 * @param ChrmNumber
	 * @return the string representation of this chromosome
	 */
	public static String convertChrmNumToString(int chrmNum){
		String out = "chr";
		
		//Consider X, Y, M
		switch (chrmNum) {
			case Definitions.chromosomeY:
				out += "Y";
				break;
			case Definitions.chromosomeX:
				out += "X";
				break;
			case Definitions.chromosomeM:
				out += "M";
				break;
			default:
				out += chrmNum;
				break;
		}
		return out;
	}//convertChrmNumToString
	
	/**
	 * Assumes you are using the Definitions File chromosomes
	 * @param ChrmNumber
	 * @return the string representation of this chromosome
	 */
	public static String convertStrandToString(int strand){
		String out = "!";
		
		if(strand == Definitions.genomicStrandNEGATIVE){
			out = "-";
		}else if(strand == Definitions.genomicStrandPOSITIVE){
			out = "+";
		}

		return out;
	}//convertStrandToString

	
}
