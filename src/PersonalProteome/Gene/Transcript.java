package PersonalProteome.Gene;


import java.util.ArrayList;
import java.util.Collections;


import PersonalProteome.Definitions;
import PersonalProteome.U;


/**
 * Transcript represents a single transcript from an annotation file.  It contains all of the information to create a single protein.
 *   It extends GeneticObject.  <GeneticObject>
 * @author David "Corvette" Thomas
 *
 */
public class Transcript extends GeneticObject {


	//id is the index of this object GENCODE_GTF_Line in the GTFlist
	private int iD;

	//General Variables
	private int transcriptLength;
	private int variantCount = 0;
	private int endExtensionLength = 0;
	private int startExtensionLength = 0;
	private boolean containPreMatureStop = false;
	private int preMatureStopCount = 0;
	private int leftOverAminoLength = 0;
	

	//Error storage variables, default to -1 if no errors are found.
	private int startErrorLoc = -1;
	private int stopErrorLoc = -1;
	

	//Sequences within this transcript
	private String transcriptSequence;
	private String protein;
	
	//Annotation based information about this transcript
	private int strand = -1;
	private ArrayList<CDS> CDS = new ArrayList<CDS>();
	private ArrayList<Double> exonSplitLocations = new ArrayList<Double>();
	private ArrayList<Seleno> Seleno = new ArrayList<Seleno>();
	
	
	/*Amplification Score*/
	private double ampScore;

	
	/*Optional consideration of modified stops and starts*/
	private boolean considerStopsAndStarts;
	
	
	/*mitochondrion DNA related variables*/
	private int chromosome;
	private boolean isMito;
	
	
	private boolean hasStartCodon = true;
	private boolean hasStopCodon = true;
	/**
	 * Transcript represents a single transcript from an annotation file.  It contains all of the information to create a single protein.
	 * @param chromosme  This is the chromosome that this transcript is on.  It is based on the chromosome list in 'Definitions.
	 * @param transcriptSequence This string is the DNA sequence from the chromsome that this transcript covers
	 * @param iD This is the index in the GTFlist of this transcript.  GTFlist.get(transcript.getiD()) can be used to look up all of this transcripts information.
	 * @param start The start location of this transcript.
	 * @param end The end location of this transcript.
	 * @param useModStartStop Boolean variable on whether or not to take into account starts and stops that are checked for accuracy.  If a non-valid stop or start is found, an attempt to find the next avaliable one will be made.
	 */
	public Transcript(int chromosome, String transcriptSequence, int strand, int iD, int start, int end, boolean useModStartStop){
		//Set start and stop in the super
		super(start, end);
		
		//Set chromosome and determine if it is mitochondrion DNA
		this.chromosome = chromosome;
		if(chromosome == Definitions.chromosomeM){
			isMito = true;
		}else{
			isMito = false;
		}
		
		//Sequence information
		this.strand = strand;
		this.transcriptLength = transcriptSequence.length();
		this.transcriptSequence = transcriptSequence;
		
		//Transcript information
		this.iD = iD;
		
		//Modification variables
		this.considerStopsAndStarts = useModStartStop;
		


	}
	
	/**
	 * Returns a string of this transcripts iD, start, and stop.
	 */
	public String toString(){
		return iD + super.toString();
	}

	
	/**
	 * createProtein Takes the information within this transcript and creates a protein.  It handles the modified starts and stops,
	 * it sorts CDS regions, finds CDS split locations, and clears out memory once finished.
	 */
	public void createProtein(){
		
		//If the memory has been cleared or there was no sequence ever added, then there is nothing to do since a protein cannot be re/created.
		if(transcriptSequence == null){
			return;
		}
		//sb holds the dnaSequence for the protein as it is being built
		StringBuffer sb = new StringBuffer();
		
		//If there are no CDS regions, then report that there are none for this protein
		if(CDS.size() == 0){
			protein = "No CDS";
			saveMemory();
			return;
		}
		
		//Ensure the CDS are in the correct order
		Collections.sort(CDS);


		//Get the dnaSequence that the body would produce from the transcript sequence based on the CDS
		for(int i = 0; i < CDS.size(); i++){
			if(CDS.size() == 1 && CDS.get(0).getLength() < 3){
				sb.append(transcriptSequence.substring(CDS.get(i).getStart() - getStart(), CDS.get(i).getStart() - getStart() + 2 + 1));
			}else{
				sb.append(transcriptSequence.substring(CDS.get(i).getStart() - getStart(), CDS.get(i).getStop() - getStart() + 1));
			}
			
		}//CDS Loop
		

		//Put the completed sequence into a string to work on.
		String dnaSequence = sb.toString();
		StringBuffer buildingProtein = new StringBuffer();

		//Without a start codon, identify the genomic phase of this protein.
		if(!hasStartCodon){
			if(strand == Definitions.genomicStrandPOSITIVE){//Positive STrand
				if(CDS.get(0).getGenomicPhase() == Definitions.genomicPhaseONE || CDS.get(0).getGenomicPhase() == Definitions.genomicPhaseTWO){
					buildingProtein.append(Definitions.genomicPhaseCHAR);
					dnaSequence = dnaSequence.substring(CDS.get(0).getGenomicPhase());
				}//Check for phase 1 or 2
			}else{//Negative Strand
				int lastCDS = CDS.size() - 1;
				if(CDS.get(lastCDS).getGenomicPhase() == Definitions.genomicPhaseONE || CDS.get(lastCDS).getGenomicPhase() == Definitions.genomicPhaseTWO){
					buildingProtein.append(Definitions.genomicPhaseCHAR);
					dnaSequence = dnaSequence.substring(0, dnaSequence.length() - CDS.get(lastCDS).getGenomicPhase());
				}//Check for phase 1 or 2
				
			}//if/else
		}//Does not have start codon

		leftOverAminoLength = dnaSequence.length() % 3;
		
		//If a protein is to short to form any amino acids, return an error
		if(dnaSequence.length() < 3){
			this.protein =  Definitions.PROTEIN_TYPE_TOO_SHORT;
			saveMemory();
			return;
		}//to short if
		
		
				//If modified stops and starts are to be considered, modify the dna
		if(considerStopsAndStarts){
				
			//Consider modified Stops
			//Positive strand
			if(strand == Definitions.genomicStrandPOSITIVE){
				//Confirm that there is some section to work on past the end of the protein
				if(this.getStop() - CDS.get(CDS.size() - 1).getStop() >= 3){
					//Number of extra acids that would not normally be encoded
					int leftOverAminoLength = (dnaSequence.length() % 3);
					//Peel off the left over acids, as they will be compensated for during the extension efforts
					dnaSequence = dnaSequence.substring(0, dnaSequence.length() - leftOverAminoLength);
					
					//Get the 
					String regionToWork = transcriptSequence.substring(CDS.get(CDS.size() - 1).getStop() - this.getStart() + 1 - leftOverAminoLength);
					int codonCount = regionToWork.length() / 3;
					
					//Verify that the first thing after the end of the protein is not a stop
					if(!isAStop(regionToWork.substring(0, 3), strand)){
						stopErrorLoc = CDS.get(CDS.size() - 1).getStop() + 1 - leftOverAminoLength;
						dnaSequence += Definitions.END_NOT_FOUND_SEQUENCE;

						//add to the sequence until another stop is found.
						for(int l = 0; l < codonCount; l++){
							String codon = regionToWork.substring(l*3, (l+1)*3);
							
							if(!isAStop(codon, strand)){		
								endExtensionLength++;
									dnaSequence += codon;
							}else{
								//A stop is found so break;
								break;
							}//else
						}//for loop
					}//If to check to see if it is a stop initially

				}//inital setup
				
				//Write another one of these for the negative strand, peeling away from the front instead of the back;	
			}else if(strand == Definitions.genomicStrandNEGATIVE){
				//Confirm that there is some section to work on past the end of the protein
				if(CDS.get(0).getStart() - this.getStart() >= 3){
					int leftOverAminoLength = (dnaSequence.length() % 3);
					dnaSequence = dnaSequence.substring(leftOverAminoLength, dnaSequence.length());

					String regionToWork = transcriptSequence.substring(0, CDS.get(0).getStart() - this.getStart() + leftOverAminoLength);
					

					int codonCount = regionToWork.length() / 3;
					
					//Verify that the first thing after the end of the protein is not a stop
					int regionLength = regionToWork.length();
					if(!isAStop(regionToWork.substring(regionLength - 3, regionLength), strand)){
						stopErrorLoc = CDS.get(0).getStart() + leftOverAminoLength;

						dnaSequence = Definitions.END_NOT_FOUND_SEQUENCE + dnaSequence;
						//add to the sequence until another stop is found.
						for(int l = 0; l < codonCount; l++){
							
							String codon = regionToWork.substring(regionLength - (l + 1)*3, regionLength - (l)*3);

							if(!isAStop(codon, strand)){
								endExtensionLength++;
									dnaSequence = codon + dnaSequence;
	
							}else{
								
								break;
							}//else
						}//for loop
					}//If to check to see if it is a stop initially

				}//Initial setup
			}//Strand	
			
			
			//Consider modified Starts
			//Positive strand
			if(hasStartCodon){
				if(strand == Definitions.genomicStrandPOSITIVE){
					
					//Number of extra acids that would not normally be encoded
					
					int leftOverAminoLength = (dnaSequence.length() % 3);
	
	
	
					//ATG is the forwards stand start 
					int mLocation = -1;
	
					
					//Iterate through string and find a M
					for(int i = 0; i < (dnaSequence.length() - leftOverAminoLength) / 3 ; i++){
						String codon = dnaSequence.substring(i*3 ,(i+1)*3);
	
						//When a start is found, lock in the location and quit
						if(isAStart(codon, strand)){
							mLocation = i*3;
							break;
						}
						
					}//for of codons
					
					//If a error has been found, report it as the startErrorLocation
					if(mLocation != 0){
						startErrorLoc = CDS.get(0).getStart();
					}
					//Check to see if no start was found.
					if(mLocation == -1){
						//For now don't worry about these proteins
						
						
						//Set the mLocation to 1 past the end of the sequence, so padding is added all the way down the length of the string
						mLocation = dnaSequence.length() - leftOverAminoLength;
					}//check if no M was found
					
					//Generate the front part of the sequence
					StringBuffer frontPadding = new StringBuffer();
					for(int v = 0; v < (mLocation) / 3; v++){
						startExtensionLength++;
						frontPadding.append(Definitions.START_NOT_FOUND_SEQUENCE);
					}//padding for
			
					//Create the new, start extended DNA sequence
					dnaSequence = frontPadding.toString() + dnaSequence.substring(mLocation, dnaSequence.length());
					
	
	
					//Check to see that a start was found on this protein somewhere
					
				}else if(strand == Definitions.genomicStrandNEGATIVE){
					//Number of extra acids that would not normally be encoded
					int leftOverAminoLength = (dnaSequence.length() % 3);
					//Peel off the left over acids, as they will be compensated for during the extension efforts
	
	
					dnaSequence = dnaSequence.substring(leftOverAminoLength , dnaSequence.length());
					
	
					
					//Opening soon.
					int mLocation = -1;
	
					for(int l = 0; l < dnaSequence.length() / 3; l++){
						
						String codon = dnaSequence.substring(dnaSequence.length() - (l + 1)*3, dnaSequence.length() - (l)*3);
						if(isAStart(codon, strand)){
							mLocation = dnaSequence.length() - (l+1) * 3;
							break;
						}
					}//for loop
					
					//If a error has been found, report it as the startErrorLocation
					if(mLocation != 0){
						startErrorLoc = CDS.get(CDS.size() - 1).getStop() - 3;
					}
					
					//Generate the front part of the sequence
					StringBuffer backPadding = new StringBuffer();
					for(int v = (mLocation + 3) / 3; v < (dnaSequence.length()) / 3; v++){
						startExtensionLength++;
						backPadding.append(Definitions.START_NOT_FOUND_SEQUENCE);
					}//padding for
					
	
					//Create the new, start extended DNA sequence
					String tempD = dnaSequence;
					dnaSequence = dnaSequence.substring(0, mLocation + 3) + backPadding.toString();
					String tempC = dnaSequence;
					
	
					
					if(tempD.length() != tempC.length()){
	
					}
				}
			}//ifHasStartCodon
		}//consider stops and starts
		


		//Clear the sequence to save on memory

		//Convert the DNA into a Protein
		char [] codon = new char[3];
		char aminoAcid;
		int mod = 0;
		
		int increment = 1;
		int startPosition = 0;

		int stopPosition = dnaSequence.length();
		boolean isForwardsStrand = true;
		int index = 0;
		if (strand == Definitions.genomicStrandNEGATIVE){
			increment = -1;
			startPosition = stopPosition - 1;
			stopPosition = -1;
			isForwardsStrand = false;
		}

		//Iterate through the sequence and convert codons.
		for (index = startPosition ; index != stopPosition; index += increment) {
			codon[mod] = dnaSequence.charAt(index);

			if (mod == 2) {

				if(containPreMatureStop){
					preMatureStopCount++;
				}
				//Signal that a delayed start has been found
				if(codon[0] == Definitions.START_NOT_FOUND_CHAR && codon[1] == Definitions.START_NOT_FOUND_CHAR && codon[2] == Definitions.START_NOT_FOUND_CHAR){
					buildingProtein.append(Definitions.START_NOT_FOUND_CHAR);
					
				//Signal that an extensions has occurred on the end of the protein.
				}else if(codon[0] == Definitions.END_NOT_FOUND_CHAR && codon[1] == Definitions.END_NOT_FOUND_CHAR && codon[2] == Definitions.END_NOT_FOUND_CHAR){
				buildingProtein.append(Definitions.END_NOT_FOUND_CHAR);
					
				
				}else if(codon[0] == Definitions.SELENO_CHAR && codon[1] == Definitions.SELENO_CHAR && codon[2] == Definitions.SELENO_CHAR && Seleno.size() > 0){
				//Seleno found so add a U and continue
					buildingProtein.append(Definitions.SELENO_AMINO_ACID);

				}else{
				//Adjust for mitochondrion dna
				if(isMito){
					aminoAcid = Definitions.mitoAminoAcidList[indexForCodonArray(codon, isForwardsStrand)];

				}else{
					aminoAcid = Definitions.aminoAcidList[indexForCodonArray(codon, isForwardsStrand)];
				}
				if(aminoAcid == '.'){
					containPreMatureStop = true;
					preMatureStopCount++;
				}
				buildingProtein.append(aminoAcid);
				}
				
				/* reset mod */
				mod = 0;
			} else {
				mod++;
			}
		}
		
		this.protein = buildingProtein.toString();

		//Find the split locations
		findCDSSplitLocations();
		saveMemory();
	}//createProtein
	
	
	/**
	 * findCDSSplitLocations determines at what indices of a protein there is a CDS split. Assumes the CDS list is sorted.
	 */
	private void findCDSSplitLocations(){
		//There are no split locations if there is only one/zero CDS
		if(CDS.size() < 2){
			return;
		}

		//Forward strand
		//The CDS are in nucleic acid lengths, but the code needs the lengths in amino acids, so each length must be divided by 3.
		
		//Get the first CDS split location
		exonSplitLocations.add(((double)(CDS.get(0).getLength()) + 1) / 3);
		
		//Add in two-last cds split locations.
		for(int i = 1; i < CDS.size() - 1; i++){
			exonSplitLocations.add((((double)(CDS.get(i).getLength()) + 1) / 3) + exonSplitLocations.get(i - 1));
			
		}

	
		if(strand == Definitions.genomicStrandNEGATIVE){

			//The reverse strand calculates the split locations backwards because the protein is formed backwards.
			double finalLocation = ((double)CDS.get(CDS.size() - 1).getLength() + 1)/3 + exonSplitLocations.get(CDS.size() - 1 - 1);
			
			//Reverse the values based on the length
			//Work backwards to populate the split temp list.
			ArrayList<Double> temp = new ArrayList<Double>();
			for(Double val: exonSplitLocations){
				temp.add(finalLocation - val);
			}

			exonSplitLocations = temp;

			//Reverse the split locations so they are ascending relative to the start of the protein.
			Collections.sort(exonSplitLocations);
		}//Reverse

	
		
		//Error check to ensure each location has the exact value that it needs.  Rounding with doubles can be funny.
		ArrayList<Double> tempAdjust = new ArrayList<Double>();
		for(int i = 0; i < exonSplitLocations.size() ; i++){
			
			 //Deal with truncated numbers becuase comp does not recognize 69.999999999999 = 70
			double val = exonSplitLocations.get(i).doubleValue();
			 int floor = (int)val;
			 int ceiling = floor + 1;
		     double  minValue = floor + .9;
		 
		     if(val < ceiling && val >= minValue ){
		    	 val = ceiling;
		     }else{
		    	 if(val < (int)val + .1){
		    		 val = (int)val;
		    	 }
		     }
			tempAdjust.add(val);
		}
		
		exonSplitLocations = tempAdjust;

	}





	/**
	 * Save memory clears the sequence, which occupies the majority of memory within this transcript.
	 */
	private void saveMemory(){
		transcriptSequence = null;
	}



	/**
	 * Method used to lookup a codons corresponding amino acid based on its nucleic sequence.
	 * @author Brian Risk
	 * @param codon
	 * @param forwards
	 * @return The index into the codon array.
	 */
	public static int indexForCodonArray(char [] codon, boolean forwards) {
		int out = indexForCodonArray(codon);
		if (out == -1) return 50; //if unknown, return STOP
		if (forwards) {
			return indexForCodonArray(codon);
		} else {
			return 63 - indexForCodonArray(codon);
		}
	}
	
	/**
	 * Method used to lookup a codons corresponding amino acid, but assumes the direction is forwards.
	 * @param codon
	 * @return
	 * @author Brian Risk
	 */
	private static int indexForCodonArray(char [] codon) {
		if (codon[0] == 'A') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 0;
				} else if (codon[2] == 'C') {
					return 1;
				} else if (codon[2] == 'G') {
					return 2;
				} else {
					return 3;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 4;
				} else if (codon[2] == 'C') {
					return 5;
				} else if (codon[2] == 'G') {
					return 6;
				} else {
					return 7;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 8;
				} else if (codon[2] == 'C') {
					return 9;
				} else if (codon[2] == 'G') {
					return 10;
				} else {
					return 11;
				}
			} else {
				if (codon[2] == 'A') {
					return 12;
				} else if (codon[2] == 'C') {
					return 13;
				} else if (codon[2] == 'G') {
					return 14;
				} else {
					return 15;
				}
			}
		} else if (codon[0] == 'C') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 16;
				} else if (codon[2] == 'C') {
					return 17;
				} else if (codon[2] == 'G') {
					return 18;
				} else {
					return 19;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 20;
				} else if (codon[2] == 'C') {
					return 21;
				} else if (codon[2] == 'G') {
					return 22;
				} else {
					return 23;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 24;
				} else if (codon[2] == 'C') {
					return 25;
				} else if (codon[2] == 'G') {
					return 26;
				} else {
					return 27;
				}
			} else {
				if (codon[2] == 'A') {
					return 28;
				} else if (codon[2] == 'C') {
					return 29;
				} else if (codon[2] == 'G') {
					return 30;
				} else {
					return 31;
				}
			}
		} else if (codon[0] == 'G') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 32;
				} else if (codon[2] == 'C') {
					return 33;
				} else if (codon[2] == 'G') {
					return 34;
				} else {
					return 35;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 36;
				} else if (codon[2] == 'C') {
					return 37;
				} else if (codon[2] == 'G') {
					return 38;
				} else {
					return 39;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 40;
				} else if (codon[2] == 'C') {
					return 41;
				} else if (codon[2] == 'G') {
					return 42;
				} else {
					return 43;
				}
			} else {
				if (codon[2] == 'A') {
					return 44;
				} else if (codon[2] == 'C') {
					return 45;
				} else if (codon[2] == 'G') {
					return 46;
				} else {
					return 47;
				}
			}
		} else if (codon[0] == 'T') {
			if (codon[1] == 'A') {
				if (codon[2] == 'A') {
					return 48;
				} else if (codon[2] == 'C') {
					return 49;
				} else if (codon[2] == 'G') {
					return 50;
				} else {
					return 51;
				}
			} else if (codon[1] == 'C') {
				if (codon[2] == 'A') {
					return 52;
				} else if (codon[2] == 'C') {
					return 53;
				} else if (codon[2] == 'G') {
					return 54;
				} else {
					return 55;
				}
			} else if (codon[1] == 'G') {
				if (codon[2] == 'A') {
					return 56;
				} else if (codon[2] == 'C') {
					return 57;
				} else if (codon[2] == 'G') {
					return 58;
				} else {
					return 59;
				}
			} else {
				if (codon[2] == 'A') {
					return 60;
				} else if (codon[2] == 'C') {
					return 61;
				} else if (codon[2] == 'G') {
					return 62;
				} else {
					return 63;
				}
			}
		} else {
			return -1; //return STOP
		}
	}
	
	
	/**
	 * addSeleno adds a Selenocysteine to this transcripts Seleno list.  This method also modifies the transcripts sequence to mark it with the selenocysteines location.
	 *  This method does not sort the Selenocysteines.
	 *  ****THIS METHOD DOES NOT WORK FOR MITOCHONDRIA****
	 */
	public void addSeleno(Seleno s){
		
		//Get the sequence ready to insert the special selenocysteine codon
		String before = transcriptSequence.substring(0, s.getStart() - getStart());
		String after = transcriptSequence.substring(s.getStop() + 1 - getStart(), transcriptLength);
		String selenoRef = transcriptSequence.substring(s.getStart() - getStart(), s.getStop() + 1 - getStart());

		//Errorcheck to ensure that a stop codon was found at the marked location.
		if (strand == Definitions.genomicStrandPOSITIVE){
			if(!selenoRef.equals("TGA")){
				return;
			}	
		}else if(strand == Definitions.genomicStrandNEGATIVE){
			if(!selenoRef.equals("TCA")){
				return;
			}
		}else{
			return;
		}

		//Convert the transcripts sequence with the selenocysteine
		transcriptSequence = before + s.getSequence() + after;
		
		//Add this seleno to the seleno objects list for use later in statistics
		Seleno.add(s);
	}
	
	/**
	 * isAStop determines if a 3 character string is a DNA stop codon.  Returns true if it is, false otherwise.  Considers mitochondrion DNA
	 */
	public boolean isAStop(String codon, int strand){
		
		if(!isMito){
			String stop1 = "TGA";
			String stop2 = "TAG";
			String stop3 = "TAA";
		
			if(strand == Definitions.genomicStrandNEGATIVE){
				stop1 = "TCA";
				stop2 = "CTA";
				stop3 = "TTA";
			}
				
			if(codon.equals(stop1) || codon.equals(stop2) || codon.equals(stop3)){
				return true;
			}
			
			
			return false;
		}else{
			String stop1 = "AGA";
			String stop2 = "TAG";
			String stop3 = "TAA";
			String stop4 = "AGG";
			
			if(strand == Definitions.genomicStrandNEGATIVE){
				stop1 = "TCT";
				stop2 = "CTA";
				stop3 = "TTA";
				stop4 = "CCT";
			}
				
			if(codon.equals(stop1) || codon.equals(stop2) || codon.equals(stop3) || codon.equals(stop4)){
				return true;
			}
			
			
			return false;
			
		}
	}
	
	/**
	 * getProteinWithSpliceMarkers returns this transcripts protein with special markers inserted at the splice junction locations of this protein.
	 * @return  Returns this proteins sequence with | inserted at splice junctions.
	 */
	public String getProteinWithSpliceMarkers(){
		StringBuffer workBench = new StringBuffer(this.protein);
		
		for(int k = exonSplitLocations.size() - 1; k > 0; k--){
			workBench.insert((int) (exonSplitLocations.get(k) /1), '|');
		}
		
		
		return workBench.toString();
	}
	
	/**
	 * isAStart determines if a 3 character string is a DNA start codon.  Returns true if it is a 3 character match, false otherwise.  Considers mitochondrion DNA
	 */
	public boolean isAStart(String codon, int strand){
		
		//Normal Genetic Code
		if(!isMito){
			//Positive strand
			String start = "ATG";
	
			//Negative strand
			if(strand == Definitions.genomicStrandNEGATIVE){
				start = "CAT";
			}
				
			if(codon.equals(start)){
				return true;
			}
			
			
			return false;
			//mitochondrion
		}else{
			//Positive Strand
			String start = "ATG";
			String start2 = "ATA";
			
			
			//Negative Strand
			if(strand == Definitions.genomicStrandNEGATIVE){
				start = "CAT";
				start2 = "TAT";
	
			}
				
			if(codon.equals(start) || codon.equals(start2)){
				return true;
			}
			
			
			return false;
			
			
		}
			
	}
	
	/**
	 * addCDS adds a CDS to this transcripts CDS list.  This method does not sort the CDS.
	 */
	public void addCDS(CDS c){
		//Add it to the end of the list if its place has not been found
		CDS.add(c);
	}
	
	/**
	 * @return iD is the index of this objects GENCODE_GTF_Line in the GTFlist.  Use the id to reference into the GTFlist to get all of the associated information for this transcript.
	 */
	public int getiD() {
		return iD;
	}


	/**
	 * @return The length of this transcript.
	 */
	public int gettranscriptLength() {
		return transcriptLength;
	}
	/**
	 * @return The previously determined score for this transcript.
	 */
	public double getAmpScore() {
		return ampScore;
	}
	/**
	 * @param Amplification Score to set for this transcript.
	 */
	public void setAmpScore(double AmplificationScore) {
		ampScore = AmplificationScore;
	}
	/**
	 * @return Returns the protein for this transcript.
	 */
	public String getProtein() {
		return protein;
	}


	/**
	 * @return The DNA sequence of this transcript.  This sequence is a raw/unmodified DNA encoded sequence.
	 */
	public String gettranscriptSequence() {
		return transcriptSequence;
	}



	/**
	 * @return An arraylist of CDS objects for this transcript.
	 */
	public ArrayList<CDS> getCDS() {
		return CDS;
	}

	/**
	 * @param Increments the variantCount
	 */
	public void incrementVariantCount() {
		variantCount++;
	}

	/**
	 * @return The number of variants for this transcript.
	 */
	public int getvariantCount() {
		return variantCount;
	}
	
	/**
	 * @return the extensionLength
	 */
	public int getEndExtensionLength() {
		return endExtensionLength;
	}
	/**
	 * @return the startExtensionLength
	 */
	public int getStartExtensionLength() {
		return startExtensionLength;
	}
	/**
	 * @return the containsPreMatureStop
	 */
	public boolean doesContainPreMatureStop() {
		return containPreMatureStop;
	}
	/**
	 * @return the startErrorLoc
	 */
	public int getStartErrorLoc() {
		return startErrorLoc;
	}

	/**
	 * @return the stopErrorLoc
	 */
	public int getStopErrorLoc() {
		return stopErrorLoc;
	}

	/**
	 * @return the preMatureStopCount
	 */
	public int getPreMatureStopCount() {
		return preMatureStopCount;
	}
	/**
	 * @return the seleno
	 */
	public ArrayList<Seleno> getSeleno() {
		return Seleno;
	}
	/**
	 * @return the leftOverAminoLength
	 */
	public int getLeftOverAminoLength() {
		return leftOverAminoLength;
	}
	/**
	 * @return the exonSplitLocations
	 */
	public ArrayList<Double> getExonSplitLocations() {
		return exonSplitLocations;
	}

	/**
	 * @return the hasStartCodon
	 */
	public boolean isHasStartCodon() {
		return hasStartCodon;
	}

	/**
	 * @param hasStartCodon the hasStartCodon to set
	 */
	public void setHasStartCodon(boolean hasStartCodon) {
		this.hasStartCodon = hasStartCodon;
	}
	/**
	 * @return the hasStopCodon
	 */
	public boolean isHasStopCodon() {
		return hasStopCodon;
	}

	/**
	 * @param hasStopCodon the hasStopCodon to set
	 */
	public void setHasStopCodon(boolean hasStopCodon) {
		this.hasStopCodon = hasStopCodon;
	}
}
