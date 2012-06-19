package PersonalProteome.PeptideAnalysisTool;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;

import Peppy.Sequence_DNA;
import PersonalProteome.Definitions;
import PersonalProteome.GENCODE_GTF_Line;
import PersonalProteome.U;
import PersonalProteome.Gene.CDS;
import PersonalProteome.Gene.Seleno;
import PersonalProteome.Gene.Transcript;
import PersonalProteome.PeptideAnalysisTool.SubClasses.Peptide;
import PersonalProteome.PeptideAnalysisTool.SubClasses.PeptideLocationLine;
/**
 * PeptideExonSpanningCheck is a tool that takes in a Peppy results file, pulls out the peptides and determines which peptides span an intron and how many.  It computes statistics and outputs that information
 * and the peptides to a specified output file.
 * 
 * PeptideAnalysisTool also can label peptides as intronic, exonic, transcriptic, intergenic, etc based on its location.  This is done using a BED File.
 * 
 * 
 * Spanning - Splice Junctions
 * Location - Label: Intronic, Exonic, etc
 * **NOTE FOR SPLICE JUNCTION CHECKING: IF A PEPTIDE OCCURS IN MORE THEN ONE LOCATION IN A PROTEIN THEN ONLY THE FIRST ONE IS CONSIDERED.  THIS IS A AMBIGUITY THAT CAN ONLY BE SOLVED IF CORRECT START/STOP LOCATIONS ARE AVAILABLE***
 * @author David "Corvette" Thomas
 *
 */
public class PeptideAnalysisTool{

	
	//Universal
	//This Arraylist contains all of the information from the annotation.
	ArrayList<GENCODE_GTF_Line> GTFlist;
	
	
	
	//Splice junction check
	//This Arraylist contains all of the peptides from the peptides file.
	ArrayList<Peptide> peptideList;
	//A list of transcripts to produce proteins
	ArrayList<Transcript> transcriptList;
	//An arrayList to store all peptides that span a exon
	ArrayList<Peptide> spanningPeptides = new ArrayList<Peptide>();
	
	//Locaton check
	ArrayList<PeptideLocationLine> peptideLocationList  = new ArrayList<PeptideLocationLine>();
	
	
	private boolean genomeHasMito;
	
	//File locations
	String annotationFile;
	String peptideListFile;
	String outputDir;
	String outputFileName;
	String chrmDir;
	String[] chrmFile = new String[25];
	
	
	//Time/Formating
	SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
	SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	Calendar cal;
	String startTime = "";
	String endTime = "";

	/**
	 * CONSTRUCTOR FOR PEPTIDE LOCATION IDENTIFICATION!!! This constructor takes in a Peptide file (Format(BEDFILE): chr4	179830523	179830547	DKPGCPQR	1000	-	209.814621	300.8298934	1	1) and determines
	 * whether each peptide is intergenic, exonic, intronic, etc.  
	 * @param annotationFile  The annotation file to use in determing peptide location.
	 * @param peptideListFile  The peptide file containing the peptides to checked locations on.
	 * @param outputDir  The directory to put output in.
	 */
	public PeptideAnalysisTool(String annotationFile,String peptideListFile, String outputDir){
		//Set default values
		this.annotationFile = annotationFile;
		this.peptideListFile = peptideListFile;
		this.outputDir = outputDir;
		
	}

	/**
	 * 
	 * CONSTRUCTOR FOR EXONS SPANNING CHECK!!! This constructor takes in a Peppy results file, pulls out the peptides and determines which peptides span a splice junction and how many.  It computes statistics and outputs that information
     * and the peptides to a specified output file.
     * 
	 * @param annotationFile  GENCODE annotation file containing all of the transcripts to be checked.
	 * @param peptideListFile Peppy results file containing information about peptides and matches.
	 * @param outputDir Directory to output the results file.
	 * @param outputFileName  File name to give to the output file.
	 * @param chrmDir Directory where the chromosome files are located.
	 */
	public PeptideAnalysisTool(String annotationFile,String peptideListFile, String outputDir, String outputFileName, String chrmDir){
		//Set default values
		this.annotationFile = annotationFile;
		this.peptideListFile = peptideListFile;
		this.outputDir = outputDir;
		this.outputFileName = outputFileName;
		this.chrmDir = chrmDir;
		
		File mitoFile = new File(chrmDir + "chrM.fa");
		genomeHasMito = mitoFile.exists();
		
		//Populate the chromosome files
		populateChrmArray(chrmDir);
	}
	 
	
	/**
	 * checkPeptideLocation uses the input files to create a list of transcripts/exons and peptides.  It then labels each peptide based on its location in the annotation.
	 * Finally checkPeptideLocation creates a set of output files in the specified output directory.  This method should only be called after its appropriate constructor.
	 */
	public void checkPeptideLocation(){
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		
		//Populate GTFlist
		U.p("Uploading annotation.");
		U.startStopwatch();
		populateGTFList(true);
		U.stopStopwatch();
		
		
		//Populate peptides list
		U.p("Populating the peptidesList.");
		U.startStopwatch();
		populatePeptidesList();
		U.stopStopwatch();
		
		U.p("Total number of protein coding Genes/Exons int annotation: " + GTFlist.size());
		U.p("Total number of peptides to check location on is: " + peptideLocationList.size());
		
		//Populate peptides list
		U.p("Determining the location of each peptide.");
		U.startStopwatch();
		labelPeptidesBasedOnLocation();
		U.stopStopwatch();
		
		//Create Location output
		U.p("Creating output based on location");
		U.startStopwatch();
		createLocationOutput();
		U.stopStopwatch();
		
		
		//Get the ending time
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.p("Finsihing up: " + endTime);
	}
	

	/**
	 * populatePeptidesList loads in Peptides from a file with peptide location and score information, and creates PeptideLocationLines from the file.
	 * This list is stored in the List peptideLocationList
	 */
	private void populatePeptidesList(){
		File peptideFile = new File(peptideListFile);
		//
		int chromosomeName = -1;
		int startLocation = -1;
		int stopLocation = -1;
		String sequence = "";
		String score = "";
		int strand = -1;
		String restOfLine = "";
		//Parse each line of the file and save it
		try{
			Scanner s = new Scanner(peptideFile);
			String token;

			//Parse each Line, one at a time
			while(s.hasNextLine()){
				 token = s.next();
				
				 /*Plug values from the fields array into variables*/
				if(token.equalsIgnoreCase("chrM")){
					chromosomeName = Definitions.chromosomeM;
				}else if(token.equalsIgnoreCase("chrX")){
					chromosomeName = Definitions.chromosomeX;
				}else if(token.equalsIgnoreCase("chrY")){
					chromosomeName = Definitions.chromosomeY;
				}else{
					chromosomeName = Integer.parseInt(token.substring(token.indexOf('r') + 1));
					
				}
					
				
				token = s.next();
				//Get the start location
				startLocation = Integer.parseInt(token);
				
				token = s.next();
				//Get the stop location
				stopLocation = Integer.parseInt(token);
				 
				//Get the amino acid sequence from the file
				token = s.next();
				sequence = token;
				
				//Get the first score from the file
				token = s.next();
				score = token;
				
				
				//Determine genomic strand
				token = s.next();
				if(token.equalsIgnoreCase("+")){
					strand = Definitions.genomicStrandPOSITIVE;
				}else {
					strand = Definitions.genomicStrandNEGATIVE;
				}
				
				//Preserve the rest of the data from the input file
				token = s.nextLine();
				restOfLine = token;
				
				//Add the newly read line into the peptideLocationList
				peptideLocationList.add(new PeptideLocationLine(chromosomeName, startLocation, stopLocation, sequence, score, strand, restOfLine));
				
				 //Reset all variables
				chromosomeName = -1;
				startLocation = -1;
				stopLocation = -1;
				sequence = "";
				score = "";
				strand = -1;
				restOfLine = "";
				
			}//while

			
		}catch(FileNotFoundException e){
			U.p("Error populating Peptide Location List: " + e);

		}//catch
	}


	
	
	/**
	 * labelPeptidesBasedOnLocation labels each peptide based on its location.  It will label a peptides start and stop regions, and determine if the start and stop occur in the same region.
	 * labelPeptidesBasedOnLocation used a custom algorithm for comparisons, and is completely memory in-place.  n = # of lines read in from the annotation file, m = # of peptides read from the input.
	 * labelPeptidesBasedOnLocation runs in O(mlogm + m) time, with memory usage of ½(1).
	 * A hash of the annotation could result in a faster runtime if m is on the same order of magnitude or greater then n.  It runs in O(n + m), but would require extra memory on the order of ½(n).  
	 * This method checks every transcript that a peptide could occur in and selects the most specific result. EX. EXONIC > INTRONIC > TRANSCRIPTIC > INTERGENIC
	 * Future Improvement: replace the iterator with a for loop and direct object manipulation.  Without this implementation memory is ½(m).
	 */
	private void labelPeptidesBasedOnLocation(){
		
		//Sort the peptide list.  Uses a merge sort and sorts in mlog(m) time.
		Collections.sort(peptideLocationList);
		
		ArrayList<PeptideLocationLine> tempOut = new ArrayList<PeptideLocationLine>();
		//Iterate through the peptides, moving the through the annotation as each peptide is checked.
		int annIndex = 0;
		for(int i = 0; i < 25; i++){
			U.p("Working on Chrm: " + (i + 1));

			
			//Iterate through the peptide list
			for(PeptideLocationLine p: peptideLocationList){
				//Ignore peptides that are not on the current chromosome
				if(p.getChromosomeName() != (i+1)){
					continue;
				}

				
				while(p.getChromosomeName() > GTFlist.get(annIndex).getChromosomeName()){
					annIndex++;
				}

				//Loop through the annotation until a unit is found that the peptide is in or beyond.
				for(annIndex = annIndex;p.getChromosomeName() == GTFlist.get(annIndex).getChromosomeName() && p.getStartLocation() > GTFlist.get(annIndex).getStopLocation(); annIndex++){
				}

				
				String geneID = GTFlist.get(annIndex).getGeneID();
				ArrayList<Integer> transList = new ArrayList<Integer>();
				transList.add(annIndex);
				//Add the additional transcirpts to the transList
				int temp = 1;
				//Search for other transcripts in this gene that the peptide occurs atleast partly in
				while(geneID.equals(GTFlist.get(annIndex+temp).getGeneID())){
					if(GTFlist.get(annIndex+temp).getFeatureType() == Definitions.featureTypeTRANSCRIPT){
						if(p.getStartLocation() < GTFlist.get(annIndex+temp).getStopLocation()){
							transList.add(annIndex+temp);
						}//Inner if
					}//Outer If
					
					temp++;
				}//While Loop

				ArrayList<PeptideLocationLine> pepPossibleClassifications = new ArrayList<PeptideLocationLine>();
				

				//Iterate through each possible transcript this peptide could occur in, and classify it.
				for(int k = 0; k < transList.size(); k++){
					int transID = transList.get(k);
					int startID = transID;
					int stopID = transID;
					
					//A gene is always identified
					//Potentially creep forward to determine if a Exon is more accurate then a transcript
					int creepIndex = transID + 1;
					//Check to determine 
					while(GTFlist.get(creepIndex).getFeatureType() != Definitions.featureTypeTRANSCRIPT){
						
						//If a peptide has its start within any exon
						if(p.getStartLocation() >= GTFlist.get(creepIndex).getStartLocation() && p.getStartLocation() <= GTFlist.get(creepIndex).getStopLocation()){
							startID = creepIndex;
						}
						//If a peptide has its stop within any exon
						if(p.getStopLocation() >= GTFlist.get(creepIndex).getStartLocation() && p.getStopLocation() <= GTFlist.get(creepIndex).getStopLocation()){
							stopID = creepIndex;
						}
						
						
						creepIndex++;
					}

					//Determine how many exons this transcript contains
					int testIndex = transID + 1;
					int exonsCount = 0;
					while(GTFlist.get(testIndex).getFeatureType() != Definitions.featureTypeTRANSCRIPT){
							exonsCount++;
							testIndex++;
					}
					
					
					GENCODE_GTF_Line transcriptUnit = GTFlist.get(transID);
					GENCODE_GTF_Line startUnit = GTFlist.get(startID);
					GENCODE_GTF_Line stopUnit = GTFlist.get(stopID);
					
					//Check situations where the start and stop occur inside or before a transcript
					if(transID == startID && transID == stopID){
						//Start and stop occur before the transcript
						if(p.getStartLocation() < transcriptUnit.getStartLocation() && p.getStopLocation() < transcriptUnit.getStartLocation()){
							p.setStartLocType(Definitions.LOCATION_TYPE_INTERGENIC);
							p.setStopLocType(Definitions.LOCATION_TYPE_INTERGENIC);
							p.setSameSubUnit(true);
						}
						//Start occurs before transcript, stop occurs inside of transcript
						if(p.getStartLocation() < transcriptUnit.getStartLocation() && p.getStopLocation() >= transcriptUnit.getStartLocation()){
							p.setStartLocType(Definitions.LOCATION_TYPE_INTERGENIC);
							p.setStopLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
							p.setSameSubUnit(false);
						}
						//Start occurs in transcript, stop occurs after transcript
						if(p.getStartLocation() >= transcriptUnit.getStartLocation() && p.getStopLocation() > transcriptUnit.getStopLocation()){
							p.setStartLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
							p.setStopLocType(Definitions.LOCATION_TYPE_INTERGENIC);
							p.setSameSubUnit(false);
						}
						//Start and stop occur inside of transcript
						if(p.getStartLocation() >= transcriptUnit.getStartLocation() && p.getStopLocation() <= transcriptUnit.getStopLocation()){
							//Since both occur in the transcript, then the peptide must be in the transcript or a intron, possibly spanning the two.
							
							//Default value for the transcript
							p.setStartLocType(Definitions.LOCATION_TYPE_INTRON);
							p.setStopLocType(Definitions.LOCATION_TYPE_INTRON);
							p.setSameSubUnit(true);
							

			
							//Determine whether either the start or stop occur in the transcript.
							if(p.getStrand() == Definitions.genomicStrandPOSITIVE){
								//if peptide occurs before first EXON
								int first = GTFlist.get(transID + 1).getStartLocation();
								if(p.getStartLocation() < first ){
									p.setStartLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
									
								}
								if(p.getStopLocation() < first){
									p.setStopLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
								}
								int last = GTFlist.get(transID + exonsCount).getStopLocation();
								if(p.getStartLocation() > last){
									p.setStartLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
									
								}
								if(p.getStopLocation() > last){
									p.setStopLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
								}
								//or after last exon, it is transcriptic
							//Negative
							}else{
								int first = GTFlist.get(transID + 1).getStopLocation();
								if(p.getStartLocation() > first && p.getStopLocation() > first){
									p.setStartLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
									p.setStopLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
								}
								int last = GTFlist.get(transID + exonsCount).getStartLocation();
								if(p.getStartLocation() < last && p.getStopLocation() < last){
									p.setStartLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
									p.setStopLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
								}
							}//Negative
							
							
						}//Start/Stop occur inside of transcript
						
						
					}//tran tran
					
					//Check the situation where the start and stop both occur in an exon
					if(startUnit.getFeatureType() == Definitions.featureTypeEXON && stopUnit.getFeatureType() == Definitions.featureTypeEXON){
						p.setStartLocType(Definitions.LOCATION_TYPE_EXON);
						p.setStopLocType(Definitions.LOCATION_TYPE_EXON);
						
						if(startID == stopID){
							p.setSameSubUnit(true);
						}else{
							p.setSameSubUnit(false);
						}
						
					}//exon exon
					//Consider the situation when a start is in a transcript, and the stop is in an exon.
					if(startID == transID && stopUnit.getFeatureType() == Definitions.featureTypeEXON){
						
						//Start occurs before the transcript, end occurs in an exon
						if(p.getStartLocation() < transcriptUnit.getStartLocation()){
							p.setStartLocType(Definitions.LOCATION_TYPE_INTERGENIC);
	
						}else if((transID + 1) == stopID){
							p.setStartLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
	
						}else{
						
						//At this point, the stop must occur in a exon that is not the first one.  Since the start is in the transcript, it must occur in an intron
						p.setStartLocType(Definitions.LOCATION_TYPE_INTRON);
					
						}
						p.setStopLocType(Definitions.LOCATION_TYPE_EXON);
						p.setSameSubUnit(false);
						
						
					}//tran exon

					
					//Consider the situation when the start is in an exon, and the stop is in a transcript
					if(startUnit.getFeatureType() == Definitions.featureTypeEXON && stopID == transID){
						
						//Stop occurs after the transcript ends
						if(p.getStopLocation() > transcriptUnit.getStopLocation()){
							p.setStopLocType(Definitions.LOCATION_TYPE_INTERGENIC);
							
							//Stop occurs after the exons
						}else if(p.getStopLocation() > GTFlist.get(transID + exonsCount).getStopLocation()){
								p.setStopLocType(Definitions.LOCATION_TYPE_TRANSCRIPT);
						}else{
								p.setStopLocType(Definitions.LOCATION_TYPE_INTRON);
						}
						
						
						
						//Otherwise stop occurs within an intron
						p.setStartLocType(Definitions.LOCATION_TYPE_EXON);
						p.setSameSubUnit(false);
					}//exon tran
					
					if(GTFlist.get(transID).getGenomicStrand() != p.getStrand()){
						p.setStartLocType(Definitions.LOCATION_TYPE_INTERGENIC);
						p.setStopLocType(Definitions.LOCATION_TYPE_INTERGENIC);
						p.setSameSubUnit(true);
					}
					
					
					//Create a copy of p to store in the possibilites list
					PeptideLocationLine placeHolder = new PeptideLocationLine(p.getChromosomeName(), p.getStartLocation(), p.getStopLocation(), p.getSequence(), p.getScore(), p.getStrand(), p.getRestOfLine());
					placeHolder.setStartLocType(p.getStartLocType());
					placeHolder.setStopLocType(p.getStopLocType());
					placeHolder.setSameSubUnit(p.isSameSubUnit());
					
					//Add a peptide to the possibilities list.
					pepPossibleClassifications.add(placeHolder);
					//Return this peptide to its default values
					p.setSameSubUnit(true);
					p.setStartLocType(Definitions.LOCATION_TYPE_UNKNOWN);
					p.setStopLocType(Definitions.LOCATION_TYPE_UNKNOWN);
				}//trans loop

				//Determine which p to keep
				PeptideLocationLine bestP = pepPossibleClassifications.get(0);
			
				//Determine which pep has the best possible score
				for(int n = 1; n < pepPossibleClassifications.size();n++){
					PeptideLocationLine currentPep = pepPossibleClassifications.get(n);
					int topScore = bestP.getStartLocType() + bestP.getStopLocType();
					
					if((currentPep.getStartLocType() + currentPep.getStopLocType()) > topScore){
						bestP = currentPep;
					}
				}
			
				tempOut.add(bestP);

				//After the peptide has been determined, move the peptides 
			}//for iterator through peptide list
		}//for loop chromosome files
		
		
		//Save all of the modified peptides
		peptideLocationList = tempOut;
		
		//Checking for unidentified peptides
		U.p("Error printing.  Any peptide with a unidentified start or stop location type will be printed here.");
		for(PeptideLocationLine p: peptideLocationList){
			if(p.getStartLocType() == Definitions.LOCATION_TYPE_UNKNOWN || p.getStopLocType() == Definitions.LOCATION_TYPE_UNKNOWN){
				U.p("Peptide: " +   PeptideLocationLine.convertLocTypeToString(p.getStartLocType()) + " " + PeptideLocationLine.convertLocTypeToString(p.getStopLocType()) + " sameSub " + p.isSameSubUnit());
				U.p(p.toString());
			}

		}
		U.p("Total number of peptides identified: " + peptideLocationList.size());
	}//labelPeptidesBasedOnLocation
	
	
	
	
	

	

	/**
	 * createLocationOuput creates a set of files containing the various classifications of peptides, as well as a file containing information about what files were used and statistics of the run.
	 */
	private void createLocationOutput() {
		
		//Writing to HDD related variables
		StringBuffer sb;
		BufferedWriter out;
		File outFile;
		
		//File names corresponds to the pepLists Arary
		ArrayList<String> FileNames = new ArrayList<String>();
		FileNames.add("IntergenicPeps.txt");
		FileNames.add("ExonicPeps.txt");
		FileNames.add("IntronicPeps.txt");
		FileNames.add("IntronExonBoundaryPeps.txt");
		FileNames.add("TranscripticPeps.txt");
		FileNames.add("AllOthers.txt");
		
		ArrayList<ArrayList<PeptideLocationLine>> pepLists = new ArrayList<ArrayList<PeptideLocationLine>>();
		
		ArrayList<PeptideLocationLine> intergenicPeps = new ArrayList<PeptideLocationLine>();
		ArrayList<PeptideLocationLine> exonicPeps = new ArrayList<PeptideLocationLine>();
		ArrayList<PeptideLocationLine> intronicPeps = new ArrayList<PeptideLocationLine>();
		ArrayList<PeptideLocationLine> intronExonPeps = new ArrayList<PeptideLocationLine>();
		ArrayList<PeptideLocationLine> transcriptPeps = new ArrayList<PeptideLocationLine>();
		ArrayList<PeptideLocationLine> allOthers = new ArrayList<PeptideLocationLine>();
		
		//Populate the lists
		for(PeptideLocationLine p: peptideLocationList){
			if(p.getStartLocType()  == Definitions.LOCATION_TYPE_INTERGENIC && p.getStopLocType() == Definitions.LOCATION_TYPE_INTERGENIC && p.isSameSubUnit()){
				intergenicPeps.add(p);
			}else if(p.getStartLocType() == Definitions.LOCATION_TYPE_EXON && p.getStopLocType() == Definitions.LOCATION_TYPE_EXON && p.isSameSubUnit()){
				exonicPeps.add(p);
			}else if(p.getStartLocType() == Definitions.LOCATION_TYPE_TRANSCRIPT && p.getStopLocType() == Definitions.LOCATION_TYPE_TRANSCRIPT && p.isSameSubUnit()){
				transcriptPeps.add(p);
			}else if(p.getStartLocType() == Definitions.LOCATION_TYPE_INTRON && p.getStopLocType() == Definitions.LOCATION_TYPE_INTRON && p.isSameSubUnit()){
				intronicPeps.add(p);
			}else if(p.getStartLocType() == Definitions.LOCATION_TYPE_INTRON && p.getStopLocType() == Definitions.LOCATION_TYPE_EXON || p.getStartLocType() == Definitions.LOCATION_TYPE_EXON && p.getStopLocType() == Definitions.LOCATION_TYPE_INTRON ){
				intronExonPeps.add(p);
				
			}else{
				allOthers.add(p);
			}
			
			
		}
		
		//Put the individual lists into the pepLists
		pepLists.add(intergenicPeps);
		pepLists.add(exonicPeps);
		pepLists.add(intronicPeps);
		pepLists.add(intronExonPeps);
		pepLists.add(transcriptPeps);
		pepLists.add(allOthers);
		
		try{
		outFile = new File( outputDir + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");
		outFile.mkdir();
		String outputDir = outFile.getAbsolutePath() + "/";
		
		for(int i = 0; i < pepLists.size(); i++){
			sb = new StringBuffer();
			out = new BufferedWriter(new FileWriter(outputDir + FileNames.get(i)));
			
			//Write the file contents to the string buffer
			sb.append("//Number of Peptides: " + pepLists.get(i).size() + "\n");
			for(PeptideLocationLine pepToWrite: pepLists.get(i)){
				sb.append(pepToWrite.toString() + "\n");
				sb.append("Start: " + PeptideLocationLine.convertLocTypeToString(pepToWrite.getStartLocType()) + "\tStop: " + PeptideLocationLine.convertLocTypeToString(pepToWrite.getStopLocType()) + "\t SameSubUnit: " + pepToWrite.isSameSubUnit() + "\n");

			}
			
			
			//Write the file and close it out.
			out.write(sb.toString());
			out.flush();
			out.close();
		}
		//Write out the stats File
		sb = new StringBuffer();
		out = new BufferedWriter(new FileWriter(outputDir + "Statistics.txt"));
		sb.append("//This output file was written at: " + sdf.format(cal.getTime()) + "\n");
		sb.append("//Annotation file used: " + annotationFile + "\n");
		sb.append("//Peptides file used: " + peptideListFile + "\n");
		sb.append("//Total # of peptides labeled: " + peptideLocationList.size() + "\n");
		sb.append("//Total # of annotation lines used: " + GTFlist.size() + "\n");
		sb.append("//************************" + "\n");
		sb.append("//Total # of intergenic peptides: " + intergenicPeps.size() + "\n");
		sb.append("//Total # of exonic peptides: " + exonicPeps.size() + "\n");
		sb.append("//Total # of intronic peptides: " + intronicPeps.size() + "\n");
		sb.append("//Total # of intron/exon boundry peptides: " + intronExonPeps.size() + "\n");
		sb.append("//Total # of transcriptic peptides: " + transcriptPeps.size() + "\n");
		sb.append("//Total # of misc. peptides: " + allOthers.size() + "\n");
			
		
		//Write the file and close it out.
		out.write(sb.toString());
		out.flush();
		out.close();
		
		}catch (IOException e){
			e.printStackTrace();
		}
		
	}//createLocationOutput
	
	/**
	 * checkPeptideSpanning is the primary action performed by PeptideExonSpanningCheck.  It displays information to the default output, and drives all of the gears inside of PeptideExonSpanningCheck.
	 * It handles loading in input, processing the data, and outputing the results. This method should only be called after its appropriate constructor.
	 */
	public void checkPeptidesSpanning(){
		
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		
		//Populate GTFlist
		U.p("Uploading annotation.");
		U.startStopwatch();
		populateGTFList(false);
		U.stopStopwatch();
		
		//Populate peptides list
		U.p("Populating the peptidesList.");
		U.startStopwatch();
		populatePeptidesListPeppyFile();
		U.stopStopwatch();
		
		//Get transcripts/proteins
		U.p("Adding transcripts to transcriptList.");
		U.startStopwatch();
		populateTranscriptsAndCDS();
		U.stopStopwatch();
		

		U.p("Comparing peptides and transcripts.");
		U.p("Total number of peptides: " + peptideList.size());
		U.p("Total number of transcripts: " + transcriptList.size());
		
		//Compare peptides to their transcripts
		U.startStopwatch();
		determineIfPeptidesAreOnSplit();
		U.stopStopwatch();
		U.p("spanning Peptides size " + spanningPeptides.size());
		
		//Output the final list of peptides that are a match
		U.p("Outputing the found peptides");
		U.startStopwatch();
		outputPeptidesSpliceJunction();
		U.stopStopwatch();
		
		
		//Get the ending time
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.p("Finsihing up: " + endTime);
	
	}
	/**
	 * Populate the peptideList with peptide objects.  This assumes a peppy output file was used.
	 * This method assumes that a transcriptID + geneID are a unique identifier for a transcript.
	 */
	private void populatePeptidesListPeppyFile(){
		File peptideFile = new File(peptideListFile);
		//
		peptideList = new ArrayList<Peptide>();
		String transcriptID = "";
		String sequence = "";
		String geneID = "";
		
		//Parse each line of the file and save it
		try{
			Scanner s = new Scanner(peptideFile);
			String token;
			//Parse each Line, one at a time
			while(s.hasNextLine()){
				 token = s.next();

				 try{
				 Integer.parseInt(token);
				 }catch (NumberFormatException e){
					 s.nextLine();
					 continue;
				 }
				 //On a valid line, so grab the sequence and transcript id
				 //Skip the garbage ahead of the sequence
				 for(int i = 0; i < 6; i++){
					 s.next();
				 }
				 
				 //Get the sequence of acids for this peptide
				 sequence = s.next();
				 
				 //Skip this peptides start and stop
				 s.next();
				 s.next();

				 transcriptID = s.next();

				 geneID = transcriptID.substring(transcriptID.indexOf('|') + 1);
				 geneID = geneID.substring(0, geneID.indexOf('|'));
				 //Get the transcript ID of this peptide
				 transcriptID = transcriptID.substring(0, transcriptID.indexOf('|'));
				 
				 s.nextLine();
				 
				 peptideList.add(new Peptide(transcriptID, geneID, sequence));

				transcriptID = "";
				sequence = "";
				geneID = "";
			}//while

			
		}catch(FileNotFoundException e){
			U.p("Error populating PeppyFile Peptide List: " + e);

		}//catch
		
		
	}//populatePeptidesListPeppyFile
	/**
	 * Creates a master list of Transcripts and CDS.  Each transcripts creates it protein after it is made.  The transcripts are stored in transcirptList.
	 * It creates transcripts based on the files loaded into the chromosome directory.  
	 */
	private void populateTranscriptsAndCDS(){
		boolean modStopStart = false;
		
		transcriptList = new ArrayList<Transcript>();
		//Variables for the identification loop
		int transcriptCount = 0;
		int i = 0;
		//A list of indices in the GTFlist of transcript objects to lookup and create
		ArrayList<Integer> transcripts = new ArrayList<Integer>();
		while(i < GTFlist.size()){
			if(GTFlist.get(i).getFeatureType() == Definitions.featureTypeTRANSCRIPT){
				transcriptCount++;
				/*Add the index of the GTFLine to create a transcript from*/
				transcripts.add(i);
				
			}//if
			i++;
		}//while

		//Variables for the addition to the transcript list loop
		int addCount = 0;
		Transcript t;	
		Sequence_DNA seq;
		//Iterate through each chromosome file, and create all of the transcripts possible from each one
		//Code is performed in this order to minize HDD reads.
		for(int k = 0; k < 25; k++){
			//Ignore the mitochondrian DNA
			if(k == 22){
				if(!genomeHasMito){
					continue;
				}
			}
			
			U.p("Working on Chrm File: " + (k + 1));

			seq = new Sequence_DNA(chrmFile[k]);
			GENCODE_GTF_Line beingTested;
			
			//Check each transcript to determine if it is using the currently open file.
			for(int j = 0; j < transcripts.size(); j++){
				
				//Get the currently used transcripts line
				beingTested = GTFlist.get(transcripts.get(j));
				
				//Check to see if a transcript is in the current file
				if(beingTested.getChromosomeName() - 1 == k){
					
					//sequence is the DNA strand from the chromosome file for this transcript
					String sequence;
					if(k != 22){
						sequence = seq.getNucleotideSequences().get(0).getSequence().substring(beingTested.getStartLocation() - 1,beingTested.getStopLocation());
					}else{
						sequence = seq.getNucleotideSequences().get(0).getSequence().substring(beingTested.getStartLocation(),beingTested.getStopLocation() + 1);
					}
					//Construct a basic transcript to be added to the list
					t = new Transcript(beingTested.getChromosomeName(), sequence, beingTested.getGenomicStrand() , transcripts.get(j), beingTested.getStartLocation(), beingTested.getStopLocation(), modStopStart);
					

					//Prime the loop
					//c simply stores each line as it is iterated through
					GENCODE_GTF_Line c = GTFlist.get(transcripts.get(j) + 1);
					//Let every transcript through since we are looking for peptides
					boolean startFound = true;
					
					for(int x = transcripts.get(j) + 1; x < GTFlist.size() && c.getFeatureType() != Definitions.featureTypeTRANSCRIPT; x++){
						c = GTFlist.get(x);

						
						//If a start codon is found, then set startFound to true
						if(c.getFeatureType() == Definitions.featureTypeSTART_CODON){
							startFound = true;
						}
						
						//If a selenocysteine is found, then add it to the the transcripts seleno list
						if(c.getFeatureType() == Definitions.featureTypeSELENOCYSTEINE){
							Seleno s = new Seleno(c.getStartLocation(), c.getStopLocation());
							t.addSeleno(s);
						}
						//if a CDS is found, add it to this transcripts CDS list
						if(c.getFeatureType() == Definitions.featureTypeCDS){
							CDS temp = new CDS(c.getStartLocation(), c.getStopLocation(), c.getGenomicPhase());
							t.addCDS(temp);

						}
						
						
					}
					//Confirm that the transcript has a start and has a CDS to code from
					if(startFound && t.getCDS().size() > 0){

						
					
						//Create t's protein and add it to the list
						t.createProtein();
						
						//Now that t has a score, has CDS regions, and has been confirmed to have a start, add it to the list of transcripts to produce proteins with
						transcriptList.add(t);
						addCount++;
					}


					
				}//if
				
			}//inner for
			//Quit the loop early when all of the transcripts have been found
			
		}//outer for

		
	}//populateTranscriptsAndCDS
	
	/**
	 * determineIfPeptidesAreOnSplit finds iterates through the peptides fetches each peptides transcript.  The transcripts are checked to contain the peptide, and if they do the number of introns
	 * spanned by the peptide is totaled.  Each peptide that spans at least one intron is added to the spanningPeptides ArrayList.
	 */
	private void determineIfPeptidesAreOnSplit(){
		HashMap<String, Transcript> transcriptHash = new HashMap<String, Transcript>(transcriptList.size() + 10);
		
		for(Transcript tran: transcriptList){
			if(!tran.getProtein().equals("Protein Too Short")){
				transcriptHash.put(GTFlist.get(tran.getiD()).getGeneID() + GTFlist.get(tran.getiD()).getTranscriptID(), tran);
			}
		}//for
		Transcript t;
		//Iterate through each peptide and find its transcript
		for(Peptide pep: peptideList){

			//Add all peptides that span a splice junction to the spanning peptides list
			t = transcriptHash.get(pep.getGeneID() + pep.getTranscriptID());
			if(t != null){
				pep.determineIfOnSplit(t);
				if(pep.getSpanCount() > 0){
					spanningPeptides.add(pep);
				}
			}

		}//for
	}//determineIfPeptidesAreOnSplit
	/**
	 * outputPeptides outputs the peptide list to the disk.  The output file consists of several lines of comments (marked with // at the start) containing stats/information about the run,
	 * followed by the peptides that span an intron.  The format for these pepties is:
	 * transcriptID|geneID peptideSequence #ofIntronsSpanned
	 */
	private void outputPeptidesSpliceJunction(){
		BufferedWriter out;
		File outFileTime;
		
		try {
			
			//Construct a folder for the output
			outFileTime = new File(outputDir + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");
			outFileTime.mkdir();
			String modOutputDir = outFileTime.getAbsolutePath() + "/";
			out = new BufferedWriter(new FileWriter(modOutputDir + outputFileName));
			
			//Determine what the maximum number of introns spanned by any peptide.
			int maxIntronsSpanned = 0;
			for(Peptide p: spanningPeptides){
				if(p.getSpanCount() > maxIntronsSpanned){
					maxIntronsSpanned = p.getSpanCount();
				}
			}
			

			//Write output to the file

			out.write("//This output file was written at: " + sdf.format(cal.getTime()));
			out.newLine();
			out.write("//Total number of transcripts found in annotation: " + transcriptList.size());
			out.newLine();
			out.write("//Annotation file used " + annotationFile);
			out.newLine();
			out.write("//Genome file used " + chrmDir);
			out.newLine();
			out.write("//Total number of peptides read in from peptide file: " + peptideList.size());
			out.newLine();
			out.write("//Peptide file used " + peptideListFile);
			out.newLine();
			out.write("//Total number of peptides that were found to span a intron: " + spanningPeptides.size());
			out.newLine();
			out.write("//Percentage of peptides that spanned a intron: " + (((double)spanningPeptides.size()/(double)peptideList.size())*100) + "%");
			out.newLine();
			out.write("//The maximum number of introns spanned by a peptide was: " + maxIntronsSpanned);
			out.newLine();
			out.write("transcriptID|geneID NumberOfIntronsSpanned sequence spliceJuncLocations protein ");
			out.newLine();
			for(Peptide p: spanningPeptides){

				String spliceJuncs = "";
				for(int k = 0; k < p.getSpliceJuncs().size(); k++){
					DecimalFormat df = new DecimalFormat("#.##");
					spliceJuncs += df.format(p.getSpliceJuncs().get(k).doubleValue());
					if(k + 1 != p.getSpliceJuncs().size()){
						spliceJuncs += ",";
					}
				}
				
				out.write(p.getTranscriptID() + "|" + p.getGeneID() + " "+ p.getSpanCount() + " " + spliceJuncs +  " " + p.getSequence() + " " );
				out.newLine();
			}
			
			
			//Flush and close the stream to ensure that all data is written to disk.
			out.flush();
			out.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	
	}//outputPeptidesSpliceJunction
	
	/**
	 * Populates the chromosome array with the file locations of each chromosome.
	 * @param chrmDir The directory containing the chromosome files.
	 */
	public void populateChrmArray(String chrmDir){
		for(int i = 0; i < 22; i++){
			chrmFile[i] = chrmDir + "chr" + (i + 1) + ".fa";
		}
		chrmFile[22] = chrmDir + "chrM.fa";
		chrmFile[23] = chrmDir + "chrX.fa";
		chrmFile[24] = chrmDir + "chrY.fa";
	}
	/**
	 * parseGTFFile takes a GTF format annotation file and populates an ArrayList, (GTFlist), of GENCODE_GTF_Line objects representing each line of the file.
	 * @param onlyTranscriptsExons is a boolean variable that only allows for transcripts and exons of protein encoding genes to be stored.
	 */
	private void populateGTFList(boolean onlyTranscriptExons){
		File GTFFile = new File(annotationFile);
		//GTFint List to work on
		GTFlist = new ArrayList<GENCODE_GTF_Line>();
		//Represents a single line
		GENCODE_GTF_Line line;
		//Working string
		String temp = null;
		//Fields of the line
		String[] mandField;
		String[] optionalField;
		//Variables for the mandatory fields
		int id = 0;
		int chromosomeName;
		int annotationSource = -1;
		int featureType = -1;
		int startLocation = -1;
		int stopLocation = -1;
		int score = -1;
		int genomicStrand = -1;
		int genomicPhase = -1;
		String geneID;
		String transcriptID;
		String gene_Type;
		int gene_Status = -1;
		String gene_Name;
		String transcript_Type;
		int transcript_Status = -1;
		String transcript_Name;
		int level = -1;
		boolean skip = false;
		boolean finished = false;
		int count = 1;
		
		//Parse each line of the file and save it
		try{
			Scanner s = new Scanner(GTFFile);

			//Parse each Line, one at a time
			while(s.hasNextLine()){
				
				if(!skip){
				temp = s.next();
				}

				//Move past commments
				if(temp.startsWith("##")){
					s.nextLine();
					
					continue;
				}
				
				/*Get Mandatory fields*/
				//Initialize variables
				mandField = new String[17];
				//Parse out the information from the line and store it in the mandFields array
				for(int i = 0; i < 8; i++){
					mandField[i] = temp;
					if(s.hasNext()){
					temp = s.next();
					}else{
						finished = true;
						if(finished){
							break;
						}
					}
				}
				if(finished){
					break;
				}
				
				/*Populate 8-16 with values*/
				for(int i = 8; i < 16; i++){
					temp = s.next();
					mandField[i] = temp.substring(temp.indexOf('"') +1 , temp.lastIndexOf('"'));
					s.next();
				}
				
				//Current token is level after the prvious loop, so add level in and start checking for optional parameters
				temp = s.next();
				mandField[16] = temp.substring(0, temp.length() - 1);
				
				
				/*Search for optional fields*/
				/*0 = tag, 1 = ccdsid, 2 = havana_gene, 3 = havana_transcript, 4 = ont*/
				optionalField = new String[5];
				
				
				//Iterate through and handle optional fields
				while(s.hasNext()){
					temp = s.next();

					if(temp.equals("tag")){
						temp = s.next();
						optionalField[0] = temp.substring(temp.indexOf('"') + 1, temp.lastIndexOf('"'));
					}else if(temp.equals("ccdsid")){
						temp = s.next();
						optionalField[1] = temp.substring(temp.indexOf('"') + 1, temp.lastIndexOf('"'));
					}else if(temp.equals("havana_gene")){
						temp = s.next();
						optionalField[2] = temp.substring(temp.indexOf('"') + 1, temp.lastIndexOf('"'));
					}else if(temp.equals("havana_transcript")){
						temp = s.next();
						optionalField[3] = temp.substring(temp.indexOf('"') + 1, temp.lastIndexOf('"'));
					}else if(temp.equals("ont")){
						temp = s.next();
						optionalField[4] = temp.substring(temp.indexOf('"') + 1, temp.lastIndexOf('"'));
					}else{
						break;
					}
					
				}//while
				skip = true;
				
				/*Plug values from the fields array into variables*/
				if(mandField[0].equalsIgnoreCase("chrM")){
					chromosomeName = Definitions.chromosomeM;
				}else if(mandField[0].equalsIgnoreCase("chrX")){
					chromosomeName = Definitions.chromosomeX;
				}else if(mandField[0].equalsIgnoreCase("chrY")){
					chromosomeName = Definitions.chromosomeY;
				}else{
					chromosomeName = Integer.parseInt(mandField[0].substring(mandField[0].indexOf('r') + 1));
					
				}
				
				

				if(mandField[1].equalsIgnoreCase("ENSEMBL")){
					annotationSource = Definitions.ENSEMBL;
				}else{
					annotationSource = Definitions.HAVANA;
				}
				
				if(mandField[2].equalsIgnoreCase("gene")){
					featureType = Definitions.featureTypeGENE;
				}else if(mandField[2].equalsIgnoreCase("transcript")){
					featureType = Definitions.featureTypeTRANSCRIPT;
				}else if(mandField[2].equalsIgnoreCase("CDS")){
					featureType = Definitions.featureTypeCDS;
				}else if(mandField[2].equalsIgnoreCase("start_codon")){
					featureType = Definitions.featureTypeSTART_CODON;
				}else if(mandField[2].equalsIgnoreCase("stop_codon")){
					featureType = Definitions.featureTypeSTOP_CODON;
				}else if(mandField[2].equalsIgnoreCase("exon")){
					featureType = Definitions.featureTypeEXON;
				}else if(mandField[2].equalsIgnoreCase("UTR")){
					featureType = Definitions.featureTypeUTR;
				}else if(mandField[2].equalsIgnoreCase("Selenocysteine")){
					featureType = Definitions.featureTypeSELENOCYSTEINE;
				}
				

				
				startLocation = Integer.parseInt(mandField[3]);
				stopLocation = Integer.parseInt(mandField[4]);

				/*Score is not used, so a place holder of -1 is put in*/
				score = -1;
				
				if(mandField[6].equalsIgnoreCase("+")){
					genomicStrand = Definitions.genomicStrandPOSITIVE;
				}else {
					genomicStrand = Definitions.genomicStrandNEGATIVE;
				}
				
				if(mandField[7].equalsIgnoreCase("0")){
					genomicPhase = Definitions.genomicPhaseZERO;
				}else if(mandField[7].equalsIgnoreCase("1")){
					genomicPhase = Definitions.genomicPhaseONE;
				}else if(mandField[7].equalsIgnoreCase("2")){
					genomicPhase = Definitions.genomicPhaseTWO;
				}else if(mandField[7].equalsIgnoreCase(".")){
					genomicPhase = Definitions.genomicPhasePERIOD;
					
				}

				geneID = mandField[8];
				transcriptID = mandField[9];
				gene_Type = mandField[10];
				
				if(mandField[11].equalsIgnoreCase("KNOWN")){
					gene_Status = Definitions.geneStatusKNOWN;
				}else if(mandField[11].equalsIgnoreCase("NULL")){
					gene_Status = Definitions.geneStatusNULL;
				}else if(mandField[11].equalsIgnoreCase("NOVEL")){
					gene_Status = Definitions.geneStatusNOVEL;
				}else if(mandField[11].equalsIgnoreCase("UNKNOWN")){
					gene_Status = Definitions.geneStatusUNKNOWN;
				}
				
				gene_Name = mandField[12].trim();
				transcript_Type = mandField[13].trim();
				
				if(mandField[14].trim().equalsIgnoreCase("KNOWN")){
					transcript_Status = Definitions.transcriptStatusKNOWN;
				}else if(mandField[14].equalsIgnoreCase("NULL")){
					transcript_Status = Definitions.transcriptStatusNULL;
				}else if(mandField[14].equalsIgnoreCase("NOVEL")){
					transcript_Status = Definitions.transcriptStatusNOVEL;
				}else if(mandField[14].equalsIgnoreCase("UNKNOWN")){
					transcript_Status = Definitions.transcriptStatusUNKNOWN;
				}

				transcript_Name = mandField[15];
				if(mandField[16].equalsIgnoreCase("1")){
					level = Definitions.GENCODELevelONE;
				}else if(mandField[16].equalsIgnoreCase("2")){
					level = Definitions.GENCODELevelTWO;
				}else if(mandField[16].equalsIgnoreCase("3")){
					level = Definitions.GENCODELevelTHREE;
				}
						

				line = new GENCODE_GTF_Line(id, chromosomeName,  annotationSource,  featureType,  startLocation,  stopLocation,  score,  genomicStrand,  genomicPhase,
						 geneID,  transcriptID,  gene_Type,  gene_Status,  gene_Name,  transcript_Type,  transcript_Status,
						 transcript_Name,  level,  optionalField[0],  optionalField[1],  optionalField[2],  optionalField[3],  optionalField[4]);
				
				
				/*Since we only care about encoding proteins, only add ones that are protein_coding*/
				/*Add the line object to the list of parsed lines*/
				if(transcript_Type.equalsIgnoreCase("protein_coding") || transcript_Type.equalsIgnoreCase("nonsense_mediated_decay")){
					/*Finished with variables, now create a GENCODE_GTF_Line object and insert it into the ArrayList*/
					if(!onlyTranscriptExons){
						GTFlist.add(line);
						id++;
					}else{
						//Only genes are exons are important for location information
						if(line.getFeatureType() == Definitions.featureTypeTRANSCRIPT || line.getFeatureType() == Definitions.featureTypeEXON){
							GTFlist.add(line);
							id++;
						}
						
					}
				}
				//reset all of the variables
				 
				 chromosomeName = -1;
				 annotationSource = -1;
				 featureType = -1;
				 startLocation = -1;
				 stopLocation = -1;
				 score = -1;
				 genomicStrand = -1;
				 genomicPhase = -1;
				 geneID = "";
				 transcriptID = "";
				 gene_Type = "";
				 gene_Status = -1;
				 gene_Name = "";
				 transcript_Type = "";
				 transcript_Status = -1;
				 transcript_Name = "";
				 level = -1;
				 count++;
			}//while

			
		}catch(FileNotFoundException e){
			U.p("Error populating GTF List: " + e);

		}//catch

		
		
	}//populateGTFList
}//Peptide Analysis tool
