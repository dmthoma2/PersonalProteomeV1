package PeppyPeptideLocationIdentifer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.NoSuchElementException;
import java.util.Scanner;

import PeppyPeptideLocationIdentifer.Subclasses.PEPPY_RESULTS_Line_FORMAT16;
import PersonalProteome.Annotation;
import PersonalProteome.Definitions;
import PersonalProteome.GENCODE_GTF_Line;
import PersonalProteome.U;
import PersonalProteome.Gene.CDS;
import PersonalProteome.Gene.Transcript;


/**
 * PeptideLocationIdentifer is a tool that uses Personal Proteome to create a proteome to use to identify a peptides location in the genome.  It outputs the coordinates of the peptide as would be used in the UCSC browser.
 * It allows for splice junctions to be identified in the UCSC browser.
 *
 * *****FIX THE OFF BY ONE ERROR ON THOSE SELECT PEPTIDES*******
 * 
 * @author David "Corvette" Thomas
 *
 */
public class PeptideLocationIdentifer {

	
	/* MARKER FOR USING MITO*/
	boolean useMito = false;
	//File storage information
	String annotationFile;
	String chrmDir;
	String outputDir;
	
	//Information for writing out files
	private SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
	private SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	private Calendar cal;
	private String startTime = "";
	
	//Storage of Input
	ArrayList<GENCODE_GTF_Line> GTFlist = new ArrayList<GENCODE_GTF_Line>();
	ArrayList<Transcript> transcriptList = new ArrayList<Transcript>();
	HashMap<String, Transcript> transcriptHash = new HashMap<String, Transcript>();
	
	
	//Storage of peppy file/output.
	ArrayList<PEPPY_RESULTS_Line_FORMAT16> PEPPYList = new ArrayList<PEPPY_RESULTS_Line_FORMAT16>();
	ArrayList<String> PEPPYListEndings = new ArrayList<String>();
	ArrayList<String> outputList = new ArrayList<String>();
	ArrayList<String> bedList = new ArrayList<String>();
	
	
	/**
	 * PeptideLocationIdentifer is a tool that uses Personal Proteome to create a proteome to use to identify a peptides location in the genome.  It outputs the coordinates of the peptide as would be used in the UCSC browser.
	 * 
	 * @param annotationFile  The annotation file to use when creating a proteome.
	 * @param chrmDir  The directory of chromosome files.
	 * @param outputDir  The directory to place output files.
	 */
	public PeptideLocationIdentifer(String annotationFile, String chrmDir, String outputDir){
		//Store the needed input and output file locaitons.
		this.annotationFile = annotationFile;
		this.chrmDir = chrmDir;
		this.outputDir = outputDir;
	}
	
	/**
	 * uploadProteome Uses PP to upload an annotation and genome, and create a proteome from it.  This will then be used to get the location of peptides from a peppy results file.
	 * This method only preps this object to find peptide locations.  It only needs to be ran once per object.
	 */
	public void uploadProteome(){
		//Get the current time
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		
		U.p("Uploading Annoation and populating the genome directory.");
		U.startStopwatch();
		Annotation a = new Annotation(useMito);
		a.populateChrmArray(this.chrmDir);
		a.populateGTFList(new File(annotationFile));
		this.GTFlist = a.getAnnotaitonLines();
		U.stopStopwatch();
		
		U.p("Uploading Transcripts and Creating Proteins");
		U.startStopwatch();
		a.populateTranscriptsAndCDS(false);
		this.transcriptList = a.getTranscripts();
		U.stopStopwatch();
		
		U.p("Hashing transcripts.");
		U.startStopwatch();
		hashTranscriptList();
		U.stopStopwatch();
		//Clean up the a to save memory
		a = null;
	}
	
	
	/**
	 * This method should be called in succession for each peppy results file to modify.  Simply pass it the peppy results file and it will create a bed file and a appended peppy results file by the same name in the output directory.
	 * 
	 * @param peppyResultsFile File location of the peppy results file.
	 */
	public void identifyLocations(String peppyResultsFile){
		
		//Store the header information so that it can be modified and inserted into the output.
		String peppyFileHeader;
		
		//Instantiate needed variables to store peptide information
		PEPPYList = new ArrayList<PEPPY_RESULTS_Line_FORMAT16>();
		PEPPYListEndings = new ArrayList<String>();
		outputList = new ArrayList<String>();
		bedList = new ArrayList<String>();
		
		
		//Read in Peppy Results file
		U.p("Uploading peppy results file.");
		U.startStopwatch();
		peppyFileHeader = uploadPeppyResults(peppyResultsFile);
		U.stopStopwatch();
		
		//Remove the newLine character from the end of the peppyFileHeader to allow for the extra column headers to be inserted.
		peppyFileHeader = peppyFileHeader.substring(0, peppyFileHeader.length() - 1);
		
		
		U.p("Number of PEPPY Lines read in: " + PEPPYList.size());

		
		//Look up each peptides location, and mark it in the peptides Line
		U.p("Getting locations for each peptide.");
		U.startStopwatch();
		peppyFileHeader = identifyPeptideLocation(peppyFileHeader);
		U.stopStopwatch();
		//print each of the peptide line objects to the file.
		
		
		U.p("Creating output file.");
		createOutputFile(peppyResultsFile, peppyFileHeader);
		
		U.p("Number of peptides in the output file is " + outputList.size());
	
	}
	
	
	
	/**
	 * Creates a peppy output file identical to the original except that it has been appended with columns about location and strand.  Also create a bed file with each line formated for the UCSC browser.
	 * This file is formated with all 12 available columns from the UCSC browser.  Its label is the peptides sequence, and its score is the peptides ID from the peppy file.
	 * 
	 */
	private void createOutputFile(String peppyFile, String headerInfo){
		
		
		//Output file try block
		try {
			//Writing to HDD related variables
			StringBuffer sb;
			BufferedWriter out;
			File outFile;
			
			//Create a unqiue folder for hte output with a timestamp for a name.
			outFile = new File(outputDir + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");
			outFile.mkdir();
			String outputDir = outFile.getAbsolutePath() + "/";
		
			//Create a buffer and a writer for creating output with.
			out = new BufferedWriter(new FileWriter(outputDir + "RESULTS_" + (peppyFile.substring(peppyFile.lastIndexOf('/') + 1))));
			sb = new StringBuffer();

			//Write the file header to the buffer.
			sb.append(headerInfo + "\n");
			
			//Write each appended peppy line to the buffer.
			for(int i = 0; i < outputList.size(); i++){
				sb.append(outputList.get(i) + "\n");
			}
			
			//Write the data to the file
			out.write(sb.toString());
			
			//Flush the pipeline to ensure every part of the file has been written and close the writer.
			out.flush();
			out.close();
		
			//Create a buffer and a writer for creating output with.
			out = new BufferedWriter(new FileWriter(outputDir + "BEDFILE_" + (peppyFile.substring(peppyFile.lastIndexOf('/') + 1)) + ".bed"));
			sb = new StringBuffer();

			
			//Write each line of the bedList to the bed file.
			for(int i = 0; i < outputList.size(); i++){
				sb.append(bedList.get(i) + "\n");
			}
			
		
			//Write the data to the file
			out.write(sb.toString());
			
			//Flush the pipeline to ensure every part of the file has been written and close the writer.
			out.flush();
			out.close();
		
			
	    //Catch exceptions from output files
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}//createOutputFile
	
	/**
	 * identifies the location of each peptide, and creates the additional information to be stored at the end of each input line.  Returns the header with new columns appened to it.
	 * @return
	 */
	private String identifyPeptideLocation(String header){
		
		//Location to declare what strings this method adds to the peppy file.
		String fieldsToAdd = "\tChromosomeName\tStrand\tStartLocation\tStopLocation\tBlockCount\tBlockSize\tBlockStart";

		
		//Variables to write out in fieldsToAdd, and their default values.
		String chromosomeName = "";
		String strand = "";
		int startLoc = -1;
		int stopLoc = -1;
		int blockCount = -1;
		String blockSize = "";
		String blockStart = "";
		
		//Special marker used for setting which peptides to output.  Setting show = true anywhere in the code means that only peptides that reach that code will be used in the output.
		boolean show = false;

		
				
		//Iterate through the peptide list and add identify the locations of each peptide
		for( int i = 0; i < PEPPYList.size(); i++){
			PEPPY_RESULTS_Line_FORMAT16 currentPeppyLine = PEPPYList.get(i);
			
			//Identify which transcript this peptide belongs to
			String geneID = "";
			String transcriptID = PEPPYList.get(i).getProteinName();
			
			//pull the geneID outof this peppy line.
			 geneID = transcriptID.substring(transcriptID.indexOf('|') + 1);
			 geneID = geneID.substring(0, geneID.indexOf('|'));
			 //Get the transcript ID of this peptide
			 transcriptID = transcriptID.substring(0, transcriptID.indexOf('|'));
			 
			 
			 Transcript beingWorkedOn = transcriptHash.get(geneID + transcriptID);
			 
			 //A valid transcript has been matched up to a peptide.
			 if(beingWorkedOn != null){
				 
				 //The genomic strand that this peptide belongs on
				 int currentStrand = GTFlist.get(beingWorkedOn.getiD()).getGenomicStrand();

				 //Peptide locations within this protein
				 int peptideLocationInProtein = currentPeppyLine.getStart();
				 int peptideStopLocationInProtein = peptideLocationInProtein + currentPeppyLine.getPeptideSequence().length();
				 
				 //Collection of splice Junctions found within this peptide
				 ArrayList<Integer> spliceJuncWithinPeptide = new ArrayList<Integer>();
				 

				 int frameShift = Definitions.genomicPhaseZERO;
				 if(beingWorkedOn.getProtein().contains("X")){

					 frameShift = beingWorkedOn.getCDS().get(0).getGenomicPhase();
					 if(currentStrand == Definitions.genomicStrandNEGATIVE){
						frameShift =  beingWorkedOn.getCDS().get(beingWorkedOn.getCDS().size() - 1).getGenomicPhase();
					 }
						frameShift = 3 - frameShift;

				 }
				 
					if(currentPeppyLine.getID() == 96955){
						U.p("Found!");
					}
				 if(peptideLocationInProtein != -1){
				 

					 
					 int splitLoc = -1;
					 for(int k = 0; k < beingWorkedOn.getExonSplitLocations().size(); k++){
						 
//						 U.p("Split Location " + k + " is: " + beingWorkedOn.getExonSplitLocations().get(k));
						 if((beingWorkedOn.getExonSplitLocations().get(k)) > peptideLocationInProtein){
							 //Determine if a splice junction is found within a peptide
							 if(beingWorkedOn.getExonSplitLocations().get(k) < (peptideStopLocationInProtein)){
								 spliceJuncWithinPeptide.add(k);

							 }//exon is less then peptide stop
							 continue;
						 }//exon is greater then peptide start

						 
						 splitLoc++;
					 }
					 int currentCDS = splitLoc + 1;
					 if(currentStrand == Definitions.genomicStrandNEGATIVE){
						 currentCDS = (beingWorkedOn.getCDS().size() - 1) - currentCDS;

					 }


					 
					 

				
					 //This section of Code calculates
					//Start location, and stop location.
					 double lengthInCDS = 0;
					 
					//Determine how many bp the peptide goes into the given CDS
				 	if(splitLoc == -1){
				 		lengthInCDS =(peptideLocationInProtein)*3;
				 	}else{
				 		lengthInCDS = ((peptideLocationInProtein) - beingWorkedOn.getExonSplitLocations().get(splitLoc))*3;
				 	}
					 
					 //Deal with truncated numbers becuase comp does not recognize 69.999999999999 = 70
					 int floor = (int)lengthInCDS;
					 double ceiling = floor + 1;
				     double  minValue = floor + .99;
				     
				     //Clean up the length region if it needs it.
				     if(lengthInCDS < ceiling && lengthInCDS >= minValue ){
				    	 lengthInCDS = ceiling;
				     }
					 	 
					 //Start location for no splits and forward strand
					 int peptideStartLoc = beingWorkedOn.getCDS().get(currentCDS).getStart() + ((int)lengthInCDS);
					 int peptideStopLoc = peptideStartLoc + (currentPeppyLine.getPeptideSequence().length() * 3);

					 /*Start and stop identified*/
					 startLoc = peptideStartLoc - 1 - frameShift;
					 stopLoc = peptideStopLoc  - 1 - frameShift;

					 //Look for negative locations
					 if(currentStrand == Definitions.genomicStrandNEGATIVE){
						 peptideStopLoc = beingWorkedOn.getCDS().get(currentCDS).getStop() - ((int)lengthInCDS);
						 peptideStartLoc = peptideStopLoc - (currentPeppyLine.getPeptideSequence().length() * 3);
						 
						 
						 /*Start and stop identified*/
						 startLoc = peptideStartLoc + frameShift;
						 stopLoc = peptideStopLoc  + frameShift;
					 }
					 
					 blockCount = 1 + spliceJuncWithinPeptide.size();
						 
					 //Handle the case where there are no exons splits
					 if(spliceJuncWithinPeptide.size() == 0){
						 //Set the Block size for peptides without a split anywhere

						 blockSize = String.valueOf(peptideStopLoc - peptideStartLoc);
						 blockStart = "0";
					//Else for No Exon Split Locations
					//This section of code handles setting blockSize, start and stop locations for Peptides that span multiple exons.
					 }else{
						 int genomicStrand = GTFlist.get(beingWorkedOn.getiD()).getGenomicStrand();
						 
						 //debug


							 
							 
						 boolean spliceMessUp = false;

						 int inCDS = currentCDS;
						 ArrayList<CDS> cdsList =  beingWorkedOn.getCDS();
						 int distanceLeft = currentPeppyLine.getPeptideSequence().length() * 3;
						 
						 //List of start locations
						 ArrayList<Integer> startLocs = new ArrayList<Integer>();
						 
						 //List of stop Locations
						 ArrayList<Integer> stopLocs = new ArrayList<Integer>();
						 
						 
						 
						 //Single split positive
						 if(genomicStrand == Definitions.genomicStrandPOSITIVE){
							 startLocs.add(startLoc);
							 stopLocs.add(cdsList.get(inCDS).getStop());
							 
							 distanceLeft -= (stopLocs.get(0) - startLocs.get(0));
							 
							 inCDS++;
							 
							 if(distanceLeft > 0){
								 //**********************************
								 //Add in information to handle multiple jumps
								 for(int u = 1; u < spliceJuncWithinPeptide.size(); u++){
									 startLocs.add(cdsList.get(inCDS).getStart() - 1);
									 stopLocs.add(cdsList.get(inCDS).getStop());
									 
									 distanceLeft -= (stopLocs.get(u) - startLocs.get(u));
									 inCDS++;
								 }
								 
								 if(cdsList.size() - 1 >= inCDS){
									 if(distanceLeft > 0){
										 startLocs.add(cdsList.get(inCDS).getStart() - 1);
										 stopLocs.add(cdsList.get(inCDS).getStart() - 1 + distanceLeft);
									 }else{
										 blockCount--;
										 spliceMessUp = true;
									 }
								 }else{
									 U.p("CUT OFF EARLY ERROR CHECK");
								 }
							 }else{
								 blockCount--;
								 spliceMessUp = true;
							 }
						//Work on the reverse strand
						 }else{
						 
							 stopLocs.add(stopLoc);
							 startLocs.add(cdsList.get(inCDS).getStart()  - 1);
							 
							 distanceLeft -= stopLocs.get(0) - startLocs.get(0);
							 
							 
								 inCDS--;
							if(distanceLeft > 0){
								 for(int u = 1; u < spliceJuncWithinPeptide.size(); u++){
									 stopLocs.add(cdsList.get(inCDS).getStop());
									 startLocs.add(cdsList.get(inCDS).getStart() - 1);
									 
									 distanceLeft -= stopLocs.get(u) - startLocs.get(u);
									 inCDS--;
								 }
								 
								 //**********************************
								 //Add in information to handle multiple jumps
								 
								 if(0 <= inCDS){
									 
									 if(distanceLeft > 0){
										stopLocs.add(cdsList.get(inCDS).getStop());
									 	startLocs.add(cdsList.get(inCDS).getStop() - distanceLeft);
									 }else{
										 blockCount--;
										 spliceMessUp = true;
									 }
								 }
							 }else{
//								 startLocs.remove(0);
//								 startLocs.add(cdsList.get(inCDS + 1).getStart());
								 blockCount--;
								 spliceMessUp = true;
							 }
							 Collections.reverse(stopLocs);
							 Collections.reverse(startLocs);			 
							 
						 }//Negative strand
						 

						 
						//Start locations, and stop locations.
						 startLoc = startLocs.get(0);
						 stopLoc = stopLocs.get(stopLocs.size() - 1);
						 //Block size and Block start.  We always start at index 0, and the count is simply 1 + number of junctions are spanned by the peptide.
						 
						 
						 
						 //DEBUG
//						 if(blockCount == 2){
//							 for(int z = 0; z < startLocs.size(); z++){
//								 if(stopLocs.get(z) - startLocs.get(z) == 0){
//									 
//								 }
//							 }
//							 
//						 }
//						 if(show){
//							 U.p("Splice Junctions are");
//							 for(int r = 0; r < spliceJuncWithinPeptide.size(); r++){
//								 U.p("Splice Junc " + beingWorkedOn.getExonSplitLocations().get(spliceJuncWithinPeptide.get(r)));
//							 }
//							 U.p("Line being worked on is: " + currentPeppyLine.toString());
//							 U.p("inCDS is: " + inCDS);
//							 U.p("CDS regions are");
//							 for(int e = 0; e < beingWorkedOn.getCDS().size(); e++){
//								 U.p("CDS " + e + " is: " + beingWorkedOn.getCDS().get(e).toString());
//							 }
//						 }

						//set the block sizes
						//set block starts
						 int rollingStartLoc = 0;
						 for(int y = 0; y < startLocs.size(); y++){
							 int size = -1;
							 
							 size = stopLocs.get(y) - startLocs.get(y);

//							 if(spliceMessUp){
//								 size--;
//							 }
							 rollingStartLoc = startLocs.get(y) - startLocs.get(0);
							 blockSize += size;
							 blockSize += ",";
							 blockStart += rollingStartLoc;
							 blockStart += ",";
							 
//							 if(size == 0){

//							 	if(show){
//								 U.p("Current Stop Location: " + stopLocs.get(y));
//								 U.p("Current Start Location: " + startLocs.get(y));
//							 	}

//							 }
							 
							 
						 }//for loop
						 
 
					 }//exon Split size > 0
					 
					 
 
					 
					 
					 /*Strand Identified*/
					 if(GTFlist.get(beingWorkedOn.getiD()).getGenomicStrand() == Definitions.genomicStrandPOSITIVE){
						 strand = "+";
					 }else{
						 strand = "-";
					 }
					 
					 /*ChromosomeName Identified*/
					 chromosomeName = Definitions.convertChrmNumToString(GTFlist.get(beingWorkedOn.getiD()).getChromosomeName());
					 
					 


//					 if(blockCount == 1){
//						if(Integer.valueOf(blockSize) != currentPeppyLine.getPeptideSequence().length() * 3){
//							U.p("Peptide size does not match block size!");
//							U.p("Peptide: " + currentPeppyLine.toString());
//						}
//					 }
					 
				 }//Valid peptide found in protein
				 
//				 if(blockCount == 2){
//
//				 }
				
				 //show = true here beings that a valid peptide was found.  Output all valid peptides.
				 show = true;
			 }//beingWorkedOn != null   Valid transcript found

			//Save that information
			String bedLine = chromosomeName  + "\t" + startLoc + "\t" + stopLoc + "\t" + currentPeppyLine.getPeptideSequence() + "\t" + currentPeppyLine.getID() + "\t" + strand +  "\t" + startLoc + "\t" + stopLoc + "\t" + "0" + "\t"+ blockCount + "\t" + blockSize + "\t" + blockStart;
				 
			String ending = chromosomeName  + "\t" + strand + "\t" + startLoc + "\t" + stopLoc + "\t" + blockCount + "\t" + blockSize + "\t" + blockStart;
			//Devug

			PEPPYListEndings.add(ending);
			
			chromosomeName = "";
			strand = "";
			startLoc = -1;
			stopLoc = -1;
			blockCount = -1;
			blockSize = "";
			blockStart = "";

			
			if(show){
				outputList.add(PEPPYList.get(i).toString() + PEPPYListEndings.get(i));
				bedList.add(bedLine);
			}
			show = false;

		}//Outmost for loop
		
		
		return header + fieldsToAdd;
	}
	
	/**
	 * hashTranscriptList creates a hashtable used for lookup by peptides after they are read in.  This creates 
	 * Storage based on geneID + transcriptID
	 */
	private void hashTranscriptList(){
		for(Transcript tran: this.transcriptList){
			if(!tran.getProtein().equals("Protein Too Short")){
				transcriptHash.put(GTFlist.get(tran.getiD()).getGeneID() + GTFlist.get(tran.getiD()).getTranscriptID(), tran);
			}
		}//for
		
	}
	
	/**
	 * uploadPeppyResults uploads a peppy results file with modifications and stores it in the PEPPYList.
	 * @param fileLocation  The location of the peppy output file to upload and use.
	 * @return  The header lines of this file
	 */
	private String uploadPeppyResults(String fileLocation){
		
		//DEBUGGING
		//Only read in lines up to the breakAt point.  This is useful for debugging.
		int breakAt = -1;
		int count = 0;
		
		String header = "";
		try {
			Scanner s = new Scanner(new File(fileLocation));
			String token;
			
			for(int i = 0; i < 4; i++){
				token = s.nextLine();
				header += token;
				header += "\n";
			}
			
			
			//Parse information for each line
			int ID = -1;
			int SpectrumID = -1;
			String MD5 = "";
			String FileName = "";
			double score = -1;
			double PrecursorMZ = -1;
			double PrecursorNeutralMass = -1;
			double eValue = -1;
			String peptideSequence = "";
			int start = -1;
			int stop = -1;
			String proteinName = "";
			int matchRank = -1;
			int RankCount = -1;
			int IonCount = -1;
			boolean labeled = false;
			int charge = -1; 
			int cleavageAcidCount = -1;
			boolean inORF = false;
			double hydrophobic = -1;
			double hydrophilic = -1;
			boolean isModified = false;
			double modMass = 0;
			int modIndex = -1;

			/*Parse each line*/
			while(s.hasNextLine()){
				
				//Once a line without any content on it is reached, consider the file over
				if(!s.hasNext()){
					break;
				}
				//Get the ID
				token = s.next();
				ID = Integer.parseInt(token);
				
				//Get the SpectrumID
				token = s.next();
				SpectrumID = Integer.parseInt(token);
				
				//Get the MD5
				token = s.next();
				MD5 = token;
				
				//Get the fileName
				token = s.next();
				FileName = token;
				
				
				//Get the score
				token = s.next();
				score = Double.parseDouble(token);
				
				//Get the PrecursorMZ
				token = s.next();
				PrecursorMZ = Double.parseDouble(token);
				
				//Get the PrecursorNeutralMass
				token = s.next();
				PrecursorNeutralMass = Double.parseDouble(token);
				
				//Get the eValue
				token = s.next();
				eValue = Double.parseDouble(token);
				
				//Get the peptideSequence
				token = s.next();
				peptideSequence = token;
				
				//Get the start
				token = s.next();
				start = Integer.parseInt(token);
				
				//Get the stop
				token = s.next();
				stop = Integer.parseInt(token);
				
				//Get the protein Name
				token = s.next();
				proteinName = token;
				
				//Get the matchRank
				token = s.next();
				matchRank = Integer.parseInt(token);
				
				//Get the RankCount
				token = s.next();
				RankCount = Integer.parseInt(token);
				
				//Get the IonCount
				token = s.next();
				IonCount = Integer.parseInt(token);
				
				
				//Get the labeled column value
				token = s.next();
				labeled = Boolean.valueOf(token);
				
				//Get the charge
				token = s.next();
				charge = Integer.parseInt(token);
				
				//Get the cleavage acid count
				token = s.next();
				cleavageAcidCount = Integer.parseInt(token);
				
				//Get the inORF
				token = s.next();
				inORF = Boolean.valueOf(token);
				
				//Get the hydrophobic
				token = s.next();
				hydrophobic = Double.valueOf(token);
				
				//Get the Hydrophilic
				token = s.next();
				hydrophilic = Double.valueOf(token);
				
				//get the isModified
				token = s.next();
				isModified = Boolean.valueOf(token);
				
				//get the modMass
				token = s.next();
				modMass = Double.valueOf(token);
				
				
				//get the mod Index
				token = s.next();
				modIndex = Integer.valueOf(token);
				
				s.nextLine();
				
				//List of the peppy input file lines
				PEPPYList.add(new PEPPY_RESULTS_Line_FORMAT16(ID, SpectrumID, MD5, FileName, score, PrecursorMZ, PrecursorNeutralMass, eValue, peptideSequence, start, stop, proteinName, matchRank, RankCount, IonCount, labeled, charge,
						        cleavageAcidCount, inORF, hydrophobic, hydrophilic, isModified, modMass, modIndex));
				
				//Reset the default values
				ID = -1;
				SpectrumID = -1;
				MD5 = "";
				FileName = "";
				score = -1;
				PrecursorMZ = -1;
				PrecursorNeutralMass = -1;
				eValue = -1;
				peptideSequence = "";
				start = -1;
				stop = -1;
				proteinName = "";
				matchRank = -1;
				RankCount = -1;
				IonCount = -1;
				labeled = false;
				charge = -1; 
				cleavageAcidCount = -1;
				inORF = false;
				hydrophobic = -1;
				hydrophilic = -1;
				isModified = false;
				modMass = 0;
				modIndex = -1;
				
				//DEBUGGING BREAK
				if(count == breakAt){
					U.p("Breaking early in the peptide list, with a break early value of " + breakAt);
					break;
				}
				
				count++;
			}
			
			
			
			
			
			
			
			
		//File is not found.
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		//Throw an error when a unparsable token is read in.
		} catch (NoSuchElementException e) {
			U.p("Invalid file.  Unable to completely parse input.  Output will not contain all peptides.");
		}	
		
		//Return the header information from the file.
		return header;
	}//uploadPeppyResults
	
}//PeptideLocationIdentifer
