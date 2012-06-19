package GENCODESubmissionTool;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Scanner;


import GENCODESubmissionTool.Subclasses.PeppyLine;
import PersonalProteome.Definitions;
import PersonalProteome.U;


/**
 * GENCODE Submission Tool is a tool designed to take Peppy input files, and create bed files for submission to GENCODE.
 * @author David "Corvette" Thomas
 *
 */
public class GENCODESubmissionTool {

	
	
	//I/O Processing Variables
	private Calendar cal;
	String inputFile;
	String outputDir;
	String FDRFile;
	
	//FDR Values
	double zero_fdr;
	double five_fdr;
	double ten_fdr;
	
	//Match Rank Filter
	int maxMatchRank = 1;
	
	
	//Peppy File Storage
	String header = "";
	ArrayList<PeppyLine> PeppyList = new ArrayList<PeppyLine>();
	ArrayList<Integer> bedScore = new ArrayList<Integer>();
	ArrayList<PeppyLine> matchRankFiltered = new ArrayList<PeppyLine>();
	boolean isModified = false;
	
	//Time/Output Variables
	private SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	
	
	
	
		/**
		 * GENCODESubmissionTool is a tool designed to take Peppy input files, and create bed files for submission to GENCODE.
		 * @param inputFile
		 * @param FDRFile
		 * @param outputDir
		 * @param maxMatchRank
		 * @param currentTime
		 */
	public GENCODESubmissionTool(String inputFile, String FDRFile, String outputDir, int maxMatchRank, Calendar currentTime){
		this.inputFile = inputFile;
		this.FDRFile = FDRFile;
		this.outputDir = outputDir;
		this.maxMatchRank = maxMatchRank;
		this.cal = currentTime;
		
//		Parse();
	}
	
	
	public void setInputFile(String inputFile){
		this.inputFile = inputFile;
	}
	public void Parse(){
		//Parse the peppyFile
		//Only keep lines with a filtered match count < User Specified
		parsePeppyFile();
	}
	
	public void Run(){

		
		
		
		//Assign unique keys to each peptide
		U.p("Assigning keys to each uploaded peptide");
		assignUniqueKey();
		
		//Remove peptides from the list that have the same Unique Key and Spectrum as one previously in the list
		removeDuplicates();

		
		U.p("PeppyList size: " + PeppyList.size());
		//Filter input
		funnel();

		//Get FDR Results
		getFDRResults();
		U.p("FDR Results: 0% " + zero_fdr + " 5% " + five_fdr + " 10% " + ten_fdr);

		
		//Get the location of all peptides in the genome
		determinePeptideCount();
		
		
		//Calculate BED score for each peppy line
		calcBedScore();
		//This is now done when creating bedFiles
		
		//Create bed files
		createBedFiles();
	}//Run
	
	/**
	 * 
	 */
	private void assignUniqueKey(){
		//Assign Unique keys to each peptide
		for(PeppyLine pl: PeppyList){
			//Key is of the format chr#:start-stop:strand
			//Example is     chr1:100-203:+
			String uniqueKey = "";
			//Set the chromosome name based on which type of file was passed in.
			if(!pl.isDna()){
				uniqueKey += Definitions.convertChrmNumToString(pl.getChromosomeName()) + ":";
				uniqueKey += pl.getStartLocation() + "-" + pl.getStopLocation() + ":";
			}else{
				uniqueKey += pl.getSequenceDescription() + ":";
				uniqueKey += pl.getStart() + "-" + pl.getStop() + ":";
				
			}//else
			
			//Add the start and stop to the uniqueKey
	//		uniqueKey += pl.getStart() + "-" + pl.getStop() + ":";
			
			//Determine the strand character
			String strand = Definitions.convertStrandToString(pl.getFirstStrand());
			
			//Finish the uniqueKey by adding the strand
			uniqueKey += strand;
			pl.setUniqueKey(uniqueKey);
		}
	}//assignUniqueKey
	
	/**
	 * Funnel narrows down the input lines.
	 */
	private void funnel(){
		//Go through the peptide list, and remove duplicates
		U.p("PeppyList "  + PeppyList.size());
		ArrayList<PeppyLine> masterDelete = new ArrayList<PeppyLine>();
		for(PeppyLine pl: PeppyList){
			//If a DNA
			if(pl.isDna()){

				//Since removing duplicates before funnel tool, then there will only be one match.
				ArrayList<PeppyLine> matches = new ArrayList<PeppyLine>();
				//Check for a match in the ppList
				for(PeppyLine pp: PeppyList){
					if(pp.getID() == pl.getID()){
						continue;
					}
					//If a match is found, add it to the matches list.
					if(pl.getMD5().equals(pp.getMD5())){
						matches.add(pp);
//						U.p("pl " + pl.toString());
//						U.p("pp " + pp.toString());
					}//matches check
				}//Iterate through possible pep
				
				//Continue if no matches are found
				if(matches.size() == 0){
					continue;
				}
				if(matches.size() >= 1){
//					U.p("Match size: " + matches.size());
				}
				
				
				
				//Now matches contains a list of all peptides that map to a single spectra.  
				matches.add(pl);
				
				//Iterate through the mathces, and find the greatest score within the set.
				
				double score = 0;
				for(int i = 0; i < matches.size(); i++){
					if(matches.get(i).getScore() > score){
						score = matches.get(i).getScore();
					}

				}
//				U.p(score);
				
				//Determine which matches to deleted
				ArrayList<PeppyLine> toDelete = new ArrayList<PeppyLine>();
				for(PeppyLine p: matches){
					if(p.getScore() < score){
						toDelete.add(p);
					}
				}
				for(PeppyLine td: toDelete){
					matches.remove(td);
				}
				
				//Now determine if there is a mix of types.
				
				boolean mix = false;
				PeppyLine first = matches.get(0);
				for(int j = 1; j < matches.size(); j++){

					if(matches.get(j).isDna() != first.isDna()){
						mix = true;
						break;
					}//if
				}//for
				
				//Label the lines as both if they are
				if(mix == true){
					for(PeppyLine p: matches){
						p.setBoth(true);
					}//for
				}//if
				
				//Remove the lower scored peptides
				for(PeppyLine tr: toDelete){
					masterDelete.add(tr);
				}//for
		
			}//if isDNA
		}//iterate through peppyList
		
		//Remove Duplicates
		// add elements to al, including duplicates
		HashSet<PeppyLine> hs = new HashSet<PeppyLine>();
		hs.addAll(masterDelete);
		masterDelete.clear();
		masterDelete.addAll(hs);
		
		

		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(this.outputDir + "/" + "FunnelRemovedLines.txt"));
			out.write("Total number of lines " + masterDelete.size() + "\n");
		//Remove peptides that have a lowered score
		for(PeppyLine td: masterDelete){
//			U.p(td.toString());
//			if(!td.isDna()){
				out.write(td.toString() + "\n");
//			}
			PeppyList.remove(td);
		}
		
		out.flush();
		out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		U.p("PeppyList "  + PeppyList.size());
	}//funnel
	/**
	 * getFDRResutls gets a peppy FDR results file and parses out the 0,5,10 percent FDR cutouff scores.
	 */
	private void getFDRResults(){
		
		try {
			Scanner s = new Scanner(new File(FDRFile));
			String token = s.nextLine();
			
			while(!token.contains("Precision")){
				token = s.nextLine();
				
			}
			
			//Skip to the zero% fdr
			s.next();
			s.next();
			token = s.next();
			
			//Grab the 0% FDR Value
//			U.p(token);
			zero_fdr = Double.valueOf(token);
			//Skip to the 5% FDR
			for(int i = 0; i < 15; i++){
				token = s.next();
			}
//			U.p(token);
			five_fdr = Double.valueOf(token);
			
			//Skip to the 10% FDR
			for(int i = 0; i < 15; i++){
				token = s.next();
			}
//			U.p(token);
			ten_fdr = Double.valueOf(token);
			
			
			//U.p("Values are " + zero_fdr + " " + five_fdr + " " + ten_fdr);
			
			
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}//Catch
		
	}//getFDRResults
	
	/**
	 * parsePeppyFile gets a files header information, and then parses out each line to the PeppyList
	 */
	private void parsePeppyFile(){
		try {
			Scanner s = new Scanner(new File(inputFile));
			String token;
			
			for(int i = 0; i < 4; i++){
				token = s.nextLine();
				header += token;
				header += "\n";
			}
			

			
			if(header.contains("ChromosomeName")){
				parseMod();
				isModified = true;
			}else{
				parseUnMod();
				isModified = false;
			}
			}catch (FileNotFoundException e){
				e.printStackTrace();
			}
			U.p("******Input Parse******");
			U.p("Total number of Peppy file lines parsed: " + (PeppyList.size() + matchRankFiltered.size()));
			U.p("Total number of match rank filtered lines: " + matchRankFiltered.size());
			U.p("Total number of valid peppy lines parsed: " + PeppyList.size());
	}//parsePeppyFile
	
	/**
	 * parseUnMOd takes in a unmodified Peppy results file, and adds each line to the PeppyList
	 */
	private void parseUnMod(){
		//Only read in lines up to the breakAt point.  This is useful for debugging.
		int breakAt = -1;
		int count = 0;
		
		try {
			Scanner s = new Scanner(new File(inputFile));
			String token = "";
			//Skip the header
			for(int i = 0; i < 4; i++){
				token = s.nextLine();
			}
			
			
			//Parse information for each line
			int index = 0;
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
			String sequenceFile = "";
			String sequenceDescription = "";
			int intronStart = -1;
			int intronStop = -1;
			int strand = -1;
			boolean isSpliced = false;
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

				
				
				//Get the sequence File
				token = s.next();
				sequenceFile = token;

				
				//Get the sequence description
				token = s.next();
				sequenceDescription = token;

				
				
				//Get the Intron-Start
				token = s.next();
				intronStart = Integer.parseInt(token);

				
				//Get the Intron-Stop
				token = s.next();
				intronStop = Integer.parseInt(token);

				
				//Get the Strand
				token = s.next();
				if(token.equals("+")){
					strand = Definitions.genomicStrandPOSITIVE;
				}else{
					strand = Definitions.genomicStrandNEGATIVE;
				}

				
				//Get the isSpliced
				token = s.next();
				isSpliced = Boolean.valueOf(token);

				
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
				
				
				if(matchRank <= maxMatchRank){
					//List of the peppy input file lines
					PeppyLine pl = new PeppyLine(ID, SpectrumID, MD5, FileName, score, PrecursorMZ, PrecursorNeutralMass, eValue, peptideSequence, start, stop, sequenceFile, sequenceDescription, intronStart, intronStop, strand, isSpliced, matchRank, RankCount, IonCount, labeled, charge,
					        cleavageAcidCount, inORF, hydrophobic, hydrophilic, isModified, modMass, modIndex);
					pl.setIndex(index);
					pl.setDna(true);
					PeppyList.add(pl);
					index++;
//					U.p(index);
				}else{
					matchRankFiltered.add(new PeppyLine(ID, SpectrumID, MD5, FileName, score, PrecursorMZ, PrecursorNeutralMass, eValue, peptideSequence, start, stop, sequenceFile, sequenceDescription, intronStart, intronStop, strand, isSpliced, matchRank, RankCount, IonCount, labeled, charge,
							        cleavageAcidCount, inORF, hydrophobic, hydrophilic, isModified, modMass, modIndex));
				}
				
//				U.p(PeppyList.get(PeppyList.size() - 1).toString());
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
				sequenceFile = "";
				sequenceDescription = "";
				intronStart = -1;
				intronStop = -1;
				strand = -1;
				isSpliced = false;
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
			}//While has next line
			
		//File is not found.
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		//Throw an error when a unparsable token is read in.
		} catch (NoSuchElementException e) {
			U.p("Invalid file.  Unable to completely parse input.  Output will not contain all peptides.");
		}//catch
		
	}//parseUnMod
	
	/**
	 * parseMod takes in a modified BedFileValidation tool file, and adds each line to the PeppyList
	 */
	private void parseMod(){
		//DEBUGGING
		//Only read in lines up to the breakAt point.  This is useful for debugging.
		int breakAt = -1;
		int count = 0;
		

		try {
				Scanner s = new Scanner(new File(inputFile));
				String token;
			
				for(int i = 0; i < 4; i++){
					token = s.nextLine();
		
				}

			
			
			//Parse information for each line
			int index = 0;
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
			int ChromosomeName = -1;
			int Strand = -1;
			int StartLocation = -1;
			int StopLocation = -1;	
			int BlockCount = -1;
			String BlockSize = "";
			String BlockStart = "";

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
				
				//Get the ChromosomeName
				token = s.next();

				/*Plug values from the fields array into variables*/
				if(token.equalsIgnoreCase("chrM")){
					ChromosomeName = Definitions.chromosomeM;
				}else if(token.equalsIgnoreCase("chrX")){
					ChromosomeName = Definitions.chromosomeX;
				}else if(token.equalsIgnoreCase("chrY")){
					ChromosomeName = Definitions.chromosomeY;
				}else{
					ChromosomeName = Integer.parseInt(token.substring(token.indexOf('r') + 1));
					
				}
				
				//Get the strand
				token = s.next();
				if(token.equalsIgnoreCase("+")){
					Strand = Definitions.genomicStrandPOSITIVE;
				}else {
					Strand = Definitions.genomicStrandNEGATIVE;
				}
				
				//Get the start Location
				token = s.next();
				StartLocation = Integer.parseInt(token);
				
				//Get the stop location
				token = s.next();
				StopLocation = Integer.parseInt(token);	
				
				//Get the block count
				token = s.next();
				BlockCount = Integer.parseInt(token);
				
				//Get the block size
				token = s.next();
				BlockSize = token;
				
				//get the blockStart
				token = s.next();
				BlockStart = token;

				s.nextLine();
				
				if(matchRank <= maxMatchRank){
				//List of the peppy input file lines
					PeppyLine pl = new PeppyLine(ID, SpectrumID, MD5, FileName, score, PrecursorMZ, PrecursorNeutralMass, eValue, peptideSequence, start, stop, proteinName, matchRank, RankCount, IonCount, labeled, charge,
					        cleavageAcidCount, inORF, hydrophobic, hydrophilic, isModified, modMass, modIndex, ChromosomeName, Strand, StartLocation, StopLocation, BlockCount, BlockSize, BlockStart );
					pl.setIndex(index);
					pl.setDna(false);
					PeppyList.add(pl);
					index++;
					
				}else{
					matchRankFiltered.add(new PeppyLine(ID, SpectrumID, MD5, FileName, score, PrecursorMZ, PrecursorNeutralMass, eValue, peptideSequence, start, stop, proteinName, matchRank, RankCount, IonCount, labeled, charge,
					        cleavageAcidCount, inORF, hydrophobic, hydrophilic, isModified, modMass, modIndex, ChromosomeName, Strand, StartLocation, StopLocation, BlockCount, BlockSize, BlockStart));
				}
				
//				U.p(PeppyList.get(PeppyList.size() - 1).toString());
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
				ChromosomeName = -1;
				Strand = -1;
				StartLocation = -1;
				StopLocation = -1;	
				BlockCount = -1;
				BlockSize = "";
				BlockStart = "";
				
				//DEBUGGING BREAK
				if(count == breakAt){
					U.p("Breaking early in the peptide list, with a break early value of " + breakAt);
					break;
				}
				
				count++;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		//Throw an error when a unparsable token is read in.
		} catch (NoSuchElementException e) {
			U.p("Invalid file.  Unable to completely parse input.  Output will not contain all peptides.");
		}//catch
	}//parseMod
	

	
	/**
	 * determinePeptideCount determines how many of each peptide occurs in the results file.  It does this by first grouping all peptides that occur in the same place, then amending 
	 * each peptide with information about how many of it occurs in the input file.
	 */
	private void determinePeptideCount(){
		HashMap<String, ArrayList<PeppyLine>> peptideHash = new HashMap<String, ArrayList<PeppyLine>>();
		
		
		for(PeppyLine pl: PeppyList){
			//Key is of the format chr#:start-stop:strand
			//Example is     chr1:100-203:+
			String key = "";
			
			//Key is the peptide
			key = pl.getPeptideSequence();
			if(pl.getPeptideSequence().endsWith(".")){
				key = pl.getPeptideSequence().substring(0, pl.getPeptideSequence().length() - 1);
			}

			
			//Attempt to pull this key out of the peptide hash
			ArrayList<PeppyLine> temp = peptideHash.get(key);
			
			//If key is null, the create a new entry with this key and peptide
			if(temp == null){
				temp = new ArrayList<PeppyLine>();
				temp.add(pl);
				peptideHash.put(key, temp);
				
				//Just add this peptide to the existing keys list of items
			}else{
				boolean pepFound = false;
				
				//check all of the items in this bucket, to ensure this peptide is not already in there
				for(int k = 0; k < temp.size(); k++){
					if(temp.get(k).getUniqueKey().equals(pl.getUniqueKey())){
						pepFound = true;
					}
				}
				
				if(!pepFound){
					//Determine if this peptide should be added.  See if it is the same as any other peptide.
					temp.add(pl);
				}
			}
			
			
		}//for loop
		
		//Create a collection and iterator to navigate through this hash map
		Collection<ArrayList<PeppyLine>> col = peptideHash.values();
		Iterator<ArrayList<PeppyLine>> iter = col.iterator();
		
		//Iterate through the peptides and tell each peptide how many copies of it there are in the genome
		ArrayList<PeppyLine> current = iter.next();
		
		//iterate through the peptides and amend each one to reflect the number of copies it has in the genome
		while(iter.hasNext()){
			
			//When a peptide is found with a unique location count, search for every other peptide like this and 
			
			//For each collection of peptides in the hash map, set all of them to reflect the appropriate size
			for(int i = 0; i < current.size(); i++){
				current.get(i).setLocCount(current.size());
				
				//for each peptide, search the list for every ohter peptide with this unqiue key and set the unique location key
				for(PeppyLine pl: PeppyList){
					//If the unique keys are the same, then set this pl's location count
					if(pl.getUniqueKey().equals(current.get(i).getUniqueKey())){
						pl.setLocCount(current.size());
					}
				}
			}
	
			//Increment the iterator
			current = iter.next();
		}
		

		
	}//determinePeptideCount
	
	
	/**
	 * removeDuplicates removes duplicate peptides from the bed file list.
	 */
	private void removeDuplicates(){
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(this.outputDir + "/" + "DuplicatesRemoved.txt"));

		
		ArrayList<Integer> markedForDeletion = new ArrayList<Integer>();
		
		int preSize = PeppyList.size();
		
		HashMap<String, PeppyLine> peptideHash = new HashMap<String, PeppyLine>();
		//Iterate through PeppyList
		for(int i = 0; i < PeppyList.size(); i++){
			PeppyLine p = PeppyList.get(i);
			if(p.getMD5().equals("593bca275f60b17d23e0356d454d380d")){// && !p.isDna()){
//				@SuppressWarnings("unused")
//				int y = 0;
			}
			PeppyLine exists = peptideHash.get(p.getUniqueKey() + p.getMD5());//p.getSpectrumID());
			
//			U.p("Key is " + p.getUniqueKey() + p.getSpectrumID());
			if(exists == null){
//				U.p("***");
				peptideHash.put(p.getUniqueKey() + p.getMD5(), p);
			}else{
				markedForDeletion.add(i);
//				p.delete = true;
//				PeppyList.remove(p);
//				if(p.isDna() && !exists.isDna() || !p.isDna() && exists.isDna()){
				if(!p.isDna()){
					out.write(p.toString() + "\n");
				}
				
//				U.p("Key is " + p.getUniqueKey() + p.getSpectrumID() + " Index " + i);
			}
			

		}//initial PeppyList
		for(int j = 0; j < markedForDeletion.size(); j++){
			PeppyList.remove(markedForDeletion.get(j) - j);
//			U.p("removing " + (markedForDeletion.get(j) - j) + " " + PeppyList.get(markedForDeletion.get(j) - j).getProteinName());
		}
	
		U.p("Total numnber of PeppyList lines before removing duplicates: " + preSize);
		U.p("Total number of Peppylist lines after removing duplicates: " + PeppyList.size());
		
		out.flush();
		out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}//removeDuplicates


	/**
	 * Calculate bed score calculates the bed score based on the following formula for each line the in the peppy line list:
	 * 1000 * ((Score - 10%FDR)/(0%FDR - 10%FDR))
	 * 
	 * This formula represents a peptides spectral match score relative to the 10% to 0% range.  A score of 0 means the peptide is below the 10% fdr, a score of 1000 means above the 
	 * 0% fdr, and any score between 0-1000 represents a score relative to the 10% to 0% range.
	 */
	private void calcBedScore(){
		
		Integer score = new Integer(0);
		//Iterate through each line of the PeppyList, and calculate its bed score
		for(PeppyLine pl: PeppyList){
			
			score = (int)(1000 * ((pl.getScore() - ten_fdr) / (zero_fdr - ten_fdr)));
			
			//Clamp the score to the range 0-1000
			if(score < 0){
				score = 0;
			}else if(score > 1000){
				score = 1000;
			}
			
			//Add this score to the 
			pl.setBedScore(score);

		}//For loop
	}//Calculate Bed Score
	/**
	 * createBedFiles creates two bed files, one filtered and one unfiltered.  Each bed file contains the standard tweleve columns + 4 additional informaiton columns.  In unmodified 
	 * input files those columns are: Peppy score , MD5, MatchRank, Location Count
	 * In modified input files those columns are: Peppy score , MD5, MatchRank, Modification Mass
	 * 
	 * Unique peptides (match rank = 1) are colored bright red, while all other peptides are colored black.
	 */
	private void createBedFiles(){
		//Writing to HDD related variables
		StringBuffer sb;
		BufferedWriter out;
		File outFile;
		
		U.p("******Output Totals******");
		
		//Keep track of the number of filtered and unfiltered bed file lines created
		int filteredCount = 0;
		int unFilteredCount = 0;
		String outputDir = "";
		String tab = "\t";
		String folder =  outputDirectoryFormat.format(cal.getTime());
		
		Collections.sort(PeppyList);
//		calcBedScore();
		
		try {
		for(int i = 0; i < 2; i++){
			/*Initialize all varibles required for creating output*/
			outFile = new File(this.outputDir + "/" + folder +  "/");
			outFile.mkdir();
			 outputDir = outFile.getAbsolutePath() + "/";
			
			
			sb = new StringBuffer();
			//This information is necessary for colors to be displayed in the UCSC browser
			sb.append("track itemRgb='On'" + "\n");
			out = null;
				for(int k = 0; k < PeppyList.size(); k++){
					PeppyLine pl = PeppyList.get(k);
					if(pl.getMD5().equals("593bca275f60b17d23e0356d454d380d")){// && !p.isDna()){
//							@SuppressWarnings("unused")
						U.p();
					}
					//Filter for the 5 percent fdr and location count < 10
					if(pl.getScore() < five_fdr){
						continue;
					}
					if(pl.getlocationCount() > 10){
						continue;
					}
					if(i == 0){
						//if 0 then only allow modified
						if(pl.isModified() == false){
							continue;
						}
						out = new BufferedWriter(new FileWriter(outputDir + "filtered_modified.txt" /*+ inputFile.substring(inputFile.lastIndexOf('/') + 1)*/));
					}else{
						//if 1 allow unmodified
						if(pl.isModified() == true){
							continue;
						}
						out = new BufferedWriter(new FileWriter(outputDir + "filtered_unmodified.txt" /*+ inputFile.substring(inputFile.lastIndexOf('/') + 1)*/));
					}
					
					
					String strand = "+";
					if(pl.getFirstStrand() == Definitions.genomicStrandNEGATIVE){
						strand = "-";
					}
					String color = "0,0,0";
					if(pl.getlocationCount() == 1){
						color = "255,0,0";
					}
					String peptideEnd = "";
					String lastPiece = "";
					if(pl.isModified()){
						lastPiece = String.valueOf(pl.getModMass());
						peptideEnd =  "-" + String.valueOf((int)pl.getModMass());
					}else{
						lastPiece = String.valueOf(pl.getlocationCount());
					}
					if(pl.isDna()){
						//IMP score (as text), MD5, MatchRank (should be all 1), Pep Count Genome
						sb.append(pl.getSequenceDescription() + tab + pl.getStart() + tab + pl.getStop() + tab + pl.getPeptideSequence() + peptideEnd + tab + pl.getBedScore() + tab + strand + tab +
								pl.getStart() + tab + pl.getStop() + tab + color + tab + "1" + tab + (pl.getStop() - pl.getStart()) + tab + "0" + tab + pl.getScore() + tab + pl.getMD5() + tab
								+ pl.getMatchRank() + tab + lastPiece +  "\n");
					//isModified = true
					}else{
						//IMP score (as text), MD5, MatchRank (should be all 1), modMass
						sb.append(Definitions.convertChrmNumToString(pl.getChromosomeName()) + tab + pl.getStartLocation() + tab + pl.getStopLocation() + tab + pl.getPeptideSequence() + peptideEnd + tab + pl.getBedScore() + tab +
								strand + tab + pl.getStartLocation() + tab + pl.getStopLocation() + tab + color + tab + pl.getBlockCount() + tab  + pl.getBlockSize() + tab + pl.getBlockStart() + tab +
								pl.getScore() + tab + pl.getMD5() + tab + pl.getMatchRank() + tab + lastPiece + "\n");
					}
					//chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
					//IMP score (as text), MD5, MatchRank (should be all 1), Pep Count Genome
					filteredCount++;
				}//for
			
			
			//Write the buffered information to disk, then flush and close the pipeline to ensure that this information is pushed to disk.
			out.write(sb.toString());
			out.flush();
			out.close();
			
			
			//Initialize all variables required for creating output
			out = null;
			sb = new StringBuffer();
			
			for(int k = 0; k < PeppyList.size(); k++){
				PeppyLine pl = PeppyList.get(k);
				
				//Filter for the 10 percent fdr and location count < 10
				if(pl.getScore() < ten_fdr){
					continue;
				}
//				if(pl.getlocationCount() > 10){
////					continue;
//				}
				if(i == 0){
					//if 0 then only allow modified
					if(pl.isModified() == false){
						continue;
					}
					out = new BufferedWriter(new FileWriter(outputDir + "unfiltered_modified.txt" /*+ inputFile.substring(inputFile.lastIndexOf('/') + 1)*/));
				}else{
					//if 1 allow unmodified
					if(pl.isModified() == true){
						continue;
					}
					out = new BufferedWriter(new FileWriter(outputDir + "unfiltered_unmodified.txt" /*+ inputFile.substring(inputFile.lastIndexOf('/') + 1)*/));
				}
				String strand = "+";
				if(pl.getFirstStrand() == Definitions.genomicStrandNEGATIVE){
					strand = "-";
				}
				String color = "0,0,0";
				if(pl.getlocationCount() == 1){
					color = "255,0,0";
				}
				
				String peptideEnd = "";
				String lastPiece = "";
				if(pl.isModified()){
					lastPiece = String.valueOf(pl.getModMass());
					peptideEnd = "-" + String.valueOf((int)pl.getModMass());
				}else{
					lastPiece = String.valueOf(pl.getlocationCount());
				}
				if(pl.isDna()){
					//Peppy score, MD5, MatchRank, Pep Count Genome
					sb.append(pl.getSequenceDescription() + tab + pl.getStart() + tab + pl.getStop() + tab + pl.getPeptideSequence() + peptideEnd + tab + pl.getBedScore() + tab + strand + tab +
							pl.getStart() + tab + pl.getStop() + tab + color + tab + "1" + tab + (pl.getStop() - pl.getStart()) + tab + "0" + tab + pl.getScore() + tab + pl.getMD5() + tab
							+ pl.getMatchRank() + tab + lastPiece + "\n");
				//isModified = true
				}else{
					//IMP score (as text), MD5, MatchRank, modMass
					sb.append(Definitions.convertChrmNumToString(pl.getChromosomeName()) + tab + pl.getStartLocation() + tab + pl.getStopLocation() + tab + pl.getPeptideSequence() + peptideEnd + tab + pl.getBedScore() + tab +
							strand + tab + pl.getStartLocation() + tab + pl.getStopLocation() + tab + color + tab + pl.getBlockCount() + tab  + pl.getBlockSize() + tab + pl.getBlockStart() + tab +
							pl.getScore() + tab + pl.getMD5() + tab + pl.getMatchRank() + tab + lastPiece +"\n");
				}
	
				unFilteredCount++;
			}//for
			
			//Write the buffered information to disk, then flush and close the pipeline to ensure that this information is pushed to disk.
			out.write(sb.toString());
			out.flush();
			out.close();
			if(i == 0){
				U.p("Modified");
			}else{
				U.p("Unmodified");
			}

			U.p("Total Number of filtered bed lines: " + filteredCount);
			U.p("Total Number of unfiltered bed lines: " + unFilteredCount);
			//*****STATS*****
			String slot = "";
			if(i == 0){
				slot = "modified";
			}else{
				slot = "unmodified";
			}
			out = new BufferedWriter(new FileWriter(outputDir + "statistics_" + slot + "_" + inputFile.substring(inputFile.lastIndexOf('/') + 1)));
			sb = new StringBuffer();
			
			SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
			cal = Calendar.getInstance();
			sb.append("This run was completed at: " + sdf.format(cal.getTime()) + "\n");
			sb.append("\n");
			sb.append("\n");
			sb.append("Input File: " + inputFile + "\n");
			sb.append("File isModified: " + isModified + "\n");
			sb.append("\n");
			sb.append("\n");
			sb.append("FDR Results: 0% " + zero_fdr + "\n");
			sb.append("             5% " + five_fdr + "\n");
			sb.append("            10% " + ten_fdr + "\n");
			sb.append("\n");
			sb.append("***Input Parse***" + "\n");
			sb.append( "Peppy file lines parsed: " + (PeppyList.size() + matchRankFiltered.size()) + "\n");
			sb.append("Match rank filtered lines: " + matchRankFiltered.size() + "\n");
			sb.append("Match rank filter value: " + maxMatchRank + "\n");
			sb.append("Valid peppy lines parsed: " + PeppyList.size() + "\n");
			sb.append("\n");
			sb.append("***Output***" + "\n");
			if(i == 0){
				sb.append("Modified" + "\n");
			}else{
				sb.append("Unmodified" + "\n");
			}
			sb.append("Filtered bed lines: " + filteredCount + "\n");
			sb.append("Unfiltered bed lines: " + unFilteredCount + "\n");
			//Write the buffered information to disk, then flush and close the pipeline to ensure that this information is pushed to disk.
			out.write(sb.toString());
			out.flush();
			out.close();
		}

		
		
		} catch (IOException e) {
			e.printStackTrace();
		}


	}//createBedFiles
}//GENCODESubmissionTool
