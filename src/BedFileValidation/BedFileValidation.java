package BedFileValidation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Scanner;

import BedFileValidation.Objects.BedFileLine;
import Peppy.Sequence_DNA;

import PersonalProteome.Definitions;
import PersonalProteome.U;

/**
 * BedFileValidation is a class that takes in a bed file, and then validates that the coordinates in that bed file create the protein specified.
 * @author David "Corvette" Thomas
 *
 */
public class BedFileValidation {

	
	//Input files
	private String bedFile;
	private String outputDir;
	private String[] chrmFile = new String[25];
	
	private boolean genomeHasMito = false;
	
	//Date/Time information
	//Information for writing out files

	private SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	private Calendar cal;

	
	
	//Storage of Input
	ArrayList<BedFileLine> BEDlist = new ArrayList<BedFileLine>();
	ArrayList<String> misMatchedList = new ArrayList<String>();
	ArrayList<String> matchedList = new ArrayList<String>();
	
	
	/**
	 * BedFileValidation is a class that takes in a bed file, and then validates that the coordinates in that bed file create the protein specified.
	 * @param bedFile
	 * @param chrmDir
	 * @param outputDir
	 */
	public BedFileValidation(String bedFile, String chrmDir, String outputDir, Calendar currentTime) {
		//Establish default values
		this.bedFile = bedFile;
		this.outputDir = outputDir;
		
		//Set the specific time for creating a output file.
		cal = currentTime;
		
		//Determine if the mitochondrion DNA is to be included
		File mitoFile = new File(chrmDir + "chrM.fa");
		genomeHasMito = mitoFile.exists();
		
		//Get the file locaitons for each chromosome file
		populateChrmArray(chrmDir);
	}
	
	
	/**
	 * validate is the main driving method in BedFileValidation.  It uploads the bed file, creates and compares the proteins from the bed file, and calculates statistic and create output.
	 */
	public void validate(){

		//Upload Bed File
		U.p("Uploading bed file.");
		U.startStopwatch();
		uploadBedFile();
		U.stopStopwatch();
		
		U.p("Number of bed file lines uploaded: " + BEDlist.size());
		
		//Create peptide List
		U.p("Creating peptide list and comparing.");
		U.startStopwatch();
		createAndCompare();
		U.stopStopwatch();
		
		
		//calcStatsAndOutput
		U.p("Calculating stats and creating output file.");
		U.startStopwatch();
		calcStatsAndCreateOutput();
		U.stopStopwatch();
		
		
		

	}
	
	/**
	 * uploadBedFile allows uploads and stores lines for each bed file.
	 */
	private void uploadBedFile(){
		File peptideFile = new File(bedFile);
		//
		int id = 0;
		int chromosomeName = -1;
		int startLocation = -1;
		int stopLocation = -1;
		String sequence = "";
		String score = "";
		int strand = -1;
		String color = "";
		String blockCount = "";
		String blockSize = "";
		String blockStart = "";
		
		
		//Parse each line of the file and save it
		try{
			Scanner s = new Scanner(peptideFile);
			String token;
			
			//Parse each Line, one at a time
			while(s.hasNextLine()){
				
				if(!s.hasNext()){
					break;
				}
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
				
				//Skip the repeated start and stop location
				token = s.next();
				token = s.next();
				
				token = s.next();
				color = token;
				
				
				token = s.next();
				blockCount = token;
				
				token = s.next();
				blockSize = token;
				
				token = s.next();
				blockStart = token;
				
				//Add the newly read line into the peptideLocationList
				BEDlist.add(new BedFileLine(id ,chromosomeName, startLocation, stopLocation, sequence, score, strand, color, blockCount, blockSize, blockStart));
				
				id++;
				 //Reset all variables
				chromosomeName = -1;
				startLocation = -1;
				stopLocation = -1;
				sequence = "";
				score = "";
				strand = -1;
				color = "";
				blockCount = "";
				blockSize = "";
				blockStart = "";
				

			}//while

		//Let the user know if the bed file is not found	
		}catch(FileNotFoundException e){
			U.p("Error populating Bed File List: " + e);

		}//catch
	}//upLoadBedFile
	
	
	/**
	 * createAndCompare iterates through each chromosome file and grabs each bed files DNA sequence.  
	 */
	private void createAndCompare(){
		
		
		//Grab each chromosome file.
		Sequence_DNA seq;
		StringBuffer sb;
		//Iterate through each chromosome file, and create all of the transcripts possible from each one
		//Code is performed in this order to minize HDD reads.
		for(int k = 0; k < 25; k++){
			//Ignore the mitochondrian DNA
			if(k == 22){
				if(!genomeHasMito){
					continue;
				}
			}
						
			seq = new Sequence_DNA(chrmFile[k]);
			sb = new StringBuffer(seq.getNucleotideSequences().get(0).getSequence());
			BedFileLine beingTested;
			
		
			//Check each transcript to determine if it is using the currently open file.
			for(int j = 0; j < BEDlist.size(); j++){
				//Current BEDlist line being worked on.
				beingTested = BEDlist.get(j);
					
				//Ignore lines that do not occur on the chromosome currently being worked on. chromosome
				if(beingTested.getChromosomeName() != k + 1){
					continue;
				}

				int blockCount = Integer.valueOf(beingTested.getBlockCount());
				String dnaSeq = sb.substring(beingTested.getStartLocation(), beingTested.getStopLocation());;
				
				//Handle the situation where there are multiple blocks to string together.
				if(blockCount > 1){
					
					String buildingSeq = "";
	
					int blckStart;
					int blckSize;
					
					Scanner starts = new Scanner(beingTested.getBlockStart());
					Scanner sizes = new Scanner(beingTested.getBlockSize());
					
					starts.useDelimiter(",");
					sizes.useDelimiter(",");
					
					for(int p = 0; p < blockCount; p++){
						blckStart = (Integer.valueOf(starts.next()));
						blckSize = (Integer.valueOf(sizes.next()));
						if(beingTested.getStrand() == Definitions.genomicStrandPOSITIVE){

						}
							
						
						buildingSeq += dnaSeq.substring(blckStart, blckStart + blckSize);
					}
					
					
					dnaSeq = buildingSeq;
				}
				
					beingTested.setDnaSeq(dnaSeq);
			}//inner for
			

		}//outer for
		

		
	}//Create and compare
	
	
	/**
	 * calcStatsAndCreateOutput translates each transcript determined by create and compare.  It then checks these proteins for variances and calculates statistics based on that information.
	 * It creates output files based on the name of the input, by adding a RESULTS tag to the front of the name.
	 */
	private void calcStatsAndCreateOutput(){
		
		int variantCountTotal = 0;
		int dnaSeqLengthMismatch = 0;
		int aminoLengthMistacth = 0;
		//Convert and compare each peptide
		/*Brian Risk*/
		for(int i = 0; i < BEDlist.size(); i++){	
		
			
			
			//Create a protein for each line in the bed file
			String protein;
			String dnaSequence = BEDlist.get(i).getDnaSeq();
			int strand = BEDlist.get(i).getStrand();
			
			if(dnaSequence.equals("")){
				continue;
			}
			
			//Convert the DNA into a Protein
			char [] codon = new char[3];
			char aminoAcid;
			int mod = 0;
			
			int increment = 1;
			int startPosition = 0;
	
			StringBuffer buildingProtein = new StringBuffer();
			
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
	
					
					//Adjust for mitochondrian dna
					if(BEDlist.get(i).getChromosomeName() == Definitions.chromosomeM){
						aminoAcid = Definitions.mitoAminoAcidList[PersonalProteome.Gene.Transcript.indexForCodonArray(codon, isForwardsStrand)];
	
					}else{
						aminoAcid = Definitions.aminoAcidList[PersonalProteome.Gene.Transcript.indexForCodonArray(codon, isForwardsStrand)];
					}
	
					buildingProtein.append(aminoAcid);
					
					
					/* reset mod */
					mod = 0;
				} else {
					mod++;
				}
			}
			
			protein = buildingProtein.toString();
			
			//Compare the created protein to the sequence already in the bed file.
			
			
			int variantCount = 0;
			//Determine which acids are different from the reference genome.
			int limit = 0;
			boolean lengthMismatch = false;
			String refProtein = BEDlist.get(i).getSequence();
			if(protein.length() > refProtein.length()){
				limit = refProtein.length();
				lengthMismatch = true;
			}else if(refProtein.length() > protein.length()){
				limit = protein.length();
				lengthMismatch = true;
			}else{
				//the proteins are the same length
				limit = protein.length();
			}
			
			//Convert any characters that are different to lower case.
			for(int n = 0; n < limit; n++){
				//Check to determine if a substitution has to occur
				if(protein.charAt(n) != refProtein.charAt(n)){
					
					//If the comparison is a delayed start ignore it
					if(protein.charAt(n) == Definitions.START_NOT_FOUND_CHAR){
						continue;
					}
					variantCount++;
					String charToLower = protein.substring(n, n + 1);
					charToLower = charToLower.toLowerCase();
					String start = protein.substring(0, n);
					String end = protein.substring(n + 1);
					protein = start + charToLower + end;
					
				}
				
			}
			
			//Determine which errors happen with the bed file line.
			if(variantCount > 0 || lengthMismatch || BEDlist.get(i).getDnaSeq().length() != protein.length()*3){
				if(variantCount > 0){
					variantCountTotal++;
				}
				if(lengthMismatch){
					aminoLengthMistacth++;
					
				}
				if(BEDlist.get(i).getDnaSeq().length() != protein.length()*3){
					dnaSeqLengthMismatch++;
				}
				
				misMatchedList.add(variantCount + "\t" + protein + "\t" + BEDlist.get(i));
			}else{
				matchedList.add(variantCount + "\t" + protein + "\t" + BEDlist.get(i));
			}
		}//iterate bed file
		
		
		//Writing to HDD related variables
		StringBuffer sb;
		BufferedWriter out;
		File outFile;
		
		
		
		outFile = new File( outputDir + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");
		outFile.mkdirs();
		String outputDir = outFile.getAbsolutePath() + "/";

		try {
			sb = new StringBuffer();
			out = new BufferedWriter(new FileWriter(outputDir + "RESULT_STATS_" + bedFile.substring(bedFile.lastIndexOf("/") + 1)));
			
			
			sb.append("The genome directory used by this run is: " + chrmFile[0].substring(0, chrmFile[0].lastIndexOf('/' ) + 1));
			
			sb.append("\n");
			
			sb.append("Total number of bedfile lines read in: " + BEDlist.size() + "\n");
			sb.append("Lines with a matched peptide: " + matchedList.size() + "\n");
			sb.append("Lines without a matched peptide: " + misMatchedList.size()+ "\n");
			sb.append("Total number of lines with a variant amino acid: " + variantCountTotal + "\n");
			sb.append("Total number of lines where dnaSeqLength*3 != amino acid length: " + dnaSeqLengthMismatch + "\n");
			sb.append("Total number of lines with a where peptide length was mismatched: " + aminoLengthMistacth + "\n");
			
			sb.append("\n\n\n");
			sb.append("Results from bed file that did not match." + "\n");
			sb.append("VaraintCount" + "\t" +"ProteinProduced" + "\t" + "BedLine" + "\n");
			for(String s: misMatchedList){
				sb.append(s + "\n");
			}
			
			sb.append("\n\n\n");
			sb.append("Results from bed file that matched." + "\n");
			sb.append("VaraintCount" + "\t" +"ProteinProduced" + "\t" + "BedLine" + "\n");
			for(String s: matchedList){
				sb.append(s + "\n");
			}
			
			//Write the file and close it out.
			out.write(sb.toString());
			out.flush();
			out.close();
			
		//Catch a file expception, specifically with file issues.
		} catch (IOException e) {
			e.printStackTrace();
		}
	}//calcStatsandCreateOutput
	
	

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
	}//populateChrmArray
	
}//BedFileValidation
