package MiscTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Scanner;

import MiscTools.GENCODEPPCompariosonTool.Translation;
import MiscTools.GENCODEPPCompariosonTool.TranslationMatch;
import PersonalProteome.Definitions;
import PersonalProteome.U;


/**
 * GencodePPComparisonTool is a proteome comparison tool designed to validate Personal Proteome.  It takes two FASTA files as input, one in the format of the published gencode annotations,
 * and the other is in the form of personal proteome output.  Each protein from the genocode translations is matched up with its corresponding Personal Proteome translation, and character by
 * character comparison is done.  Statistics about differences are stored in a separate output file, as well as information about matches/mismatches.
 * It currently outputs 4 files: Statistics, Matches with Variants, Unmatched Gencode translations, and Unmatched PP translations(Unfinished).
 * 
 * Expected Input lines: PP Three lines of >, with the third being transID|GENEID|, followed on the next lines by a protein
 * 						GENCODE One line of >, with the first two parts being transID|GENEID|, followed on the next line by a protein.
 * @author David "Corvette" Thomas
 *
 */
public class GencodePPComparisonTool {

	//Look at maximum memory usage 
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	
	public static void main(String[] args){
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();

		
		/*  hello! */
		printGreeting();
		
		String args0a = "/Users/davidthomas/Peppy/synthesis/PPGencodeComparisonTool/gencode.v11.pc_translations.fa";
		String args0b = "/Users/davidthomas/Peppy/synthesis/PPGencodeComparisonTool/PersonalProteomeOutput.fasta";
		String args1 = "/Users/davidthomas/Peppy/synthesis/PPGencodeComparisonTool/output/";
		String args2 = "results.txt";

		
		/*setup and do some work!*/
		GencodePPComparisonTool tool = new GencodePPComparisonTool(args0a, args0b, args1, args2);

		tool.Compare();
		
		/* i'm finished! */
		printFarewell();
		
		
	}
	
	public static void printGreeting() {		
		U.p("Proteome Comparison at the speed of a turbo-Vette!");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}

	//Time
	SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
	SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	Calendar cal;
	String startTime = "";
	String endTime = "";
	
	//File Location Information
	private String gencodeFile;
	private String ppFile; 
	private String outputDir; 
	private String outputFileName;
	
	//Storage of Translations
	ArrayList<Translation> ppTransList = new ArrayList<Translation>();
	ArrayList<Translation> gencodeTransList = new ArrayList<Translation>();
	Hashtable<String, TranslationMatch>  matchHash;
	
	/**
	 * GencodePPComparisonTool is a proteome comparison tool designed to validate Personal Proteome.  It takes two FASTA files as input, one in the format of the published gencode annotations,
	 * and the other is in the form of personal proteome output.  Each protein from the genocode translations is matched up with its corresponding Personal Proteome translation, and character by
	 * character comparison is done.  Statistics about differences are stored in a separate output file, as well as information about matches/mismatches.
	 * @param gencodeFile Gencode translations in FASTA format.
	 * @param ppFile  PP output file.
	 * @param outputDir  Directory to store output information.
	 * @param outputFileName Name to assign to output fields.
	 */
	public GencodePPComparisonTool(String gencodeFile, String ppFile, String outputDir, String outputFileName){
		this.gencodeFile = gencodeFile;
		this.ppFile = ppFile;
		this.outputDir = outputDir;
		this.outputFileName = outputFileName;
	}
	
	/**
	 * Compare is the primary method of GENCODEPPComparison.  Once called, this method handles uploading translations, comparing them, and creating output.
	 */
	public void Compare(){
		//Get the current time
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		
		U.startStopwatch();
		U.p("Loading Personal Proteome file.");
		loadPP();
		U.stopStopwatch();
		
		U.startStopwatch();
		U.p("Loading Gencode file.");
		loadGencode();
		U.stopStopwatch();
		
		U.startStopwatch();
		U.p("Creating GENCODE based hash table.");
		createGencodeHash();
		U.stopStopwatch();
		
		U.startStopwatch();
		U.p("Adding PP to the hash table.");
		addPPToGencode();
		U.stopStopwatch();
		
		
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.startStopwatch();
		U.p("Comparing PP and Gencode, and creating output files");
		compareAndCreateOutput();
		U.stopStopwatch();
		
		
		//Get the ending time

		U.p("Finishing up: " + endTime);
	}
	
	/**
	 * loadPP uploads and stores the information in the constructor specified Personal Proteome output.
	 */
	private void loadPP(){
		File ppFile = new File(this.ppFile);
		StringBuffer sequence = new StringBuffer();
		String transID = "";
		String geneID = "";

		try {
			Scanner s = new Scanner(ppFile);
			String token;
			int arrowCount = 0;
			while(s.hasNextLine()){
				//Ensure that a valid line is coming up next.  THIS IS IMPORTANT****  OTHERWISE IF THE LAST PROTEIN CONTAINS A EVEN NUMBER OF LINES THIS METHOD WILL THROW AN EXCEPTION.
				if(s.hasNext()){
					token = s.next();
				}else{
					break;
				}

				//Skip header information that is not important
				if(token.equals(">")){
					if(arrowCount < 2){
						arrowCount++;
						s.nextLine();
						continue;
					}
					arrowCount = 0;
				}
	
				token = s.nextLine();
				token = token.trim();
				
				transID = token.substring(0, token.indexOf("|"));
				geneID = token.substring(transID.length() + 1, token.length() - 1);

				//Loop through and pull out the protein from the fasta file.
				token = s.nextLine();
				while(token.length() > 0){
					sequence.append(token);
					token = s.nextLine();
				}		
				
				ppTransList.add(new Translation(geneID, transID, sequence.toString()));
				 	
				sequence = new StringBuffer();
				transID = "";
				geneID = "";
			}
			
			U.p("Total number of P.P. Translations read in: " + ppTransList.size());
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		
	}
	
	/**
	 * loadGencode uploads and stores the proteins from the constructor specified GENCODE translations.
	 */
	private void loadGencode(){
		File gencodeFile = new File(this.gencodeFile);
		StringBuffer sequence = new StringBuffer();
		String transID = "";
		String geneID = "";
		try {
			Scanner s = new Scanner(gencodeFile);
			String token;
			while(s.hasNextLine()){
				token = s.nextLine();
				transID = token.substring(1, token.indexOf("|"));
				geneID = token.substring(transID.length() + 1 + 1, token.indexOf("|", 1 + transID.length() + 1));
				
				token = s.nextLine();
				sequence.append(token);
				
				gencodeTransList.add(new Translation(geneID, transID, sequence.toString()));
				
				
				sequence = new StringBuffer();
				transID = "";
				geneID = "";
			}
			
			U.p("Total number of gencode Translations read in: " + gencodeTransList.size());
			
		}catch (FileNotFoundException e){
			e.printStackTrace();
		}
	}
	
	/**
	 * createGencodeHash stores the uploaded GENCODE translations in a hashtable for optimization purposes.  It allows for constant seek time for a GENCODE translation, so when going through
	 * the PP translations (size N) a total time for comparisons is O(N).
	 */
	private void createGencodeHash(){
		//Values are to set a high enough inital capacity and load factor as to never need a resize
		matchHash = new Hashtable<String, TranslationMatch>(gencodeTransList.size() + 100, 1);
		TranslationMatch matchTrans;
		
		for(int i = 0; i < gencodeTransList.size(); i++){
			Translation gencodeT = gencodeTransList.get(i);
			matchTrans = new TranslationMatch();
			matchTrans.setGencodeTrans(gencodeT);
			matchHash.put((gencodeT.getTranscriptID() + gencodeT.getGeneID()), matchTrans);
		}
		
		U.p("Number of successfully hashed gencode trasnlations: " + matchHash.size());
	}
	
	/**
	 * addPPToGencode matches up each GENCODE translation with its accompanying PP translation.
	 */
	private void addPPToGencode(){
		
		
		for(int i = 0; i < ppTransList.size(); i++){
			Translation ppT = ppTransList.get(i);
			TranslationMatch tm = matchHash.get(ppT.getTranscriptID() + ppT.getGeneID());
			
			if(tm != null){
				//Mark these Translations as matched
				ppT.setMatchType(Definitions.MATCH_TYPE_MATCHED);
				tm.getGencodeTrans().setMatchType(Definitions.MATCH_TYPE_MATCHED);
	
				//Add the personal proteome translation to its corresponding GENCODE translation
				tm.setPpTrans(ppT);
			}
			
			
			
		}//for
		
		int count = 0;
		for(Translation t: ppTransList){
			if(t.getMatchType() == 1){
				count++;
			}
		}
		
		
		U.p("Number of successfully matched P.P. translations: " + count);
		U.p("Hash map size " + matchHash.size());
	}
	
	/**
	 * compareAndCreateOutput compares all matches and calculates statistics about this information.  It then generates about files based on the data.
	 */
	private void compareAndCreateOutput(){
		Collection<TranslationMatch> matchesCol = matchHash.values();
		Iterator<TranslationMatch> iterator;
		StringBuffer proteinsToWriteOut = new StringBuffer();
		StringBuffer unMatchedGencode = new StringBuffer();
		iterator = matchesCol.iterator();
		//Don't do anything if there are 
		int lineLength = 80;
		
		
		if(matchesCol.size() == 0){
			U.p("No Matches");
			return;
		}
		
		
		TranslationMatch token;
		
		int totalMatchCount = 0;
		int ppNotMatchedCount = 0;
		int gencodeNotMatchedCount = 0;
		
		
		//This goes through the Collection 1 match at a time.  Stats should be collected here.
		while(iterator.hasNext()){
			token = iterator.next();
			//Determine the total number of matches
			if(token.getGencodeTrans().getMatchType() == Definitions.MATCH_TYPE_MATCHED){
				totalMatchCount++;
			}
			
			//Determine the number of unmatched translations
			if(token.getGencodeTrans().getMatchType() == Definitions.MATCH_TYPE_UNMATCHED){
				gencodeNotMatchedCount++;
				if(gencodeNotMatchedCount < 10){
					U.p("Unmatched gencode protein: " + token.getGencodeTrans().toString());
				}
			}

			
		}//while
		
		for(Translation t: ppTransList){
			if(t.getMatchType() == Definitions.MATCH_TYPE_UNMATCHED){
				ppNotMatchedCount++;
			}
		}//for loop
		

		iterator = matchesCol.iterator();

		int misMatch = 0;
		int misMatchLength = 0;
		int hasVariance = 0;
		
		int ppLonger = 0;
		int gencodeLonger = 0;
		
		ArrayList<Integer> variants = new ArrayList<Integer>();
		ArrayList<Double> variantAsPercentageOfLength = new ArrayList<Double>();
		//Mark in lower case the variances that PP has from the GENCODE annotation.
		while(iterator.hasNext()){
			token = iterator.next();
			
			Translation genTran = token.getGencodeTrans();
			Translation ppTran = token.getPpTrans();
			//Ignore unamtched translations.  Handle them in a separate file
			if(ppTran == null){
				
				//Write the GENCODE to the gencode file
				int linesOfProtein = genTran.getSequence().length() / lineLength;
				unMatchedGencode.append(">GENCODE|" + genTran.getTranscriptID() + "|" + genTran.getGeneID() + "|" + "VariantCount: " + genTran.getVariants() + "|" + "Length: " + genTran.getSequence().length() + "|" + "\n");
				
				for(int k = 0; k < linesOfProtein + 1; k++){
					//If it is the last line, just write all of them
					if(k == linesOfProtein){
							unMatchedGencode.append(genTran.getSequence().substring(k*lineLength));
							unMatchedGencode.append("\n");
						
					}else{
					//Write 80 characters
						unMatchedGencode.append(genTran.getSequence().substring(k*lineLength, (k+1)*lineLength));
						unMatchedGencode.append("\n");
					}
					
				}//for

				unMatchedGencode.append("\n");

				continue;
			}
			int shortestTranLength = 0;
			
			if(genTran.getSequence().length() > ppTran.getSequence().length()){
				shortestTranLength = ppTran.getSequence().length();
			}else{
				shortestTranLength = genTran.getSequence().length();
			}//else
			
			//Convert any characters that are different to lower case.
			for(int n = 0; n < shortestTranLength; n++){
				//Check to determine if a substitution has to occur
				if(ppTran.getSequence().charAt(n) != genTran.getSequence().charAt(n)){
					
					//If the comparison is a delayed start ignore it
					if(ppTran.getSequence().charAt(n) == Definitions.START_NOT_FOUND_CHAR){
						continue;
					}
					ppTran.incrementVariants();
					String charToLower = ppTran.getSequence().substring(n, n + 1);
					charToLower = charToLower.toLowerCase();
					String start = ppTran.getSequence().substring(0, n);
					String end = ppTran.getSequence().substring(n + 1);
					ppTran.setSequence(start + charToLower + end);
					
				}//if
				
			}//for
			
			//Ignore proteins where the variantCount and length are the same
			if(ppTran.getSequence().length() == genTran.getSequence().length() && ppTran.getVariants() == genTran.getVariants()){
				continue;
			}
			

			if(!(genTran.getSequence().length() - ppTran.getSequence().length() == 1 &&  ppTran.getVariants() == genTran.getVariants())){
//				continue;
			}
			
			if(ppTran.getVariants() == 0){
//				continue;
			}
			
			
			//pp longer
			if(genTran.getSequence().length() < ppTran.getSequence().length()){
				ppLonger++;
			}else if(genTran.getSequence().length() > ppTran.getSequence().length()){
				gencodeLonger++;
			}
			//gencode Longer
			//Keep track of how many matched translations had a variance
			misMatch++;
			if(ppTran.getSequence().length() != genTran.getSequence().length()){
				misMatchLength++;
			}
			if(ppTran.getVariants() != genTran.getVariants()){
				hasVariance++;
				variants.add(ppTran.getVariants());
				variantAsPercentageOfLength.add(((double)ppTran.getVariants())/ppTran.getSequence().length());
			}
			
			int linesOfProtein = ppTran.getSequence().length() / lineLength;
			proteinsToWriteOut.append(">PersonalProteome|" + ppTran.getTranscriptID() + "|" + ppTran.getGeneID() + "|" + "VariantCount: " + ppTran.getVariants() + "|" + "Length: " + ppTran.getSequence().length() + "|" + "\n");
			
			for(int k = 0; k < linesOfProtein + 1; k++){
				//If it is the last line, just write all of them
				if(k == linesOfProtein){
						proteinsToWriteOut.append(ppTran.getSequence().substring(k*lineLength));
						proteinsToWriteOut.append("\n");
					
				}else{
				//Write 80 characters
					proteinsToWriteOut.append(ppTran.getSequence().substring(k*lineLength, (k+1)*lineLength));
					proteinsToWriteOut.append("\n");
				}
				
			}//for

			proteinsToWriteOut.append("\n");
			linesOfProtein = genTran.getSequence().length() / lineLength;
			
			proteinsToWriteOut.append(">GENCODE|" + genTran.getTranscriptID() + "|" + genTran.getGeneID() + "|" + "VariantCount: " + genTran.getVariants() + "|"+ "Length: " + genTran.getSequence().length() + "|" + "\n");
			
			for(int k = 0; k < linesOfProtein + 1; k++){
				//If it is the last line, just write all of them
				if(k == linesOfProtein){
						proteinsToWriteOut.append(genTran.getSequence().substring(k*lineLength));
						proteinsToWriteOut.append("\n");
					
				}else{
				//Write 80 characters
					proteinsToWriteOut.append(genTran.getSequence().substring(k*lineLength, (k+1)*lineLength));
					proteinsToWriteOut.append("\n");
				}
				
			}//for
			
			                                                                
			proteinsToWriteOut.append("\n");
			
			
			
			
		}//iterator
		
		
		try {
			
			File outFile = new File(this.outputDir + outputDirectoryFormat.format(cal.getTime()) +  "/");
			outFile.mkdir();
			String outputDir = outFile.getAbsolutePath() + "/";
			
			BufferedWriter gencodeUNMATCHED	 = new BufferedWriter(new FileWriter(outputDir + "gencodeUNMATCHED"+ outputFileName));;
			BufferedWriter ppUnmatched = new BufferedWriter(new FileWriter(outputDir + "ppUNMATCHED"+ outputFileName));
			BufferedWriter proteinOutput= new BufferedWriter(new FileWriter(outputDir + "MATCHEDbutVariant"+ outputFileName));
			proteinOutput.append(proteinsToWriteOut.toString());
			
			

			
			if(ppNotMatchedCount > 0){
		
				U.p("CODE FOR pp unmatched proteins is not written!*************************");
				
				
			}
			
			if(gencodeNotMatchedCount > 0){
				gencodeUNMATCHED.append(unMatchedGencode.toString());
			}
			
			BufferedWriter statsOut= new BufferedWriter(new FileWriter(outputDir + "STATS" + outputFileName));
			statsOut.write("This run was started at " + startTime + "\n");
			statsOut.write("This run was completed at " + endTime + "\n");
			statsOut.newLine();
			statsOut.write("The gencode proteome file is " + gencodeFile + "\n");
			statsOut.write("The Personal Proteome file is " + ppFile + "\n");
			statsOut.newLine();
			statsOut.write("********STATS**********" + "\n");
			statsOut.write("Total number of GENCODE Translations is: " + gencodeTransList.size() + "\n");
			statsOut.write("Total number of P.P. Translations read in: " + ppTransList.size() + "\n");
			statsOut.newLine();
			statsOut.write("The total number of matched translations is: " + totalMatchCount + "\n");
			statsOut.write("The total number of Personal Proteome unmatched translations is: " + ppNotMatchedCount + "\n");
			statsOut.write("The total number of GENCODE unmatched translations is:           " + gencodeNotMatchedCount + "\n");
			
			statsOut.newLine();
			statsOut.write("The total number of translations that were matched based on id, but have a variance or differing length in their sequences is: " + misMatch + "\n");
			statsOut.write("The total number of translations that match completely is: " + (totalMatchCount - misMatch) + "\n");
			
			statsOut.newLine();
			statsOut.write("The total number of translations that are mismatched based on length is: " + misMatchLength + "\n");
			statsOut.write("The total number of translations that are mismatched because of variance in sequence is: " + hasVariance + "\n");
			int onlyLength = misMatch - hasVariance;
			int onlyVariant = misMatch - misMatchLength;
			int hasBoth = misMatch - onlyLength - onlyVariant;
			statsOut.newLine();
			statsOut.write("The number of mismatched translations that are mismatched only on length is: " + onlyLength +  "\n");
			statsOut.write("GENCODE had the longest Length: "  + gencodeLonger + " \n");
			statsOut.write("PP had the longest Length: "  + ppLonger + " \n");
			statsOut.write("The number of mismatched translations that are mismatched only on variance is: " + onlyVariant + "\n");
			statsOut.write("The number of mismatched translations that are mismatched on both length and variance is: " + hasBoth + "\n");
			statsOut.newLine();
			Collections.sort(variants);
			int totalVariantValue = 0;
			int maxVariant = 0;
			for(Integer inte: variants){
				totalVariantValue += inte;
				if(inte > maxVariant){
					maxVariant = inte;
				}
			}
			if(variants.size() > 0){
			statsOut.write("The median value of the number of variants in a translation with at least 1 variant is: " + variants.get(variants.size() / 2) + "\n");
			
			statsOut.write("The average value of the number of variants in a translation with at least 1 variant is: " + (((double)totalVariantValue)/ variants.size()) + "\n");
			statsOut.write("The maximum number of variants a single translation has is: " + maxVariant + "\n");
			}else{
				statsOut.write("There are 0 varaints between GENCODE and PP.");
			}
			statsOut.newLine();
			Collections.sort(variantAsPercentageOfLength);
			double maxPercent = 0;
			double percentTotal = 0;
			for(Double d: variantAsPercentageOfLength){
				percentTotal += d;
				if(d > maxPercent){
					maxPercent = d;
				}
			}
			if(variantAsPercentageOfLength.size() > 0){
			statsOut.write("The median value of the percent of variants as compared to length (#variants/length) is: " + variantAsPercentageOfLength.get(variantAsPercentageOfLength.size() / 2)*100 + "%" + "\n");
			statsOut.write("The average value of the percent of variants as compared to length (#variants/length) is: " + (percentTotal/variantAsPercentageOfLength.size())*100 + "%" + "\n");
			statsOut.write("The maximum value of the percentage of variants as compared to length (#variants/length) is: " + (maxPercent*100) + "%" + "\n");
			}
			statsOut.write("********STATS**********" + "\n");
			
			if(gencodeNotMatchedCount > 0){
				gencodeUNMATCHED.flush();
				gencodeUNMATCHED.close();
				
			}
			
			
			proteinOutput.flush();
			proteinOutput.close();
			statsOut.flush();
			statsOut.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
		
		
	}//compareAndCreateOutput
	
}//class
