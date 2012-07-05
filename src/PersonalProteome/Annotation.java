package PersonalProteome;
//import Peppy.Sequence_DNA;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Calendar;
import java.text.SimpleDateFormat;






import Peppy.Sequence_DNA;
import PersonalProteome.Definitions;
import PersonalProteome.U;
import PersonalProteome.Gene.CDS;
import PersonalProteome.Gene.DataPoint;
import PersonalProteome.Gene.Mutation;
import PersonalProteome.Gene.RegionOfInterest;
import PersonalProteome.Gene.Seleno;

import PersonalProteome.Gene.Transcript;

/**
 *Annotation loads in a location for an annotation and a directory of chromosomes, and a regions of interest file.  It then uploads the annotation file and creates a list of transcripts
 * to produce proteins from.  These transcirpts are stored, compared, put into an output file and have data/statistics run on them.  The synthesis method has several parameters that
 * can be set to affect how to  and performs protein synthesis.
 * @author David "Corvette" Thomas
 *
 */
public class Annotation{
	
	//Default behavior variables
	private boolean modStopStart = Properties.useModifiedStopsAndStarts;
	private boolean displayInfo = false;
	private boolean preFilterStarts = Properties.preFilterOutTranWithoutStartCodon;
	
	//The length of the line in the output
	int lineLength = 80;
	
	//File Locations
	//File location of the regions of interest
	private String regionsFile;
	//File Location of Gencode Annotation to read from.
	private String annotationFile;
	//Dir to store output files
	private String proteinFASTAOUTPUTDIR;
	//Array of Chromosome file locations
	private String refchrmDir;
	private String chrmDir;
	private String[] chrmFile = new String[25];
	
	//Variable to store whether or not this annotation needs to take mitochondrian DNA into account.
	private boolean genomeHasMito = false;
	
	//Name for the output Files
	private String outputFileName;
	private String statsOutputFile;
	private String dataPointsFileName;
	
	//ArrayList of lines from the GTF File
	//This Arraylist contains all of the information from the annotation.
	private ArrayList<GENCODE_GTF_Line> GTFlist;
	
	

	//Special regions of interest that have a score associated with their likelihood of mutations
	private ArrayList<RegionOfInterest> regionsOfInterest = new ArrayList<RegionOfInterest>();
	
	//A list of transcripts to produce proteins
	private ArrayList<Transcript> transcriptList;

	//A reference list to compare the list being tested to.  Typically HG19
	private ArrayList<Transcript> referenceTranscriptList;
	
	//Date/Time information
	//Information for writing out files
	private SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
	private SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	private Calendar cal;
	private Calendar endCal;
	private String startTime = "";
	private String endTime = "";

	/**
	 * Use this when the populateGTFList method is desired to be used.  It will allow a GTF file to be read in using the Personal Proteome method.
	 */
	public Annotation(boolean hasMito){
		genomeHasMito = hasMito;
	}
	
	
	/**
	 * For use with Proteome Lite, this simply allows for a single proteome to be created.
	 * @param annotationFile
	 * @param genomeDirectory
	 * @param chrmDirectory
	 * @param outputDir
	 */
	public Annotation(String annotationFile, String chrmDirectory, String outputDir){
		this.proteinFASTAOUTPUTDIR = outputDir;
		this.annotationFile = annotationFile;
		this.chrmDir = chrmDirectory;
		
		File mitoFile = new File(chrmDirectory + "chrM.fa");
		genomeHasMito = mitoFile.exists();
	}
	/**
	 * Annotation loads in a location for an annotation and a directory of chromosomes, and a regions of interest file.  It then uploads the annotation file and creates a list of transcripts
	 * to produce proteins from.  These transcripts are stored, compared, put into an output file and have data/statistics run on them.  The synthesis method has several parameters that
	 * can be set to affect how to  and performs protein synthesis.
	 * 
	 */
	public Annotation(String annotationFile, String refchrmDirectory,String chrmDirectory, String proteinFASTAOutputDir, String regionsOfInterest, String outputFileName, String statsOutputFile, String dataPointsFileName){
		this.proteinFASTAOUTPUTDIR = proteinFASTAOutputDir;
		this.annotationFile = annotationFile;
		this.refchrmDir = refchrmDirectory;
		this.chrmDir = chrmDirectory;
		this.regionsFile = regionsOfInterest;
		this.outputFileName = outputFileName;
		this.statsOutputFile = statsOutputFile;
		this.dataPointsFileName = dataPointsFileName;
		
		File mitoFile = new File(chrmDirectory + "chrM.fa");
		genomeHasMito = mitoFile.exists();
		
	}
	/**
	 * Proteome Lite simply creates a proteome, with no extra information or comparisons done.
	 * @param debug
	 * @param useModifiedStopsStarts
	 * @param preFilterForStartCodon
	 */
	public void proteomeLite(boolean debug , boolean useModifiedStopsStarts, boolean preFilterForStartCodon){
		//Set variables
		displayInfo = debug;
		modStopStart = useModifiedStopsStarts;
		preFilterStarts = preFilterForStartCodon;
		
		//Get the current time
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		//Read in the annotation file, create the transcripts and proteins, the create the output.
		//This is the stopwatch used when not debugging
		if(!displayInfo){
			U.startStopwatch();
		}
		if(displayInfo){
			U.p("Uploading annotation.");
			U.startStopwatch();
		}
		//Parse the GTF file to get every line converted into a GENCODE_GTF_Line object
		populateGTFList(annotationFile);
		U.p("Total number of protein encoding lines parsed: " + GTFlist.size());
		if(displayInfo){
			U.p("Annotation and regions parsed.");
			U.stopStopwatch();
		}
		
		populateChrmArray(chrmDir);
		
		if(displayInfo){
			U.startStopwatch();
		}
		U.p("Assembling transcripts and creating proteins.");
		//Parse the GTF file to get every line converted into a GENCODE_GTF_Line object
		populateTranscriptsAndCDS(modStopStart);
		
		if(displayInfo){
			U.p("Finished creating transcripts and proteins.");
			U.stopStopwatch();
		}

		if(displayInfo){
			U.startStopwatch();
		}

		
		U.p("Comparing proteins and generating output file.");
		createLiteOutput();
		if(displayInfo){
			U.stopStopwatch();
			U.p("Protein files completed");
		}
		if(!displayInfo){
			U.stopStopwatch();
		}
		
		//Get the ending time
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.p("Finishing up: " + endTime);
		
		
	}
	/**
	 * Synthesis uses information about files passed in through the constructor to synthesize proteins.
	 * Synthesis adds one protein into the output folder for each transcript with a start and coding regions passed in the annotation file.
	 * This runs with the default options turned on.  Minimal information is displayed and the synthesis is performed mirroring nature.
	 */
	public void synthesize() {
		//Call synthesize with the default value for modStopStart
		synthesize(displayInfo, modStopStart, preFilterStarts);

	}
	
	/**
	 * Synthesis uses information about files passed in through the constructor to synthesize proteins.
	 * Synthesis adds one protein into the output folder for each transcript with a start and coding regions passed in the annotation file.
	 * This allows for the debugging information to be turned on or off, but the synthesis occurs with default options (Mirror real life).
	 * @param debug displays additional information for each individual component of the method.
	 */
	public void synthesize(boolean debug) {
		//Call synthesize with the default value for modStopStart
		synthesize(debug, modStopStart, preFilterStarts);

	}
	/**
	 * Synthesis uses information about files passed in through the constructor to synthesize proteins.
	 * Synthesis adds one protein into the output folder for each transcript with a start and coding regions passed in the annotation file.
	 * @param debug displays additional information for each individual component of the method.
	 * @param useModifiedStopsStarts this parameter sets whether to use a exploration for stops and starts that are not in the location based on the index.  This mirrors how it would occur in the natural world.
	 * @param preFilterForStopCodon This option determines whether to throw out transcripts that do not have a start 
	 */
	public void synthesize(boolean debug , boolean useModifiedStopsStarts, boolean preFilterForStartCodon) {
		//Set variables
		displayInfo = debug;
		modStopStart = useModifiedStopsStarts;
		preFilterStarts = preFilterForStartCodon;
		
		//Get the current time
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		//Read in the annotation file, create the transcripts and proteins, the create the output.
		//This is the stopwatch used when not debugging
		if(!displayInfo){
			U.startStopwatch();
		}
		if(displayInfo){
			U.p("Uploading annotation and regions of interest.");
			U.startStopwatch();
		}
		//Parse the GTF file to get every line converted into a GENCODE_GTF_Line object
		populateGTFList(annotationFile);
		loadRegionsOfInterest(regionsFile);
		U.p("Total number of protein encoding lines parsed: " + GTFlist.size());
		U.p(regionsOfInterest.size() + " regions of interest.");
		if(displayInfo){
			U.p("Annotation and regions parsed.");
			U.stopStopwatch();
		}
		
		populateChrmArray(refchrmDir);
		
		if(displayInfo){
			U.startStopwatch();
		}
		U.p("Creating reference list of transcripts.");
		//Parse the GTF file to get every line converted into a GENCODE_GTF_Line object
		populateTranscriptsAndCDS(modStopStart);
		
		if(displayInfo){
			U.p("Finished creating transcripts and CDS.");
			U.stopStopwatch();
		}
		
		//Save the results for comparison
		referenceTranscriptList = new ArrayList<Transcript>();
		referenceTranscriptList = transcriptList;
		if(displayInfo){
			U.p("Reference List complete.");
		}
		
		if(displayInfo){
			U.p("Writing the reference genome to disk.");
		}
		//Output the reference genome
		createOutputFiles("ref" + this.outputFileName, this.dataPointsFileName, this.statsOutputFile);
		
		
		//Now recompile the chrmDir
		populateChrmArray(chrmDir);
		
		if(displayInfo){
			U.startStopwatch();
		}
		U.p("Creating transcripts and CDS.");
		//Parse the GTF file to get every line converted into a GENCODE_GTF_Line object
		populateTranscriptsAndCDS(modStopStart);
		
		if(displayInfo){
			U.p("Finished creating transcripts and CDS.");
			U.stopStopwatch();
		}


		
		if(displayInfo){
			U.startStopwatch();
		}

		
		U.p("Comparing proteins and generating output file.");
		createOutputFiles(this.outputFileName, this.dataPointsFileName, this.statsOutputFile);
		if(displayInfo){
			U.stopStopwatch();
			U.p("Protein files completed");
		}
		if(!displayInfo){
			U.stopStopwatch();
		}
		
		//Get the ending time
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.p("Finishing up: " + endTime);
		

	}
	/**
	 * parseGTFFile takes a GTF format annotation file and populates an ArrayList, (GTFlist), of GENCODE_GTF_Line objects representing each line of the file.
	 * @param The path to the Annotation file.
	 */
	private void populateGTFList(String GTFFileLocation){
		File f = new File(GTFFileLocation);
		populateGTFList(f);
	}//populateGTFList
	
	/**
	 * parseGTFFile takes a GTF format annotation file and populates an ArrayList, (GTFlist), of GENCODE_GTF_Line objects representing each line of the file.
	 * @param A file object representing the Annotation file.
	 */
	public void populateGTFList(File GTFFile){
		//GTFint List to work on
		//2588001 is the size of the complete v11 annotation
		GTFlist = new ArrayList<GENCODE_GTF_Line>(2588001);
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

					GTFlist.add(line);
					id++;
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
			GTFlist.trimToSize();
			
		}catch(FileNotFoundException e){
			U.p("Error populating GTF List: " + e);

		}//catch
		
	}//populateGTFList
	
	/**
	 * loadRegionsOfInterest loads in a file containing a score based on a regions amplification.  These regions are stored, and are used to score proteins as they are created.
	 * @param The path of the file that contains the regions of interest.
	 * 
	 * @author Brian Risk, modified by David "Corvette" Thomas
	 */
	private void loadRegionsOfInterest(String regionsFile) {
		try {
			
			//Don't upload any regions of interest if there is not a file input
			if(regionsFile.equals("noregions")){
				return;
			}
			BufferedReader br = new BufferedReader(new FileReader(regionsFile));
			
			//Prime the loop
			String line = br.readLine();
			
			//Iterate through the file as long as lines are found.
			while (line != null) {
				//Delimit with tabs
				String [] chunks = line.split("\t");
				
				//A variable to store the chromosome file.  These are the names used in 'Definitions'
				int chromosomeName;
				//Convert the chromosome given in each line into the integer representation in the definitions file
				if(chunks[0].equalsIgnoreCase("chrM") || chunks[0].equalsIgnoreCase("M")){
					chromosomeName = Definitions.chromosomeM;
				}else if(chunks[0].equalsIgnoreCase("chrX") || chunks[0].equalsIgnoreCase("X")){
					chromosomeName = Definitions.chromosomeX;
				}else if(chunks[0].equalsIgnoreCase("chrY") || chunks[0].equalsIgnoreCase("Y")){
					chromosomeName = Definitions.chromosomeY;
				}else{
					chromosomeName = Integer.parseInt(chunks[0].substring(chunks[0].indexOf('r') + 1));
					
				}
				//Store the region into the regionsOfInterest List
				regionsOfInterest.add(
						new RegionOfInterest(
							chromosomeName,
							Integer.parseInt(chunks[1]),
							Integer.parseInt(chunks[2]),
							Double.parseDouble(chunks[3])
						)
					);
				// Don't forget to load in that next line!
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			U.p("No regions of interest file was specified.");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}//Load regions of interest

	/**
	 * Creates a master list of Transcripts and CDS.  Each transcripts creates it protein after it is made.  The transcripts are stored in transcirptList.
	 * It creates transcripts based on the files loaded into the chromosome directory.  
	 */
	public void populateTranscriptsAndCDS(boolean useModStartStop){

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
			if(displayInfo){
				U.p("Working on " + Definitions.convertChrmNumToString(k + 1) + ".fa") ;
			}
			seq = new Sequence_DNA(chrmFile[k]);
			
//			StringBuffer sb = new StringBuffer();
//			sb.append(upLoadFile(chrmFile[k]));
			
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
//						sequence = sb.substring(beingTested.getStartLocation() - 1,beingTested.getStopLocation());
					}else{
						sequence = seq.getNucleotideSequences().get(0).getSequence().substring(beingTested.getStartLocation(),beingTested.getStopLocation() + 1);
//						sequence = sb.substring(beingTested.getStartLocation(),beingTested.getStopLocation() + 1);
					}
					
					//Construct a basic transcript to be added to the list
					t = new Transcript(beingTested.getChromosomeName(), sequence, beingTested.getGenomicStrand() , transcripts.get(j), beingTested.getStartLocation(), beingTested.getStopLocation(), useModStartStop);


					//Prime the loop
					//c simply stores each line as it is iterated through
					GENCODE_GTF_Line c = GTFlist.get(transcripts.get(j) + 1);
					boolean startFound = false;
					boolean stopFound = false;
					
					//preFilter starts determines whether to throw out transcripts that do not contain a start codon
					
					//Loops through each line after this transcript, until another transcript is encountered
					for(int x = transcripts.get(j) + 1; x < GTFlist.size() && c.getFeatureType() != Definitions.featureTypeTRANSCRIPT; x++){
						c = GTFlist.get(x);

						
						//If a start codon is found, then set startFound to true
						if(c.getFeatureType() == Definitions.featureTypeSTART_CODON){
							startFound = true;
						}
						//If a start codon is found, then set startFound to true
						if(c.getFeatureType() == Definitions.featureTypeSTOP_CODON){
							stopFound = true;
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
					//Confirm that the transcript has a CDS to code from
					if(t.getCDS().size() > 0){

						
						//Calculate the amplification Score and add it to the transcript
						int proteinStart = t.getCDS().get(0).getStart();
						int proteinEnd = t.getCDS().get(t.getCDS().size() - 1).getStop();
						double highestScore = 0;
						for (RegionOfInterest regionOfInterest: regionsOfInterest) {
							if(regionOfInterest.getChromosome() - 1 != k){
								continue;
							}
							
							//Ensure that if a region of Interest is being compared against this transcript, they are on the same chromosome
							assert (beingTested.getChromosomeName() == regionOfInterest.getChromosome());

							if(proteinStart >= regionOfInterest.getStart()  && proteinStart <= regionOfInterest.getStop()){
								if(Math.abs(regionOfInterest.getScore()) > Math.abs(highestScore)){
									highestScore = regionOfInterest.getScore();
								}
							}
							if(proteinEnd >= regionOfInterest.getStart()  && proteinEnd <= regionOfInterest.getStop()){
								if(Math.abs(regionOfInterest.getScore()) > Math.abs(highestScore)){
									highestScore = regionOfInterest.getScore();
								}
							}
							if(proteinEnd >= regionOfInterest.getStop()  && proteinStart <= regionOfInterest.getStart()){
								if(Math.abs(regionOfInterest.getScore()) > Math.abs(highestScore)){
									highestScore = regionOfInterest.getScore();
								}
							}
							
						}
						
						//Set the amplification score
						t.setAmpScore(highestScore);
					
						
						t.setHasStartCodon(startFound);
						t.setHasStopCodon(stopFound);
						
						//Create t's protein and add it to the list
						t.createProtein();

						

						

						//Now that t has a score, has CDS regions, and has been confirmed to have a start, add it to the list of transcripts to produce proteins with
						transcriptList.add(t);
						addCount++;
					}


					
				}//if
				
			}//inner for
			
		}//outer for
		
	}//populateTranscriptsAndCDS
		
	
	/**
	 * populateTranscriptsAndCDSNoAmino is a very fast method for converting a read in GTF file into a bare bones arrayList of Transcript objects.  This
	 * allows for use of the transcript objects methods without all of the overhead of creating proteins.
	 */
	public void populateTranscriptsAndCDSNoProtein(){

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
		//Iterate through each chromosome file, and create all of the transcripts possible from each one
		//Code is performed in this order to minize HDD reads.

		GENCODE_GTF_Line beingTested;
		
		//Check each transcript to determine if it is using the currently open file.
		for(int j = 0; j < transcripts.size(); j++){
			
			//Get the currently used transcripts line
			beingTested = GTFlist.get(transcripts.get(j));
			

				
				//sequence is the DNA strand from the chromosome file for this transcript
				String sequence;
				StringBuffer sb = new StringBuffer();
				for(int k = 0; k < (beingTested.getStopLocation() + 1 - beingTested.getStartLocation()); k++){
					sb.append(Definitions.NULL_CHAR);
				}//for
				
				
				sequence = sb.toString();
				
				
				//Construct a basic transcript to be added to the list
				t = new Transcript(beingTested.getChromosomeName(), sequence, beingTested.getGenomicStrand() , transcripts.get(j), beingTested.getStartLocation(), beingTested.getStopLocation(), false);


				//Prime the loop
				//c simply stores each line as it is iterated through
				GENCODE_GTF_Line c = GTFlist.get(transcripts.get(j) + 1);
				boolean startFound = false;
				boolean stopFound = false;
				
				//preFilter starts determines whether to throw out transcripts that do not contain a start codon
				
				//Loops through each line after this transcript, until another transcript is encountered
				for(int x = transcripts.get(j) + 1; x < GTFlist.size() && c.getFeatureType() != Definitions.featureTypeTRANSCRIPT; x++){
					c = GTFlist.get(x);

					
					//If a start codon is found, then set startFound to true
					if(c.getFeatureType() == Definitions.featureTypeSTART_CODON){
						startFound = true;
					}
					//If a start codon is found, then set startFound to true
					if(c.getFeatureType() == Definitions.featureTypeSTOP_CODON){
						stopFound = true;
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
				//Confirm that the transcript has a CDS to code from
				if(t.getCDS().size() > 0){
		
					
					t.setHasStartCodon(startFound);
					t.setHasStopCodon(stopFound);
					
					//Create t's protein and add it to the list

					t.createProtein();

					

					//Now that t has a score, has CDS regions, and has been confirmed to have a start, add it to the list of transcripts to produce proteins with
					transcriptList.add(t);
					addCount++;
				}


				
		
			
		}//inner for
		

		
	}//populateTranscriptsAndCDS
	
	/**
	 * createOutputFiles adds each genes protein to a fasta format file and computes statistics about the transcripts.  The proteins to output are stored in transcriptList.  The proteins are 
	 * compared to the proteins from the referenceTranscriptList for differences, and those variances are included in the header information for each FASTA entry.
	 * @param outputFileName Name of the file(including file type) to store the output in.
	 * @param dataPointsFileName Name of the file(including file type) to store the data points.  Each 'DataPoint' represents a protein with a non-zero amplification score and non-zero variant count.
	 * @param statsOutputFile 
	 */
	private void createOutputFiles(String outputFileName, String dataPointsFileName, String statsOutputFile){
		try {
			

				
	
				
				/*Prep the output to compare the difference between the reference and test genome*/
	
				//List of data points to output to the data file
				ArrayList<DataPoint> variantAmpData;
				
				
				//Writing to HDD related variables
				StringBuffer sb;
				BufferedWriter out;
				File outFile;
				
				//Data holding variables
				Transcript t;
				String geneName;
				String transcriptName;
				String protein;
				int numOfVariants;
				double score;
				int endExtensionLength;
				int startExtensionLength;
				
				//Editing related variables
				int linesOfProtein;
				
				

				
				//Instantiate variables used in the outputing of files
				variantAmpData = new ArrayList<DataPoint>();
				outFile = new File(proteinFASTAOUTPUTDIR + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");
				outFile.mkdir();
				String outputDir = outFile.getAbsolutePath() + "/";
				out = new BufferedWriter(new FileWriter(outputDir + outputFileName));

				sb = new StringBuffer();
				
				//Iterate through all of the transcripts and clean up/output those proteins.
				for(int i = 0; i < transcriptList.size(); i++){
	
					t = transcriptList.get(i);
					/*Compare the newly computed proteins with the reference protein.  Set any differences to lower case*/
					protein = t.getProtein();
					String refProtein = referenceTranscriptList.get(i).getProtein();
					
					
					//If a GENE produces a protein that has 0 length, then ignore it.
					if(protein.equals("PROTEIN TOO SHORT")){
						continue;
					}
					
					//Determine which acids are different from the reference genome.
					int limit = 0;
					if(protein.length() > refProtein.length()){
						limit = refProtein.length();
					}else if(refProtein.length() > protein.length()){
						limit = protein.length();
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
							//Increase this proteins variant count
							t.incrementVariantCount();
							//Insert the new lower case character
							String charToLower = protein.substring(n, n + 1);
							charToLower = charToLower.toLowerCase();
							String start = protein.substring(0, n);
							String end = protein.substring(n + 1);
							protein = start + charToLower + end;
							
						}
						
					}
					
					
					protein = removeModChar(protein);
					
					if(protein.length() == 0){
						//A start has ate up the entire protein, so don't bother adding it to the output
						continue;			
					}
					
					//If transcripts without starts are being prefiltered out, then skip ones without starts
					if(preFilterStarts){
						if(!t.isHasStartCodon()){
							continue;
						}
					}
					
					
					
					//Assemble information about a transcript before writing it to file
					geneName = GTFlist.get(t.getiD()).getGene_Name();
					transcriptName = GTFlist.get(t.getiD()).getTranscript_Name();
					score = t.getAmpScore();
					numOfVariants = t.getvariantCount();
					endExtensionLength = t.getEndExtensionLength();
					startExtensionLength = t.getStartExtensionLength();

					//Calculate how many lines this will take up in the fasta file
					linesOfProtein = protein.length() / lineLength;
		
					//Collect information about this transcript so it can be used later in a data file
					DataPoint currentDP = new DataPoint("> " + geneName + "|" + transcriptName +"|", score, numOfVariants);
					if(currentDP.getVariantCount() > 0){
						variantAmpData.add(currentDP);
					}
					
	
					//Prime the file with header information for this transcript
					sb.append("> " + geneName + "|" + transcriptName +"|" + protein.length() + "|" + "AmpScore: " + score +"|" + "#ofVariants: "+ numOfVariants + "|");
					//Append information about modified stops and starts if it is enabled
					if(modStopStart){
						 sb.append("StartEx: " + startExtensionLength + "|" + "EndEx: " + endExtensionLength + "|");
					}
					sb.append("\n");
					sb.append("> " + "HasPrematureStop: " + t.doesContainPreMatureStop() + "|" + "PrematureStopLength: " + t.getPreMatureStopCount() + "|");
					sb.append("\n");
					sb.append("> " + GTFlist.get(t.getiD()).getTranscriptID() + "|" + GTFlist.get(t.getiD()).getGeneID() + "|");
					sb.append("\n");
					//Write this transcript to the file.
					for(int k = 0; k < linesOfProtein + 1; k++){
						//If it is the last line, just write all of them
						if(k == linesOfProtein){
								sb.append(protein.substring(k*lineLength));
								sb.append("\n");
							
						}else{
						//Write 80 characters
						sb.append(protein.substring(k*lineLength, (k+1)*lineLength));
						sb.append("\n");
						}
						
					}//for
						sb.append("\n");
						
						
				}//Large for	
				
				out.write(sb.toString());
				
			//Flush and close the buffered writer
			out.flush();
			out.close();
			

			
			//Output the data points file
			//Sort it from smallest to greatest
			Collections.sort(variantAmpData);
			//Create a output file
			BufferedWriter dataPointsOut= new BufferedWriter(new FileWriter(outputDir + dataPointsFileName));
			for(DataPoint dp: variantAmpData){
				dataPointsOut.write(dp.toString());
				dataPointsOut.newLine();
			}

			//Calculate statistics
			//Close the output file
			dataPointsOut.flush();
			dataPointsOut.close();
			
			
			//Calculate the statistics
			/* ***** STATS ***** */
			//Error Variables
			int refDotCount = 0;
			int refMCount = 0;
			int bothTotal = 0;
			
			//Odd occurances in protein variables
			int preMatureRefStop = 0;
			int refToShort = 0;
			
			//Variables for Starts
			int refStartExtCount = 0;
			int refStartExtTotal = 0;
			ArrayList<Integer> refStartMax = new ArrayList<Integer>();
			
			//The total Number of deleted proteins is this list size
			ArrayList<String> deletedRefProteinNames = new ArrayList<String>();
			
			int refEndExtCount = 0;
			int refEndExtTotal = 0;
			ArrayList<Integer> refEndMax = new ArrayList<Integer>();
			
			int refSelenoCount = 0;
			int refTotalSelenos = 0;
			int refProteinWithSelenoCount = 0;
			
			int refLeftOverCount = 0;
			
			int refNonZeroAmpCount = 0;
			
			//Mutations are stored in a linked list in reverse order, because errors are likely to occur back to back, so it will save look up time when comparing them to keep the most recent one at the front.
			LinkedList<Mutation> refMutations = new LinkedList<Mutation>();
			
			
			//Collect error information on the reference genome
			for(Transcript tran: referenceTranscriptList){
				
				if(tran.getAmpScore() != 0){
					refNonZeroAmpCount++;
				}
				
				String tempRefProtein = tran.getProtein();
				if(tempRefProtein.equals("PROTEIN TOO SHORT")){
					refToShort++;
					continue;
				}
				
				if(tran.getLeftOverAminoLength() > 0){
					refLeftOverCount++;
				}

				//Clean up a protein
				tempRefProtein = removeModChar(tempRefProtein);
				//Statistics about modified starts and stops
				
				//Look to see if a protein was deleted by a start
				if(tempRefProtein.length() == 0){
					//A Start has been extended all the way to the end of the protein
					deletedRefProteinNames.add(GTFlist.get(tran.getiD()).getGene_Name() + " " + GTFlist.get(tran.getiD()).getTranscript_Name());
				}
				
				if(tran.getSeleno().size() > 0){
					refProteinWithSelenoCount++;
					refTotalSelenos+= tran.getSeleno().size();
					boolean alreadyExists = false;
					for(int i = 0; i < tran.getSeleno().size() && alreadyExists == false; i++){
						for(Mutation m: refMutations){
							if(m.equals(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_SELENOCYSTEINE, tran.getSeleno().get(i).getStart())){
									alreadyExists = true;
									break;
							}
						}
					}
					if(!alreadyExists){
						refSelenoCount += tran.getSeleno().size();
						//Add this mutation to the list of mutations
						for(int i = 0; i < tran.getSeleno().size(); i++){
							refMutations.addFirst(new Mutation(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_SELENOCYSTEINE, tran.getSeleno().get(i).getStart()));
						}
					}
				}
				
				//Look at the situation where a start has been extended
				if(tran.getStartExtensionLength() > 0){
					boolean alreadyExists = false;
					for(Mutation m: refMutations){
						if(m.equals(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_DELAYED_START, tran.getStartErrorLoc())){
								alreadyExists = true;
								break;
						}
					}
					if(!alreadyExists){
						refStartExtCount++;
						refStartExtTotal += tran.getStartExtensionLength();
						refStartMax.add(tran.getStartExtensionLength());
						//Add this mutation to the list of mutations
						refMutations.addFirst(new Mutation(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_DELAYED_START, tran.getStartErrorLoc()));
					}
				}
				
				//Look at the situation where a stop has been extended
				if(tran.getEndExtensionLength() > 0 && !tran.doesContainPreMatureStop()){
					
					boolean alreadyExists = false;
					for(Mutation m: refMutations){
						if(m.equals(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_EXTENDED_STOP, tran.getStartErrorLoc())){
								alreadyExists = true;
								break;
						}
					}
					if(!alreadyExists){
						refEndExtCount++;
						refEndExtTotal += tran.getEndExtensionLength();
						refEndMax.add(tran.getEndExtensionLength());
						//Add this mutation to the list of mutations
						refMutations.addFirst(new Mutation(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_EXTENDED_STOP, tran.getStartErrorLoc()));
					}
				}
				//Count premature stops
				if(tran.doesContainPreMatureStop()){
					preMatureRefStop++;
				}
				
				//Error Data
				if(tempRefProtein.contains(".")){
					refDotCount++;
				}
				
				//ensure that the length does not cause a issue in deleted proteins.
				if(tempRefProtein.length() != 0){
					if(tempRefProtein.charAt(0) != 'M'){
						refMCount++;
					}
					if(tempRefProtein.charAt(0) != 'M' && tempRefProtein.contains(".")){
						bothTotal++;
					}
				}
				
			}//For each loop

			//Variance Variables
			int isVariantCount = 0;
			int variantCount = 0;
			int oneCount = 0;
			int max = 0;
			int nonZeroAmpCount = 0;
			
			//Error variables
			int mCount = 0;
			int dotCount = 0;
			int error = 0;
			
			//Odd occurences variables
			int preMatureStop = 0;
			int toShort = 0;
			
			//Variables for Starts
			int startExtCount = 0;
			int startExtTotal = 0;
			ArrayList<Integer> startMax = new ArrayList<Integer>();
			
			//The total Number of deleted proteins is this list size
			ArrayList<String> deletedProteinNames = new ArrayList<String>();
			
			int endExtCount = 0;
			int endExtTotal = 0;
			ArrayList<Integer> endMax = new ArrayList<Integer>();
			
			int proteinWithSelenoCount = 0;
			int selenoCount = 0;
			int totalSelenos = 0;
			
			int leftOverCount = 0;
			
			//Mutations are stored in a linked list in reverse order, because errors are likely to occur back to back, so it will save look up time when comparing them to keep the most recent one at the front.
			LinkedList<Mutation> mutations = new LinkedList<Mutation>();
			
			//Collect statistics on the non-reference genome
			for(Transcript tran: transcriptList){
				int vc = tran.getvariantCount();
				//isVariant
				if(vc > 0){
					isVariantCount++;
				}
				
				//Total number of variants
				variantCount += vc;
				
				//number of proteins with one variant
				if(vc == 1){
					oneCount++;
				}
				
				//Maximum variant
				if(vc > max){
					max = vc;
				}
				
				if(tran.getAmpScore() != 0){
					nonZeroAmpCount++;
				}
				
				if(tran.getProtein().equals(Definitions.PROTEIN_TYPE_TOO_SHORT)){
					toShort++;
					continue;
				}
				

				if(tran.getLeftOverAminoLength() > 0){
					leftOverCount++;
				}
				
				//Clean up a protein
				String tempProtein = removeModChar(tran.getProtein());
				
				//Look to see if a protein was deleted by a start
				if(tempProtein.length() == 0){
					//A Start has been extended all the way to the end of the protein
					deletedProteinNames.add(GTFlist.get(tran.getiD()).getGene_Name() + " " + GTFlist.get(tran.getiD()).getTranscript_Name());
				}
				
				//Count Selenocysteine
				if(tran.getSeleno().size() > 0){
					proteinWithSelenoCount++;
					totalSelenos += tran.getSeleno().size();
					boolean alreadyExists = false;
					for(int i = 0; i < tran.getSeleno().size() && alreadyExists == false; i++){
						for(Mutation m: mutations){
							if(m.equals(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_SELENOCYSTEINE, tran.getSeleno().get(i).getStart())){
									alreadyExists = true;
									break;
							}
						}
					}
					if(!alreadyExists){
						selenoCount += tran.getSeleno().size();
						//Add this mutation to the list of mutations
						for(int i = 0; i < tran.getSeleno().size(); i++){
							mutations.addFirst(new Mutation(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_SELENOCYSTEINE, tran.getSeleno().get(i).getStart()));
						}
					}
				}
				
				//Consider modified starts
				if(tran.getStartExtensionLength() > 0 && !tran.doesContainPreMatureStop()){
					boolean alreadyExists = false;
					for(Mutation m: mutations){
						if(m.equals(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_DELAYED_START, tran.getStartErrorLoc())){
								alreadyExists = true;
								break;
						}
					}
					if(!alreadyExists){
						startExtCount++;
						startExtTotal += tran.getStartExtensionLength();
						startMax.add(tran.getStartExtensionLength());
						//Add this mutation to the list of mutations
						mutations.addFirst(new Mutation(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_DELAYED_START, tran.getStartErrorLoc()));
					}

				}
				
				//Look at the situation where a stop has been extended
				if(tran.getEndExtensionLength() > 0){
					boolean alreadyExists = false;
					for(Mutation m: mutations){
						if(m.equals(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_EXTENDED_STOP, tran.getStartErrorLoc())){
								alreadyExists = true;
								break;
						}
					}
					if(!alreadyExists){
						endExtCount++;
						endExtTotal += tran.getEndExtensionLength();
						endMax.add(tran.getEndExtensionLength());
						mutations.addFirst(new Mutation(GTFlist.get(tran.getiD()).getGene_Name(), Definitions.ERROR_EXTENDED_STOP, tran.getStartErrorLoc()));
					}

				}
				
				if(tran.doesContainPreMatureStop()){
					preMatureStop++;
				}

				
				//Error Data
				if(tempProtein.contains(".")){
					dotCount++;
				}
				//ensure that the length does not cause a issue in deleted proteins.
				if(tempProtein.length() != 0){
					if(!(tempProtein.charAt(0) == 'M' || tempProtein.charAt(0) == 'X')){
						mCount++;
						error++;
					}
					if(tempProtein.contains(".")){
						error++;
					}
				}
			}//For each loop
			
			//Sort all of the extensions and starts
			Collections.sort(refStartMax);
			Collections.sort(refEndMax);
			Collections.sort(startMax);
			Collections.sort(endMax);
			
			/*Mark the finishing time*/
			endCal = Calendar.getInstance();

			//Output the statsitics to a file
			BufferedWriter statsOut= new BufferedWriter(new FileWriter(outputDir + statsOutputFile));
			//Calculate statistics
			statsOut.write("Start time: " + startTime + "\n");
			statsOut.write("End time:   " + sdf.format(endCal.getTime()) + "\n");
			statsOut.newLine();
			statsOut.write("File Information" + "\n");
			statsOut.write("Annotation file:     " + Properties.annotationFile + "\n");
			statsOut.write("Reference Genome:    " + Properties.refGenomeDir + "\n");
			statsOut.write("Genome Directory:    " + Properties.genomeDir + "\n");
			statsOut.write("Regions of interest: " + Properties.regionsOfInterestFile + "\n");
			statsOut.write("Prefiltering Option: ");
			if(Properties.preFilterOutTranWithoutStartCodon){
				statsOut.write("On");
			}else{
				statsOut.write("Off");
			}
			statsOut.newLine();
			
			statsOut.write("Start/Stop Modification: ");
			if(Properties.useModifiedStopsAndStarts){
				statsOut.write("On");
			}else{
				statsOut.write("Off");
			}
			statsOut.newLine();
			
			
			statsOut.newLine();
			statsOut.write("*******************   STATS   ***********************" + "\n");
			statsOut.write("Totals Information" + "\n");
			statsOut.write("Annotation lines associated with protein coding: " + GTFlist.size() + "." + "\n");
			statsOut.write("Total number of proteins created:     " + (transcriptList.size() - refToShort) + " out of " + transcriptList.size() + "." + "\n");
			statsOut.write("Regions of interest count:            " + regionsOfInterest.size() + "." + "\n");

			statsOut.newLine();
			statsOut.newLine();
			statsOut.write("Output file error information" + "\n");
			statsOut.write("Output errors total:          " + error + "\n");
			statsOut.write("X or M not found at start:    " + mCount + "\n");
			statsOut.write("Stop codon inside of protein: " + dotCount + "\n");
			statsOut.newLine();
			statsOut.newLine();
			statsOut.write("Premature stop totals"  + "\n");
			statsOut.write("Reference genome: " + preMatureRefStop + "\n");
			statsOut.write("Input genome:    " + preMatureStop  + "\n");
			statsOut.newLine();
			statsOut.newLine();
			statsOut.write("Regions of Interest"  + "\n");
			statsOut.write("Reference genome non zero score total: " + refNonZeroAmpCount + "\n");
			statsOut.write("Input genome non zero score total:     " + nonZeroAmpCount  + "\n");
			statsOut.newLine();
			statsOut.newLine();
			statsOut.write("Selenocysteine" + "\n"); 
			statsOut.write("Reference Genome****" + "\n");
			statsOut.write("Locations:           " + refSelenoCount + "\n");
			statsOut.write("Total found:         " + refTotalSelenos + "\n");
			statsOut.write("Proteins containing: " + refProteinWithSelenoCount + "\n");
			statsOut.write("Genome**************" + "\n");
			statsOut.write("Locations:            " + selenoCount + "\n");
			statsOut.write("Total found:          " + totalSelenos + "\n");
			statsOut.write("Proteins containing: " + proteinWithSelenoCount + "\n");
			statsOut.newLine();
			statsOut.write("NOTE: While computing these statistics, if a error occurred on more then one protein but stemmed from the same location in the DNA, it is only counted the first time.");
			statsOut.newLine();
			if(modStopStart){
				statsOut.newLine();
				statsOut.write("Start/Stop Modifications.");
				statsOut.newLine();
				statsOut.newLine();
				statsOut.write("Start Modifications" + "\n");
				statsOut.write("Reference Genome****" + "\n");
				statsOut.write("Delayed start total:       " + refStartExtCount + "\n");
				statsOut.write("Total amino acids skipped: " + refStartExtTotal + "\n");
				if(refStartMax.size() > 0){
					statsOut.write("Maxed amino acids skipped: " + refStartMax.get(refStartMax.size() - 1) + "\n");
					statsOut.write("Median skipped amino acids: " + refStartMax.get(refStartMax.size() / 2) + "\n");
				}
				statsOut.write("Average number of skipped amino acids: " + ((double)refStartExtTotal / (double)refStartExtCount) + "\n");
				if(deletedRefProteinNames.size() > 0){
					statsOut.write("Deleted proteins because of no start total: " + deletedRefProteinNames.size());
					for(int i = 0; i < deletedRefProteinNames.size(); i++){
						if( i % 5 == 0){
							statsOut.newLine();
						}
						if(i == deletedRefProteinNames.size() - 1){
							statsOut.write(deletedRefProteinNames.get(i) + ".");
							statsOut.newLine();
							break;
						}
						statsOut.write(deletedRefProteinNames.get(i) + ", ");
						
					}//for
				}//deletedNames if
				statsOut.write("Genome**************" + "\n");
				statsOut.write("Delayed start total:       " + startExtCount + "\n");
				statsOut.write("Total amino acids skipped: " + startExtTotal + "\n");
				if(startMax.size() > 0){
					statsOut.write("Maxed amino acids skipped: " + startMax.get(startMax.size() - 1) + "\n");
					statsOut.write("Median skipped amino acids: " + startMax.get(startMax.size() / 2) + "\n");
				}
				statsOut.write("Average number of skipped amino acids: " + ((double)startExtTotal / (double)startExtCount) + "\n");
				if(deletedProteinNames.size() > 0){
					statsOut.write("Deleted proteins because of no start total: " + deletedProteinNames.size());
					for(int i = 0; i < deletedProteinNames.size(); i++){
						if( i % 5 == 0){
							statsOut.newLine();
						}
						if(i == deletedProteinNames.size() - 1){
							statsOut.write(deletedProteinNames.get(i) + ".");
							statsOut.newLine();
							break;
						}
						statsOut.write(deletedProteinNames.get(i) + ", ");
						
					}//for
				}//deletedNames if
				statsOut.newLine();
				statsOut.newLine();
				statsOut.write("Stop modifications" + "\n");
				statsOut.write("Reference Genome****" + "\n");
				statsOut.write("Extended protein total: " + refEndExtCount + "\n");
				statsOut.write("Skipped amino acid total: " + refEndExtTotal + "\n");
				if(refEndMax.size() > 0){
					statsOut.write("Max amino acids skipped:  " + refEndMax.get(refEndMax.size() - 1) + "\n");
					statsOut.write("Median skipped amino acids: " + refEndMax.get(refEndMax.size() / 2) + "\n");
				}
				statsOut.write("Average number of skipped amino acids: " + ((double)refEndExtTotal / (double)refEndExtCount) + "\n");
				//Median
				statsOut.write("Genome**************" + "\n");
				statsOut.write("Extended protein total: " + endExtCount + "\n");
				statsOut.write("Skipped amino acid total " + endExtTotal + "\n");
				if(endMax.size() > 0){
					statsOut.write("Max amino acids skipped:  " + endMax.get(endMax.size() - 1) + "\n");
					statsOut.write("Median skipped amino acids: " + endMax.get(endMax.size() / 2) + "\n");
				}
				statsOut.write("Average number of skipped amino acids: " + ((double)endExtTotal / (double)endExtCount) + "\n");

			}
			statsOut.newLine();
			statsOut.newLine();
			statsOut.write("Variances in the input genome are changed to lower case." + "\n");
			statsOut.write("Variant protein total:               " + isVariantCount + "\n");
			statsOut.write("Variations total:                    " + variantCount + "\n");
			statsOut.write("Proteins with 1 variance total:      " + oneCount + "\n");
			statsOut.write("Max variations in a single protein:  "  + max + "\n");
			statsOut.write("Average # of variants in proteins with a variant: " + ((double)variantCount/(double)isVariantCount) + "\n");
			statsOut.write("Average # of variants in proteins with more then 1 variant: " +  (((double)(variantCount - oneCount) / (double)(isVariantCount - oneCount))) + "\n");
			double totalProteins = transcriptList.size() - toShort;
			statsOut.write("Percent of proteins with a variance: " + ((isVariantCount / totalProteins) * 100) + "%" + "\n");
			statsOut.newLine();
			statsOut.newLine();
			statsOut.write("*******************   STATS   ***********************" + "\n");
			
			statsOut.flush();
			statsOut.close();
				
		} catch (IOException e) {
			e.printStackTrace();
		}

	
	}//createOutputFiles
	private void createLiteOutput(){
		//Writing to HDD related variables
		BufferedWriter out;
		
		try {
			U.p(this.chrmDir);
			String temp = this.chrmDir.substring(0, this.chrmDir.length() - 1);
			int tempIndex = temp.lastIndexOf('/');
			temp = temp.substring(tempIndex + 1);
			out = new BufferedWriter(new FileWriter(this.proteinFASTAOUTPUTDIR + temp + ".txt"));
			
			
			out.write("##Total number of proteins: " + transcriptList.size() + "\n");
			
			for(int i = 0; i < transcriptList.size(); i++){
				Transcript t = transcriptList.get(i);
				GENCODE_GTF_Line g = GTFlist.get(t.getiD());
				out.write(">" + g.getTranscript_Name() + "|" + g.getGene_Name() + "|" + "Length: " + t.getProtein().length() + "|" + Definitions.convertChrmNumToString(g.getChromosomeName()) + "|" + "\n");
			
				int linesOfProtein = t.getProtein().length() / lineLength;
				//Write this transcript to the file.
				for(int k = 0; k < linesOfProtein + 1; k++){
					//If it is the last line, just write all of them
					if(k == linesOfProtein){
							out.write(t.getProtein().substring(k*lineLength));
							out.write("\n");
						
					}else{
					//Write 80 characters
					out.write(t.getProtein().substring(k*lineLength, (k+1)*lineLength));
					out.write("\n");
					}
					
				}//for
			
				out.write("\n");
			
			}//for
			
			out.flush();
			out.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	/**
	 * removeModChar cleans up proteins with the modification delimiting marks found in 'Definitions'.
	 * It cleaves for premature stops, removes extended end characters, and padding for a delayed start.
	 * @param modProtein The protein with the modifications.
	 * @return The protein minus the modifications.
	 */
	private String removeModChar(String modProtein){
		//If modified stops are starts are disabled, just ignore this call
		if(!Properties.useModifiedStopsAndStarts){
			return modProtein;
		}
		//If a premature stop is found, remove everything after that point
		int prematureStopLoc;
		if((prematureStopLoc = modProtein.indexOf(Definitions.STOP_CODON_STRING)) != -1){
			modProtein = modProtein.substring(0, prematureStopLoc);
		}
		
		//Remove the + from any protein before putting it into the output file
		if(modProtein.contains(Character.toString(Definitions.END_NOT_FOUND_CHAR))){
			String before = modProtein.substring(0, modProtein.indexOf(Definitions.END_NOT_FOUND_CHAR));
			String after = modProtein.substring(modProtein.indexOf(Definitions.END_NOT_FOUND_CHAR) + 1);
			modProtein = before + after;
		}
		
		//Remove the & from any protein before putting it into the output file
		//Default assume no start will be found
		int firstNonStartDelay = modProtein.length();
		for(int i = 0; i < modProtein.length(); i++){
			if(modProtein.charAt(i) != Definitions.START_NOT_FOUND_CHAR){
				firstNonStartDelay = i;
				break;
			}
			
		}
		//Trim off the leading start not found characters
		modProtein = modProtein.substring(firstNonStartDelay, modProtein.length());
		
		//Ensure that the finished protein has a start of M, otherwise to much was trimmed off the front
		if(modProtein.length() > 0){
			assert (modProtein.charAt(0) == 'M' || modProtein.charAt(0) == 'm'); 
		}
		
		return modProtein;
	}
	/**
	 * Populates the chromosome array with the file locations of each chromosome.'
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
	 * 
	 * @return the lines of the annotation file.
	 */
	public ArrayList<GENCODE_GTF_Line> getAnnotaitonLines(){
		return GTFlist;
	}

		
	/**
	 * 
	 * @return the lines of the annotation file.
	 */
	public ArrayList<Transcript> getTranscripts(){
		return transcriptList;
	}
	
}
