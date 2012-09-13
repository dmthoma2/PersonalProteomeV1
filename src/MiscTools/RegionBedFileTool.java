package MiscTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Scanner;

import MiscTools.RegionBedFileToolObjects.BedLine;
import PersonalProteome.Definitions;
import PersonalProteome.U;

/**
 * RegionBedFileTool is a bed file creation tool, that is designed to create a output file that highlights areas of amplification through color.  It takes a directory containing bed
 * files as input (EX FORMAT: chr1	23	52	0.71  '#' denotes a comment).  THe line information is chromosome name, start location, stop location, and amplification score.  After
 * reading in all of the input lines, information is calculated based on the range of amplification scores and each line is assigned a color.  Output lines are of the format:
 * ChromosomeName StartLoc	StopLoc	FileName-Amplification Score Strand StartLoc, StopLoc, Color
 * EXAMPLE: chr1	35907778	41545288	WHIM2-Amplification	1000	+	35907778	41545288	0,80,0
 *                                             Amplification Score
 *    <----------------------"-"-------------------------0----------------------------"+"------------------------>
 * Bright red                                          Black                                                Bright green
 * 
 * @author David "Corvette" Thomas
 *
 */
public class RegionBedFileTool {

	
	
	public static void main(String[] args){
		
		
		String args0 = "/Users/davidthomas/Peppy/synthesis/MiscTools/BedFileCreator/";
		String args1 = "/Users/davidthomas/Peppy/synthesis/MiscTools/BedFileCreator/output/";
		
		//Get a list of all bed files in the input directory
		File jobsDir = new File(args0);
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		if (potentialJobsFiles != null) {
			for (int i = 0; i < potentialJobsFiles.length; i++) {
				if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".bed")) {
					jobFiles.add(potentialJobsFiles[i]);
				}//if	
			}//for
		}//if
		
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
		String startTime = "";
		String endTime = "";
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		
		
		//Generate a output file for each input file
		for(int j = 0; j < jobFiles.size(); j++){
			U.p("Working on file " + (j + 1) + " of " + jobFiles.size());
			RegionBedFileTool rbft = new RegionBedFileTool(jobFiles.get(j).getAbsolutePath(), args1);
		}
		
	
		//Get the ending time
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.p("Finishing up: " + endTime);
		
		
	}//Main

	
	//File information
	String inputFile;
	String outputDir;
	
	
	ArrayList<BedLine> bedLines;
	
	
	/**
	 * RegionBedFile forms a bed file based on the input files.  It creates color coden regions based on the amplification score passed in.
	 * @param inputFile  tab delimited bed file of the format(chromosomeName startLoc stopLoc amplificationScore)
	 * @param outputDir  Directory to output the finished bed file to.
	 */
	public RegionBedFileTool(String inputFile, String outputDir){
		this.inputFile = inputFile;
		this.outputDir = outputDir;
		
		createBedFile();
	}
	
	/**
	 * createBedFile is the method that handles uploading input, calculating statistics, and creating output files.
	 */
	private void createBedFile(){
		
		//Read Input in as a object
		uploadInput();
		//Create output file
		createOutput();
		
	}
	
	
	/**
	 * uploadInput uploads and stores bed files lines from the input file.
	 */
	private void uploadInput(){
		bedLines = new ArrayList<BedLine>();
		int chrmName = -1;
		int start = -1;
		int stop = -1;
		double log2 = -99;
		//Name is based off of the input files name.
		String name = inputFile.substring(inputFile.lastIndexOf('/') + 1, inputFile.indexOf('.', inputFile.lastIndexOf('/'))).toUpperCase() + "-Amplification";
		
		try{
			Scanner s = new Scanner(new File(inputFile));
			String token;
			//Parse each Line, one at a time
			while(s.hasNextLine()){
				 token = s.next();

				 //Skip lines with comments
				 if(token.contains("#")){
					 s.nextLine();
					 continue;
				 }
				 //Parse in the chromosome name
					if(token.equalsIgnoreCase("chrM")){
						chrmName = Definitions.chromosomeM;
					}else if(token.equalsIgnoreCase("chrX")){
						chrmName = Definitions.chromosomeX;
					}else if(token.equalsIgnoreCase("chrY")){
						chrmName = Definitions.chromosomeY;
					}else{
						chrmName = Integer.parseInt(token.substring(token.indexOf("r") + 1));
						
					}
				//Parse in the start location, stop location, and score
				 token = s.next();
				 start = Integer.parseInt(token);
				 token = s.next();
				 stop = Integer.parseInt(token);
				 token = s.next();
				 log2 = Double.parseDouble(token);
			 
				 //Create a bed file line.
				 bedLines.add(new BedLine(chrmName, start, stop, log2, name));
				  
				 //Clear the rest of the line.
				 s.nextLine();
			 
				chrmName = -1;
				start = -1;
				stop = -1;
				log2 = -99;

			}

	
			
		}catch(Exception e){
			U.p("Error populating PeppyFile Peptide List: " + e);
			U.p("*****");
		}//catch
		
		
	}//uploadInput
	
	/**
	 * createOutput calculates the maximum and minimum scores for each line, and assigns a color brightness based on the linear proportion each score has to the maximum.
	 * So a score of 6, with a maximum of 12, would get a (6/12) * max green color = .5 * 255 ~ 128.
	 * So a score of -5, with a minimum of -15, would get a (-5/-15) * max red color = .33 * 255 = 85.  
	 */
	private void createOutput(){
		double maxLog2Score = 0;
		double minLog2Score = 0;

		//Determine maximum and minimum score
		for(BedLine bl: bedLines){
			if(bl.getLog2() > maxLog2Score){
				maxLog2Score = bl.getLog2();
			}
			if(bl.getLog2() < minLog2Score){
				minLog2Score = bl.getLog2();
			}
		}

		//Set each lines color based on its proportion to the maximum/minimum score.
		for(BedLine bl: bedLines){
			if(bl.getLog2() > 0){
				bl.setColor(0, (int)((bl.getLog2() / maxLog2Score) * 255), 0);
			}else{
				bl.setColor((int)((bl.getLog2() / minLog2Score) * 255), 0, 0);
			}
		}

		//Writing to HDD related variables
		StringBuffer sb;
		BufferedWriter out;


		try {
			sb = new StringBuffer();
			out = new BufferedWriter(new FileWriter(outputDir + inputFile.substring(inputFile.lastIndexOf('/')) + ".bed"));
			
			//This information is necessary for colors to be displayed in the UCSC browser
			sb.append("track itemRgb='On'" + "\n");
			//Write each bed line to the string buffer
			for(BedLine bl: bedLines){
				sb.append(bl.toString() + "\n");
			}
			
			
			//Write the string buffer to the file.
			out.write(sb.toString());
			//Flush and close the stream to ensure all data was written to disk.
			out.flush();
			out.close();
			
			
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}//Catch
		
		
	}//createOutput



}//RegionBedFileTool
