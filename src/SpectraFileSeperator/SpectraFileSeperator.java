package SpectraFileSeperator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Scanner;

import PersonalProteome.U;

/**
 * This is a simple tool that takes in file containing more then one spectra, and separates them.  The format for these files is:
 * ================"FileName"=======================
 * double double
 * double double
 * double double
 * double double
 * 
 * ==============="NextFileName"===================
 * double double
 * double double
 * ...
 * @author David "Corvette" Thomas
 *
 */
public class SpectraFileSeperator {

	private String outputDir;
	private String inputFile;
	private Calendar cal;
	
	
	private SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	
	/**
	 * Once created this object splits up the spectra file and puts each indiviual file into the output.
	 * @param inputFile
	 * @param outputDir
	 */
	public SpectraFileSeperator(String inputFile, String outputDir, Calendar cal){
		this.outputDir = outputDir;
		this.inputFile = inputFile;
		this.cal = cal;
		splitFiles();
	}
	
	
	/**
	 * splitFiles breaks apart larger joined files into smaller ones.  It creates these files automatically as it goes along.
	 */
	private void splitFiles(){
		//Writing to HDD related variables
		StringBuffer sb;
		BufferedWriter out;
		File outFile;
		
		
		//Create a output file specific to the time of this run.
		outFile = new File( outputDir + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");// + inputFile.substring(inputFile.lastIndexOf("/") + 1, inputFile.length() - 4) + "/");
		outFile.mkdir();
		String outputDir = outFile.getAbsolutePath() + "/";
		
		outFile = new File( outputDir + "/" + inputFile.substring(inputFile.lastIndexOf("/") + 1, inputFile.length() - 4) + "/");
		outFile.mkdir();

		outputDir = outFile.getAbsolutePath() + "/";
		U.p("Output dir " + outputDir);
		try {
		
			Scanner s = new Scanner(new File(inputFile));
			
			String token;
			s.next();
			//Loop and create files based on the parsed file name and information
			while(s.hasNextLine()){
		
			
				token = s.next();
				sb = new StringBuffer();
				out = new BufferedWriter(new FileWriter(outputDir  +  token.substring(1, token.length() - 1)));
			
				
				token = s.next();
				token = s.next();
				
				while(!token.contains("=")){
					sb.append(token + " " + s.next() + "\n");
					if(s.hasNext()){
					token = s.next();
					}else{
						break;
					}
				}
				
				out.write(sb.toString());
				out.flush();
				out.close();
				
				if(!s.hasNext()){
					break;
				}
			}
			
		
		//Default exception and action
		}catch (IOException e) {
			e.printStackTrace();
		}//catch
		
	}//SplitFiles
}//SpectraFileSeparater
