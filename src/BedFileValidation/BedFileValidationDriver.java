package BedFileValidation;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;


import PersonalProteome.Definitions;
import PersonalProteome.U;

public class BedFileValidationDriver {

	//Look at maximum memory usage 
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	
	public static void main(String[] args){
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();

		/*  hello! */
		printGreeting();
		
		//bedFile
//		String args0 = "/Users/davidthomas/Peppy/synthesis/BedFileValidation/BEDFILE_GM12878_SDS_Cyto_dta_1331232525773_report.txt.bed";
//		String args0 = "/Users/davidthomas/Peppy/synthesis/BedFileValidation/BEDFILE_GM12878_WCL_dta_1331305355847_report.txt.bed";
		String args0 = "/Users/davidthomas/Peppy/synthesis/BedFileValidation/";
		//chrm Dir
	
//		String args1 = "/Users/davidthomas/Peppy/synthesis/chromosome/WHIM2_xeno_and_germline/";
//		String args1 = "/Users/davidthomas/Peppy/synthesis/chromosome/WHIM16_germline/";
		String args1 =  "/Users/davidthomas/Peppy/synthesis/chromosome/hg19/";
//		String args1 =  "/Users/davidthomas/Peppy/synthesis/chromosome/hg18/";
		//Output Dir
		String args2 = "/Users/davidthomas/Peppy/synthesis/BedFileValidation/output/";
		
		
		/*Brian Risk*/
		File jobsDir = new File(args0);
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		if (potentialJobsFiles != null) {
			for (int i = 0; i < potentialJobsFiles.length; i++) {
				if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".bed")) {
					jobFiles.add(potentialJobsFiles[i]);
				}	
			}
		}

		
		
		//Run a validation based on each input file
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
		String startTime = "";
		String endTime = "";
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		for(int j = 0; j < jobFiles.size(); j++){
			U.p("Working on file " + (j + 1) + " of " + jobFiles.size());
			BedFileValidation bfv = new BedFileValidation(jobFiles.get(j).getAbsolutePath(), args1, args2, cal );
			
			
			/*Do some work!*/
			bfv.validate();
		}
		
		//Get the ending time
		cal = Calendar.getInstance();
		endTime = sdf.format(cal.getTime());
		U.p("Finishing up: " + endTime);
		
		
		/* i'm finished! */
		printFarewell();
		
	}
	
	public static void printGreeting() {		
		U.p("Bed file validation at the speed of a turbo-Vette!");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}
}
