package GENCODESubmissionTool;


import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

import PersonalProteome.Definitions;
import PersonalProteome.U;
/**
 * Driver simply starts up and runs GENCODESubmissionTool.  
 * @author David "Corvette" Thomas
 */
public class Driver {
	
	//Look at maximum memory usage 
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();

		String args0 = "/Users/davidthomas/Peppy/GENCODE/input/peppyResults";
		String args1 = "/Users/davidthomas/Peppy/GENCODE/input/FDR-20.0-500.0.txt";
		String args2 = "/Users/davidthomas/Peppy/GENCODE/output";
		int maxMatchRank = 1;
		
		
		/*  hello! */
		printGreeting();
		
		/*Brian Risk*/
		File jobsDir = new File(args0);
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		if (potentialJobsFiles != null) {
			for (int i = 0; i < potentialJobsFiles.length; i++) {
				if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".txt")) {
					jobFiles.add(potentialJobsFiles[i]);
				}	
			}
		}
		
		//Complete a run based on each input file
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
		String startTime = "";
		GENCODESubmissionTool gst;
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		if(jobFiles.size() == 0){
			U.p("No Jobs Uploaded!");
			printFarewell();
			return;
		}
		gst = new GENCODESubmissionTool(jobFiles.get(0).getAbsolutePath(), args1, args2, maxMatchRank, cal);
		for(int j = 0; j < jobFiles.size(); j++){
			U.p("Working on file " + (j + 1) + " of " + jobFiles.size());
			gst.setInputFile(jobFiles.get(j).getAbsolutePath());
			gst.Parse();
		}
		
		gst.Run();
		
		
		/* i'm finished! */
		printFarewell();
	}
	
	public static void printGreeting() {		
		U.p("GENCODE Submission at the speed of a turbo-Vette!");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}

}//Driver
