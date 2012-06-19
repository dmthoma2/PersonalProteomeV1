package SpectraFileSeperator;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

import BedFileValidation.BedFileValidation;
import PersonalProteome.Definitions;
import PersonalProteome.U;

public class SpectraFileSeperatorDriver {
	//Look at maximum memory usage 
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	
	public static void main(String[] args){
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();

		/*  hello! */
		printGreeting();
		
		String args0 = "/Users/davidthomas/Peppy/synthesis/SpectraFileSeperator/";
		String args1 = "/Users/davidthomas/Peppy/synthesis/SpectraFileSeperator/output/";
		
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
		
		
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
		String startTime = "";
		String endTime = "";
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		

		
		for(int j = 0; j < jobFiles.size(); j++){
			U.p("Working on file " + (j + 1) + " of " + jobFiles.size());
			SpectraFileSeperator sfs = new SpectraFileSeperator(jobFiles.get(j).getAbsolutePath(), args1, cal);
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
