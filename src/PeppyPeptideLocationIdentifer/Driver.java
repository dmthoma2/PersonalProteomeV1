package PeppyPeptideLocationIdentifer;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.ArrayList;

import MiscTools.GencodePPComparisonTool;
import PersonalProteome.Annotation;
import PersonalProteome.Properties;
import PersonalProteome.U;

public class Driver {

	
	//Look at maximum memory usage 
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	/**
	 * 
	 * @author David "Corvette" Thomas
	 * @param args The first argument should be the file location of the properties file.
	 */
	public static void main(String[] args) {
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();

		
		/*  hello! */
		printGreeting();
		
		String args0 = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v11.annotation.gtf";
//		String args0 = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v11.chrm12.gtf";
//		String args0 = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v11.chrm1.gtf";
//		String args0 = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v11.chrm15.gtf";
		String args1 = "/Users/davidthomas/Peppy/synthesis/chromosome/hg19/";
//		String args2 = "/Users/davidthomas/Peppy/synthesis/PeptideLocationFinder/GM12878_SDS_Cyto_dta_1330973336733_report.txt";
		String args2 = "/Users/davidthomas/Peppy/synthesis/PeptideLocationFinder/mutations/";
//		String args2 = "/Users/davidthomas/Peppy/synthesis/PeptideLocationFinder/";
		String args3 = "/Users/davidthomas/Peppy/synthesis/PeptideLocationFinder/output/";
//		String args4 = "results.txt";
		
		File jobsDir = new File(args2);
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		if (potentialJobsFiles != null) {
			for (int i = 0; i < potentialJobsFiles.length; i++) {
				if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".txt")) {
					jobFiles.add(potentialJobsFiles[i]);
				}	
			}
		}

		
		/*setup and do some work!*/
		PeptideLocationIdentifer locFinder = new PeptideLocationIdentifer(args0, args1, args3);

		locFinder.uploadProteome();
		
		for(File f: jobFiles){
			locFinder.identifyLocations(f.getAbsolutePath());
		}
		/* i'm finished! */
		printFarewell();
	}
		
	public static void printGreeting() {		
		U.p("Peptide Location Identification at the speed of a turbo-Vette!");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}
}
