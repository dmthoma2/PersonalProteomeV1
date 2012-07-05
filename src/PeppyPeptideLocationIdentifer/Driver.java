package PeppyPeptideLocationIdentifer;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.ArrayList;

import PeppyPeptideLocationIdentifer.Properties;
import MiscTools.GencodePPComparisonTool;
import PersonalProteome.Annotation;
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
		
		if(args.length != 1){
			System.out.println("The properties file should be the only parameter.");
			return;
		}
		Properties.loadProperties(new File(args[0]));
		
		File jobsDir = new File(Properties.peppyFileDirectory);
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
		PeptideLocationIdentifer locFinder = new PeptideLocationIdentifer(Properties.annotationFile, Properties.outputDir);

		locFinder.uploadAnnotation();
		
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
