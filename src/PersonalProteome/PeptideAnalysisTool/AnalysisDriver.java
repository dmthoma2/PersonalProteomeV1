package PersonalProteome.PeptideAnalysisTool;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;



import PersonalProteome.U;

/**
 * Runs the Peptide Analysis Tool.
 * @author David "Corvette" Thomas
 */
public class AnalysisDriver {

	
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
		
		if(args.length < 1){
			U.p("The properties file location must be the first parameter.");
			return;
		}
		/*Load properties*/

		AnalysisProperties.loadProperties(new File(args[0]));
		
		/*setup and do some work!*/
		if(AnalysisProperties.checkSplice){
			U.p("Checking for splice junctions.");
			PeptideAnalysisTool pepResults = new PeptideAnalysisTool(AnalysisProperties.annotationFile, AnalysisProperties.peptideListFile, AnalysisProperties.outputDir, AnalysisProperties.outputFileName, AnalysisProperties.chrmDir);
			pepResults.checkPeptidesSpanning();
		}else{
			U.p("Checking for peptide locations.");
			PeptideAnalysisTool pepResults = new PeptideAnalysisTool(AnalysisProperties.annotationFile, AnalysisProperties.peptideListFile, AnalysisProperties.outputDir);
			pepResults.checkPeptideLocation();
		}
		
		
		/* i'm finished! */
		printFarewell();
	}
		
	public static void printGreeting() {		
		U.p("Peptide Analysis at the speed of a turbo-Vette!");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}
}



