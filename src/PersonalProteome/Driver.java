package PersonalProteome;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;

import PersonalProteome.U;

/**
 * Runs the personal proteome.
 * @author David "Corvette" Thomas
 */
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
		
		/*Load properties*/
		Properties.loadProperties(new File(args[0]));
		
		/*setup*/
		Annotation a = new Annotation(Properties.annotationFile, Properties.refGenomeDir, Properties.genomeDir, Properties.outputDir, Properties.regionsOfInterestFile,
										Properties.outputFileName, Properties.statsOutputFile, Properties.dataPointsFileName);
		

		
		/*Do some work!*/
		if(Properties.proteomeLite){
			a = new Annotation(Properties.annotationFile, Properties.genomeDir, Properties.outputDir);
			a.proteomeLite(Properties.debug, Properties.useModifiedStopsAndStarts, Properties.preFilterOutTranWithoutStartCodon);
		}else{
			a.synthesize(Properties.debug, Properties.useModifiedStopsAndStarts, Properties.preFilterOutTranWithoutStartCodon);
		}
		
		
		/* i'm finished! */
		printFarewell();
	}
		
	public static void printGreeting() {		
		U.p("Proteome creation at the speed of a turbo-Vette! V1.0");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}
}
