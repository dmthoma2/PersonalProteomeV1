package MiscTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Scanner;

import PersonalProteome.Annotation;
import PersonalProteome.Definitions;
import PersonalProteome.GENCODE_GTF_Line;
import PersonalProteome.U;

/**
 * UnitSizeAnalysisTool is a simple application that reads in a .GTF format annotation file, and it calculates statistics about the size of various parts of protein coding genes.
 * It determines information about genes, transcripts, exons, and CDS entries.  Statistics collected include: Size, Mean, Median, Standard Deviation, Minimum, Maximum.
 * @author David "Corvette" Thomas
 *
 */
public class UnitSizeAnalysisTool {

	//Look at maximum memory usage 
	static MemoryUsage memoryUsage;
	static long maxMemoryUsed = 0;
	
	public static void main(String[] args){
		/* track memory */
		MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
		memoryUsage = mbean.getHeapMemoryUsage();

		
		/*  hello! */
		printGreeting();
		
		String args0a = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v11.annotation.gtf";
//		String args0a = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v7.chrm1.gtf";
		String args1 = "/Users/davidthomas/Peppy/synthesis/UnitSizeAnalysis/output/";
		String args2 = "results.txt";

		
		/*setup and do some work!*/
		UnitSizeAnalysisTool tool = new UnitSizeAnalysisTool(args0a, args1, args2);

		tool.checkUnitSize();
		
		/* i'm finished! */
		printFarewell();
		
		
	}
	
	public static void printGreeting() {		
		U.p("Annotation Sub Unit Size Analysis at the speed of a turbo-Vette!");
		U.p("max available memory: " + (double) memoryUsage.getMax() / (1024 * 1024 * 1024) + " gigabytes");
	}
	
	public static void printFarewell() {
		U.p("Until next time, signing off...");
	}

	//Time
	SimpleDateFormat sdf = new SimpleDateFormat(Definitions.DATE_FORMAT);
	SimpleDateFormat outputDirectoryFormat = new SimpleDateFormat("Mdyyyykm");
	Calendar cal;
	String startTime = "";
	String endTime = "";
	
	//File Location Information
	private String annotationFile;
	private String outputDir; 
	private String outputFileName;
	
	//Size list
	ArrayList<Integer> exonSizes = new ArrayList<Integer>();
	ArrayList<Integer> geneSizes = new ArrayList<Integer>();
	ArrayList<Integer> transcriptSizes = new ArrayList<Integer>();
	ArrayList<Integer> CDSSizes = new ArrayList<Integer>();
	
	//Standard Deviation
	
	
	//Mean/Median
	double exonSizeMean = 0;
	int exonSizeMedian = 0;
	double geneSizeMean = 0;
	int geneSizeMedian = 0;
	double transcriptSizeMean = 0;
	int transcriptSizeMedian = 0;
	double CDSSizeMean = 0;
	int CDSSizeMedian = 0;
	
	//Standard Deviation
	double exonStanDev = 0;
	double geneStanDev = 0;
	double transcriptStanDev = 0;
	double CDSStanDev = 0;
	
	
	//Max//Min
	int exonSizeMax = 0;
	int exonSizeMin = 100000000;
	int geneSizeMax = 0;
	int geneSizeMin = 100000000;
	int transcriptSizeMax = 0;
	int transcriptSizeMin = 100000000;
	int CDSSizeMax = 0;
	int CDSSizeMin = 100000000;
	
	//List to store all of the transcirpts
	ArrayList<GENCODE_GTF_Line> GTFlist;
	public UnitSizeAnalysisTool(String annotationFile, String outputDir, String outputFileName){
		this.annotationFile = annotationFile;
		this.outputDir = outputDir;
		this.outputFileName = outputFileName;
	}
	
	
	/**
	 * checkUnitSize is the driving method in UnitSizeAnalysisTool.  It uploads the annotation, calculates statistics, and creates a output file.
	 */
	public void checkUnitSize(){
		
		//Get the current time
		cal = Calendar.getInstance();
		startTime = sdf.format(cal.getTime());
		U.p("Starting up: " + startTime);
		
		
		U.p("Uploading the annotation.");
		U.startStopwatch();
		Annotation a = new Annotation(true);
		a.populateGTFList(new File(annotationFile));
		this.GTFlist = a.getAnnotaitonLines();
		U.stopStopwatch();
		
		U.p("Annotation Size " + GTFlist.size());
		
		U.p("Calculating stats.");
		U.startStopwatch();
		calculateStats();
		U.stopStopwatch();
		
		U.p("Writing Output files.");
		U.startStopwatch();
		createOutput();
		U.stopStopwatch();
		
	}
	
	
	public void createOutput(){
		StringBuffer sb;
		BufferedWriter out;
		File outFile;
		try {
		sb = new StringBuffer();
		outFile = new File(this.outputDir + "/" + outputDirectoryFormat.format(cal.getTime()) +  "/");
		outFile.mkdir();
		String outputDir = outFile.getAbsolutePath() + "/";
		out = new BufferedWriter(new FileWriter(outputDir + outputFileName));
		
		sb.append("This run was completed at " + endTime + "\n");
		sb.append("The annotation file is " + annotationFile + "\n");
		sb.append("Analysis tool only considers Units that are 'protein_coding' or 'nonsense_mediated_decay'");
		sb.append("\n");
		sb.append("\n");
		sb.append("\n");
		
		sb.append("EXON" + "\n");
		sb.append("Total:              " + exonSizes.size() + "\n");
		sb.append("Minimum:            " + exonSizeMin + "\n");
		sb.append("Maximum:            " + exonSizeMax + "\n");
		sb.append("Mean:               " + exonSizeMean + "\n");
		sb.append("Median:             " + exonSizeMedian + "\n");
		sb.append("Standard Deviation: " + exonStanDev + "\n");
		sb.append("\n");
		sb.append("\n");
		
		sb.append("GENE" + "\n");
		sb.append("Total:              " + geneSizes.size() + "\n");
		sb.append("Minimum:            " + geneSizeMin + "\n");
		sb.append("Maximum:            " + geneSizeMax + "\n");
		sb.append("Mean:               " + geneSizeMean + "\n");
		sb.append("Median:             " + geneSizeMedian + "\n");
		sb.append("Standard Deviation: " + geneStanDev + "\n");
		sb.append("\n");
		sb.append("\n");
		
		sb.append("TRANSCRIPT" + "\n");
		sb.append("Total:              " + transcriptSizes.size() + "\n");
		sb.append("Minimum:            " + transcriptSizeMin + "\n");
		sb.append("Maximum:            " + transcriptSizeMax + "\n");
		sb.append("Mean:               " + transcriptSizeMean + "\n");
		sb.append("Median:             " + transcriptSizeMedian + "\n");
		sb.append("Standard Deviation: " + transcriptStanDev + "\n");
		sb.append("\n");
		sb.append("\n");
		
		sb.append("CDS" + "\n");
		sb.append("Total:              " + CDSSizes.size() + "\n");
		sb.append("Minimum:            " + CDSSizeMin + "\n");
		sb.append("Maximum:            " + CDSSizeMax + "\n");
		sb.append("Mean:               " + CDSSizeMean + "\n");
		sb.append("Median:             " + CDSSizeMedian + "\n");
		sb.append("Standard Deviation: " + CDSStanDev + "\n");
		sb.append("\n");
		sb.append("\n");
		out.write(sb.toString());
		out.flush();
		out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}//create output
	
	public void calculateStats(){

		//OPTIONAL WEEDING OUT 
//		boolean weedOut = false;
//		weedOut = true;
//		int falseCount = 0;
//		if(weedOut){
//			ArrayList<GENCODE_GTF_Line> temp = new ArrayList<GENCODE_GTF_Line>();
//			//Remove transcripts that do not have a start codon.  Same as filter option available to 
//			U.p(GTFlist.size());
//			for(int i = 0; i < GTFlist.size(); i++){
//				boolean add = true;
//				if(GTFlist.get(i).getFeatureType() == Definitions.featureTypeTRANSCRIPT){
//					boolean startFound = false;
//					int k = i + 1;
//					while(GTFlist.get(k).getFeatureType() != Definitions.featureTypeTRANSCRIPT){
//						if(GTFlist.get(k).getFeatureType() == Definitions.featureTypeCDS){
//							startFound = true;
//						}
//						k++;
//						if(k >= GTFlist.size()){
//							break;
//						}
//					}
//					if(startFound == false){
//						falseCount++;
//					}
//					add = startFound;
//				}
//				
//				if(add){
//					temp.add(GTFlist.get(i));
//				}
//			}
//			GTFlist = temp;
//			
//			U.p("False Count is: " + falseCount);
////			temp = new ArrayList<GENCODE_GTF_Line>();
//			U.p(GTFlist.size());
//		}
		
		//Populate the data arrays, and figure out the max and minimum size of each region
		for(GENCODE_GTF_Line g: GTFlist){
			
			//Min and Max
			if(g.getFeatureType() == Definitions.featureTypeEXON){
				exonSizes.add(g.getStopLocation() - g.getStartLocation() + 1);
				int size = exonSizes.get(exonSizes.size() - 1);
				if(size > exonSizeMax){
					exonSizeMax = size;
				}
				if(size < exonSizeMin){
					exonSizeMin = size;
				}
				exonSizeMean += size;
				continue;
			}

			if(g.getFeatureType() == Definitions.featureTypeCDS){
				CDSSizes.add(g.getStopLocation() - g.getStartLocation() + 1);
				int size = CDSSizes.get(CDSSizes.size() - 1);
				if(size > CDSSizeMax){
					CDSSizeMax = size;
				}
				if(size < CDSSizeMin){
					CDSSizeMin = size;
				}
				CDSSizeMean += size;
				continue;
			}
			if(g.getFeatureType() == Definitions.featureTypeTRANSCRIPT){
				transcriptSizes.add(g.getStopLocation() - g.getStartLocation() + 1);
				int size = transcriptSizes.get(transcriptSizes.size() - 1);
				if(size > transcriptSizeMax){
					transcriptSizeMax = size;
				}
				if(size < transcriptSizeMin){
					transcriptSizeMin = size;
				}
				transcriptSizeMean += size;
				continue;
			}
			if(g.getFeatureType() == Definitions.featureTypeGENE){
				geneSizes.add(g.getStopLocation() - g.getStartLocation() + 1);
				int size = geneSizes.get(geneSizes.size() - 1);
				if(size > geneSizeMax){
					geneSizeMax = size;
				}
				if(size < geneSizeMin){
					geneSizeMin = size;
				}
				geneSizeMean += size;
				continue;
			}
			

		}//for
		
		//Mean and Median
		exonSizeMean /= exonSizes.size();
		CDSSizeMean /= CDSSizes.size();
		transcriptSizeMean /= transcriptSizes.size();
		geneSizeMean /= geneSizes.size();
		
		Collections.sort(exonSizes);
		Collections.sort(CDSSizes);
		Collections.sort(transcriptSizes);
		Collections.sort(geneSizes);
		exonSizeMedian = exonSizes.get(exonSizes.size() /2);
		CDSSizeMedian = CDSSizes.get(CDSSizes.size() /2);
		transcriptSizeMedian = transcriptSizes.get(transcriptSizes.size() /2);
		geneSizeMedian = geneSizes.get(geneSizes.size() /2);
		
		
		//Standard Deviation
		exonStanDev = calculateStandardDeviation(exonSizes, exonSizeMean);
		geneStanDev = calculateStandardDeviation(geneSizes, geneSizeMean);
		transcriptStanDev = calculateStandardDeviation(transcriptSizes, transcriptSizeMean);
		CDSStanDev = calculateStandardDeviation(CDSSizes, CDSSizeMean);
		
	}//calculateStats
	
	
	/**
	 * Calculates the standard deviation of a data set, based on the data sets list of size.  ***THIS METHOD IS UNTESTED, so confirm it works next RUN if used****
	 * @param sizes
	 * @param mean
	 * @return the standard deviation
	 */
	private double calculateStandardDeviation(ArrayList<Integer> sizes){
		double avg = 0;
		
		for(Integer i: sizes){
			avg += i;
		}
		 
		avg /= sizes.size();
		
		return calculateStandardDeviation(sizes, avg);
	}
	/**
	 * Calculates the standard deviation of a data set, based on the data sets list of sizes, and a mean
	 * @param sizes
	 * @param mean
	 * @return the standard deviation
	 */
	private double calculateStandardDeviation(ArrayList<Integer> sizes, double mean){
		
		double standDev = 0;
		
		ArrayList<Double> deviations = new ArrayList<Double>((int) (sizes.size() * 2));
		for(Integer i: sizes){
			deviations.add(i - mean);
		}
		
		ArrayList<Double> temp = new ArrayList<Double>((int) (sizes.size() * 2));
		for(int k = 0; k < deviations.size(); k++){
			Double workingOn = deviations.get(k);
			workingOn = Math.pow(workingOn, 2);
			temp.add(workingOn);
		}
		deviations = temp;
		
		
		double squaresTotal = 0;
		for(Double d: deviations){
			squaresTotal += d;
		}
		
		
		standDev = Math.sqrt((squaresTotal / (deviations.size() - 1)));
		
		
		return standDev;
	}//Standard Deviation


}//UnitSizeAnalysisTool
