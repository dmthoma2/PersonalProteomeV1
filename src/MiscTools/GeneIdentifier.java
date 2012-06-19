package MiscTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import MiscTools.GeneIdentifer.Region;
import PersonalProteome.Annotation;
import PersonalProteome.Definitions;
import PersonalProteome.GENCODE_GTF_Line;
import PersonalProteome.U;




/**
 * Takes in a tab delimited file with a chromosome name, start, and stop with 1 line of header information.
 * Outputs a list of genes that occur with the input regions.  
 * @author David "Corvette" Thomas
 *
 */
public class GeneIdentifier {

	//Input
	String annotation = "/Users/davidthomas/Peppy/synthesis/annotation/gencode.v11.annotation.gtf";
	String regions = "/Users/davidthomas/Peppy/synthesis/MiscTools/GeneIdentifierInRegion/WHIM16Regions.txt";
	String output = "/Users/davidthomas/Peppy/synthesis/MiscTools/GeneIdentifierInRegion/output/";
	
	
	//List of regions
	ArrayList<Region> regionslist = new ArrayList<Region>();
	
	//GTF master list
	ArrayList<GENCODE_GTF_Line> GTFlist = new ArrayList<GENCODE_GTF_Line>();
	
	public static void main(String[] args){
		
		GeneIdentifier gi = new GeneIdentifier();
		
		gi.getGenes();
	}
	
	
	public void getGenes(){
		//UploadAnnotation
		Annotation a = new Annotation(true);
		a.populateGTFList(new File(annotation));
		GTFlist = a.getAnnotaitonLines();
		a = null;
		
		//Upload regions
		uploadRegions();
		
		//Create Output
		createOutput();
	}
	
	
	/**
	 * uploadRegions uploads the comma delmitied
	 */
	public void uploadRegions(){
		try {
			Scanner s = new Scanner(new File(regions));
			
			String token = s.nextLine();
			
			while(s.hasNextLine()){
				
				int chromosome = -1;
				int start = -1;
				int stop = -1;
				
				//Chromosome
				token = s.next();
				if(token.equalsIgnoreCase("M")){
					chromosome = Definitions.chromosomeM;
				}else if(token.equalsIgnoreCase("X")){
					chromosome = Definitions.chromosomeX;
				}else if(token.equalsIgnoreCase("Y")){
					chromosome = Definitions.chromosomeY;
				}else{
					chromosome = Integer.parseInt(token);
					
				}
				
				//Start Location
				token = s.next();
				start = Integer.valueOf(token);
				
				//StopLocation
				token = s.next();
				stop = Integer.valueOf(token);
				
				regionslist.add(new Region(chromosome, start, stop));
				
				U.p(regionslist.get(regionslist.size() -1 ).toString());
				
				chromosome = -1;
				start = -1;
				stop = -1;
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		
	}//uploadRegions
	
	/**
	 * Create Output checks to see if a gene is within a output region, and if it is it writes it to the output file.
	 */
	public void createOutput(){
		
		
		StringBuffer sb;
		BufferedWriter out;
		
		try {
			out = new BufferedWriter(new FileWriter(output + "GeneList.txt"));
			sb = new StringBuffer();
			
			
			int count = 0;
			for(int i = 0; i < GTFlist.size(); i++){
				
				if(GTFlist.get(i).getFeatureType() == Definitions.featureTypeGENE){
					if(GTFlist.get(i).getTranscript_Type().equalsIgnoreCase("protein_coding")){
						
						
//						for(int k = 0; k < regionslist.size(); k++){
//							int geneStart = GTFlist.get(i).getStartLocation();
//							int geneStop = GTFlist.get(i).getStopLocation();
//							
//							if(regionslist.get(k).getStart() <= geneStart && regionslist.get(k).getStop() >= geneStop && GTFlist.get(i).getChromosomeName() == regionslist.get(k).getChromosome()){
//								sb.append(GTFlist.get(i).getGene_Name() + "\n");
//								count++;
//								U.p(geneStart + " " + geneStop + " is inside " + regionslist.get(k).getStart() + " " + regionslist.get(k).getStop());
//							}//within a region
//							
//							
//						}//regionsLoop
							
						if(GTFlist.get(i).getChromosomeName() == 8){
							sb.append(GTFlist.get(i).getGene_Name() + "\t" + GTFlist.get(i).getStartLocation() + "\t" + GTFlist.get(i).getStopLocation() +  "\n");
							count++;
						}
							
					}//protein_coding
				}//gene
				
				
			}//for
						
			out.write("Total number of genes matched: " + count + "\n");
			out.write(sb.toString());
			out.flush();
			out.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}//createOutput
	
	
}//GeneIdentifier
