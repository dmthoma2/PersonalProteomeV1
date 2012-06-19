package PersonalProteome.PeptideAnalysisTool.SubClasses;


import java.util.ArrayList;
import PersonalProteome.Gene.Transcript;

/**
 * Peptide is a simple object that represents a single peptide from a Peppy results file.  This peptide knows basic information about it self such as its transcript ID, its gene ID, and its 
 * amino acid sequence. Each peptide can determine if it occurs over a transcirpts introns, and if it does how many.
 * 
 * @author David "Corvette" Thomas
 *
 */
public class Peptide {

	//debug
	private static int count = 0;
	
	//
	private String transcriptID;
	private String sequence;
	private String geneID;
	private int spanCount = -1;
	
	private ArrayList<Double> juncs = new ArrayList<Double>();
	
	


	/**
	 * Peptide is a simple object that represents a single peptide from a Peppy results file.  This peptide knows basic information about it self such as its transcript ID, its gene ID, and its 
	 * amino acid sequence. Each peptide can determine if it occurs over a transcirpts introns, and if it does how many.
	 * 
	 * @param transcriptID  The id from the transcript this peptide is located on.
	 * @param geneID The id from the gene this peptide is located on.
	 * @param sequence The amino acid sequence of this peptide.
	 */
	public Peptide(String transcriptID,String geneID, String sequence){
		this.transcriptID = transcriptID;
		this.geneID = geneID;
		this.sequence = sequence;
	}

	/**
	 * @return the transcriptID
	 */
	public String getTranscriptID() {
		return transcriptID;
	}
	/**
	 * @return the geneID
	 */
	public String getGeneID() {
		return geneID;
	}
	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}
	
	/**
	 * @return the length of this Peptide
	 */
	public int getLength(){
		return sequence.length();
	}
	/**
	 * @return the spanCount.  -1 if it is never compared.  0 if it is compared and nothing is found, otherwise returns the number of exons splits this peptide covers.
	 */
	public int getSpanCount() {
		return spanCount;
	}
	
	/**
	 * 
	 * @return The splice junctions this peptide occurs over.
	 */
	public ArrayList<Double> getSpliceJuncs(){
		return juncs;
	}
	


	/**
	 * Determine if on split checks to see if this peptide occurs on a splice junction within a transcript.  It stores which splice junctions it occurs over in spliceLocation.
	 * 
	 * @param t  Transcript to check which splice junction this occurs over.
	 */
	public void determineIfOnSplit(Transcript t){
		//It has been compared so mark it as having no spans.
		spanCount = 0;
		int startLocation = t.getProtein().indexOf(this.sequence);
		int test = t.getProtein().substring(startLocation + 1).indexOf(this.sequence);
		if(test != -1){
//			U.p("Double peptide found!!");
			count++;
//			U.p(count);
		}
		
		if(startLocation == -1){
			return;
		}//sequenceNotFound
		
		int endLocation = startLocation + this.sequence.length();
		for(Double splitLocation: t.getExonSplitLocations()){
			
			//Error check the split locations
			double splitLoc = splitLocation;
			int floor = (int)splitLoc;
			 double ceiling = floor + 1;
		     double  minValue = floor + .99;
		     
		     if(splitLoc < ceiling && splitLoc > minValue ){
		    	 splitLoc = ceiling;
		     }
			
		     
		     //If the peptide occurs on this split location add it to the list.
			if(startLocation < splitLoc && endLocation > splitLoc){
				juncs.add(splitLocation);
				spanCount++;
			}//if
		}//for
		
		if(spanCount > 0){
			this.sequence += "  " + t.getProteinWithSpliceMarkers();
		}
		
	}//determineIfOnSplit
}
