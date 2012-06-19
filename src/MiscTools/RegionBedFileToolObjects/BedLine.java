package MiscTools.RegionBedFileToolObjects;

import PersonalProteome.Definitions;

/**
 * BedLine represents a single bed file line with a color
 * @author David "Corvette" Thomas
 *
 */
public class BedLine {
	int chrmName = -1;
	int start = -1;
	int stop = -1;
	double log2 = -99;
	String name;
	int rValue = 0;
	int gValue = 0;
	int bValue = 0;
	int score = 1000;
	
	/**
	 * 
	 * @param chrmName
	 * @param start
	 * @param stop
	 * @param log2 The amplification score associated with this bed line.  This is specific to RegionBedFileTool.
	 * @param name
	 */
	public BedLine(int chrmName, int start, int stop, double log2, String name) {
		super();
		this.chrmName = chrmName;
		this.start = start;
		this.stop = stop;
		this.log2 = log2;
		this.name = name;
	}

	
	

	
	public String toString(){
		
		String tab = "\t";
		return Definitions.convertChrmNumToString(chrmName) + tab + start + tab + stop + tab + name + tab + score + tab + "+" + tab + start + tab + stop + tab + ( rValue + "," + gValue + "," + bValue);
		
	}
	
	public void setScore(int score){
		this.score =score;
	}

	public void setColor(int r, int g, int b){
		rValue = r;
		gValue = g;
		bValue = b;
	}
	
	/**
	 * @return the chrmName
	 */
	public int getChrmName() {
		return chrmName;
	}


	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}


	/**
	 * @return the stop
	 */
	public int getStop() {
		return stop;
	}


	/**
	 * @return the log2
	 */
	public double getLog2() {
		return log2;
	}


	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	
	
}
