package MiscTools.GeneIdentifer;

public class Region {

	int chromosome = -1;
	int start = -1;
	int stop = -1;
	
	
	
	/**
	 * @param chromosome
	 * @param start
	 * @param stop
	 */
	public Region(int chromosome, int start, int stop) {
		super();
		this.chromosome = chromosome;
		this.start = start;
		this.stop = stop;
	}
	
	public String toString(){
		return chromosome + "\t" + start + "\t" + stop + "\t";
	}
	/**
	 * @return the chromosome
	 */
	public int getChromosome() {
		return chromosome;
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

}
