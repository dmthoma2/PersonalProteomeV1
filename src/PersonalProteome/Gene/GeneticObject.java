package PersonalProteome.Gene;

/**
 * Genetic object represents the basic information each genetic object has.  
 * Each object contains a start and a stop, and can be sorted based on them.
 * It implements Comparable<Genetic Object>.
 * 
 * @author David "Corvette" Thomas
 *
 */
public class GeneticObject implements Comparable<GeneticObject>{

	
	private int start;
	private int stop;
	
	public GeneticObject(){
	}
	
	/**
	 * Constructor for a genetic object.  It contains the minimal information for a genetic object.  It is not enough information for more complex objects
	 * such as 'Transcript', 'CDS', 'RegionOfInterest', etc.  This constructor should only be called by those objects directly.
	 * @param 
	 */
	public GeneticObject(int start, int stop){
		this();
		this.start = start;
		this.stop = stop;
	}

	/**
	 * toString simply returns the start and stop of this object.
	 */
	public String toString(){
		return "Start: " + getStart() + " Stop: " + getStop();
	}
	
	/**
	 * @return The start location of this object
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @return The stop location of this object.
	 */
	public int getStop() {
		return stop;
	}

	/**
	 * @return the length of this object, based on the start and stop locations
	 */
	public int getLength(){
		return stop - start;
	}
	/**
	 * Compares this object to g.  compareTo Sorts on the start, and if they are equal on the stop.
	 * Returns -1 if this object is less then g, -1 if it is greater, and 0 if they have the same starts and stops.
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(GeneticObject g) {
		if(getStart() < g.getStart()){
			return -1;
		}
		if(getStart() > g.getStart()){
			return 1;
		}
		if(getStart() == g.getStart()){
			if(getStop() < g.getStop()){
				return -1;
			}
			if(getStop() > g.getStop()){
				return 1;
			}
		}

		return 0;
	}

}
