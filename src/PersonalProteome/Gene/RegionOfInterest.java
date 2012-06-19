package PersonalProteome.Gene;


/**
 * RegionOfInterest is a object representing just a single region of interest.  This region contains information about its location, what chromsome it is on, and what score it has
 * associated with it.  It extends GeneticObject.  <GeneticObject>
 * 
 * @author Originally Brian Risk; Modified by David "Corvette" Thomas
 */
public class RegionOfInterest extends GeneticObject{
	
	//
	private int chromosome;
	private double score;
	
	/**
	 * A region of interest is simply defined by its chromosome, its start and stop, and its score.  
	 * 
	 */
	public RegionOfInterest(int chromosome, int start, int stop, double score) {
		super(start, stop);
		this.chromosome = chromosome;
		this.score = score;
	}

	/**
	 * @return The chromosome this region is in.  The returned value is from 'Definitions'.
	 */
	public int getChromosome() {
		return chromosome;
	}

	/**
	 * @return The score of this region.
	 */
	public double getScore() {
		return score;
	}



}
