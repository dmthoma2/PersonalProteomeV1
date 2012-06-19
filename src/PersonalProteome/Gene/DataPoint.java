package PersonalProteome.Gene;

/**
 * DataPoint represents a single point of data comprised of: A Name, A Amplification Score, and a Variant Count.
 * DataPoint can easily be used in an list to compile information about various patterns in variance and amplification score.  
 * @author David "Corvette" Thomas
 *
 */
public class DataPoint implements Comparable<DataPoint> {

	private double ampScore;
	private int variantCount;
	private String name;
	/**
	 * DataPoint represents a single point of data comprised of: A Name, A Amplification Score, and a Variant Count.
	 * DataPoint can easily be used in an list to compile information about various patterns in variance and amplification score
	 * @param name The name of this DataPoint, usually the gene or transcript to which it belongs.
	 * @param ampScore The amplification score of this point.
	 * @param variantCount The number of variations for this data point.
	 */
	public DataPoint(String name, double ampScore, int variantCount){
		this.name = name;
		this.ampScore = ampScore;
		this.variantCount = variantCount;
	}

	/**
	 * The objects are compared based on amplification score. If they are equal then variant count is used.  The objects are considered to be equal if they have the same amplification score and variant count.
	 * @param dp A DataPoint to compare with this object
	 * 
	 */
	public int compareTo(DataPoint dp) {
		if(this.getAmpScore()  < dp.getAmpScore()){
			return -1;
		}
		if(this.getAmpScore() > dp.getAmpScore()){
			return 1;
		}
		if(this.getAmpScore() == dp.getAmpScore()){
			if(this.getVariantCount() < dp.getVariantCount()){
				return -1;
			}
			if(this.getVariantCount() > dp.getVariantCount()){
				return 1;
			}
		}
		return 0;
	}


	
	
	/**
	 * @return A tab delimited string of this objects name, amplification score, and variant count.
	 */
	public String toString(){
		return name + "\t" + getAmpScore() + "\t" + getVariantCount();
	}

	/**
	 * @return the ampScore
	 */
	public double getAmpScore() {
		return ampScore;
	}

	/**
	 * @return the variantCount
	 */
	public int getVariantCount() {
		return variantCount;
	}
}
