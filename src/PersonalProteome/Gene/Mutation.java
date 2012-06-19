package PersonalProteome.Gene;


/**
 * 'Mutation' is a representation of a single mutation in the genetic code.  The mutation types may be found in 'Definitions'.
 * Each mutation knows the name of the gene it came from, its location, and its error type.  
 * @author David "Corvette" Thomas
 *
 */
public class Mutation{

	private int errorType;
	private String gene;
	private int location;

	public Mutation(String geneName, int errorType, int location){
		this.errorType = errorType;
		this.gene = geneName;
		this.location = location;
	}
	
	public boolean equals(String geneName, int errorType, int location){
		
		if(this.errorType == errorType && this.gene.equals(geneName) && this.location == location){
			return true;
		}
		
		return false;
	}
	
	/**
	 * @return the errorType
	 */
	public int getErrorType() {
		return errorType;
	}
	/**
	 * @return the gene
	 */
	public String getGene() {
		return gene;
	}
	/**
	 * @return the location
	 */
	public int getLocation() {
		return location;
	}
	

}
