/**
 * 
 */
package PersonalProteome.Gene;

/**
 * CDS represents a CDS in a annotation.  A CDS may be within several transcripts, and one or more CDS combine together within a transcript to create a
 * DNA sequence to be encoded into a protein.  It extends GeneticObject.  <GeneticObject>
 * @author David "Corvette" Thomas
 *
 */
public class CDS extends GeneticObject{

	private int genomicPhase;
	/**
	 * A CDS has to know its start and stop location, and is comparable based on its starting value.
	 */
	public CDS( int startLocation, int stopLocation, int genomicPhase){
		super(startLocation, stopLocation);
		
		this.genomicPhase = genomicPhase;

	}

	/**
	 * @return the genomicPhase
	 */
	public int getGenomicPhase(){
		return genomicPhase;
	}
}
