package PersonalProteome.Gene;

import PersonalProteome.Definitions;

/**
 * Seleno represents a selenocysteine.  This amino acid is encoded by UGA (typically a stop codon) and requires special consideration when encoding codons.
 * It contains simple start and stops, as well as the Selenocysteine Sequence from the 'Definitions' folder.  It extends GeneticObject.  <GeneticObject>
 * @author David "Corvette" Thomas
 *
 */

public class Seleno extends GeneticObject {


		

		private String sequence = Definitions.SELENO_SEQUENCE;
		
		/**
		 * Seleno represents a selenocysteine.  This amino acid is encoded by UGA (typically a stop codon) and requires special consideration when encoding codons.
		 */
		public Seleno(int startLocation, int stopLocation){
			super(startLocation, stopLocation);
		}

		/**
		 * @return the sequence
		 */
		public String getSequence() {
			return sequence;
		}

	

}
