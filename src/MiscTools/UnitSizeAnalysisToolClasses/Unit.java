package MiscTools.UnitSizeAnalysisToolClasses;

import PersonalProteome.Gene.GeneticObject;


/**
 * Unit is a generic object that stores a single annotation entry unit.  Examples: CDS, GENE, TRANSCRIPT, etc
 * It only stores information about a units location and type.
 * @author David "Corvette" Thomas
 *
 */
public class Unit extends GeneticObject{
	private int featureType;
	private int length;
	
		public Unit(int startLoc, int stopLoc, int featureType){
			super(startLoc, stopLoc);
			this.featureType = featureType;
			this.length = stopLoc - startLoc + 1;
		}
		/**
		 * @return the featureType
		 */
		public int getFeatureType() {
			return featureType;
		}

		/**
		 * @return the length
		 */
		public int getLength() {
			return length;
		}
}
