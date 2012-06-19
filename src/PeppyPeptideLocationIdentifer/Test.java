package PeppyPeptideLocationIdentifer;

import Peppy.Sequence_DNA;
import PersonalProteome.U;

public class Test {

	
	public static void main(String[] args){
		
		Test t = new Test();
		t.populateChrmArray("/Users/davidthomas/Peppy/synthesis/chromosome/hg19/");
		
		int chrmNum = 0;
		int startIndex = 236385224 + 0 - 0;
		int stopIndex = 236385229 + 0 - 0;
		U.p("Getting sequence on Chromosome " + (chrmNum + 1) + " with locations: " + startIndex + " " + stopIndex);
		U.p("Sequence " + t.getBaseSequence(chrmNum, startIndex, stopIndex));
		
	}
	
	
	private String getBaseSequence(int chrmNum, int startIndex, int stopIndex){
		String sequence = "";
		
		
		Sequence_DNA seq = new Sequence_DNA(chrmFile[chrmNum]);
		
		sequence = seq.getNucleotideSequences().get(0).getSequence().substring(startIndex - 1, stopIndex);
		return sequence;
	}
	
	
	private String[] chrmFile = new String[25];
	/**
	 * Populates the chromosome array with the file locations of each chromosome.'
	 * @param chrmDir The directory containing the chromosome files.
	 */
	public void populateChrmArray(String chrmDir){
		for(int i = 0; i < 22; i++){
			chrmFile[i] = chrmDir + "chr" + (i + 1) + ".fa";
//			U.p(chrmFile[i]);
		}
		chrmFile[22] = chrmDir + "chrM.fa";
//		U.p(chrmFile[22]);
		chrmFile[23] = chrmDir + "chrX.fa";
//		U.p(chrmFile[23]);
		chrmFile[24] = chrmDir + "chrY.fa";
//		U.p(chrmFile[24]);
	}
}
