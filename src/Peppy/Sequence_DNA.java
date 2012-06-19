package Peppy;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;



/**
 * A "sequence" is a DNA sequence, typically in a FASTA file
 * This farms out the digestion of the nucleotide sequences to
 * SequenceDigestionTread objects
 * @author Brian Risk
 *
 */
public class Sequence_DNA{

	private ArrayList<Nucleotides> nucleotideSequences = null;
	File sequenceFile = null;
	
	
	public Sequence_DNA(File sequenceFile) {
		this.sequenceFile = sequenceFile;
	}
	
	public Sequence_DNA(String sequenceFileName) {
		this(new File(sequenceFileName));
	}
	
	
	/**
	 * in one FASTA file there may be many sequences.
	 * this method returns a ArrayList containing all of the sequences
	 * @return
	 */
	public ArrayList<Nucleotides> getNucleotideSequences() {
		if (nucleotideSequences == null) {
			nucleotideSequences = new ArrayList<Nucleotides>();
			try {
				BufferedReader br = new BufferedReader(new FileReader(sequenceFile));
				String line = br.readLine();

				String sequenceDescription = null;
				StringBuffer sequence = new StringBuffer();
				while (line != null) {
					
					/* lines that begin with ">" are comments which name the sequence */
					if (line.startsWith(">")) {
						
						/* if sequenceDescription has not been defined, that means it is the first in the file */
						if (sequenceDescription != null) {
							nucleotideSequences.add(new Nucleotides(sequenceDescription, sequence.toString().toUpperCase(), this));
							sequence = new StringBuffer();
						}
						sequenceDescription = line;
						line = br.readLine(); 
						continue;	
					}
					
					/* lines that begin with ";" are comments which should be ignored */
					if (line.startsWith(";")) {
						line = br.readLine(); 
						continue;	
					}
					sequence.append(line);
					line = br.readLine();
				}
				String acidSequence = sequence.toString().toUpperCase();
				nucleotideSequences.add(new Nucleotides(sequenceDescription, acidSequence, this));
				
				//close out our stream.  It's the courteous thing to do!
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return nucleotideSequences;
		} else {
			return nucleotideSequences;
		}
	}//getNucleotideSequences

}//sequence_DNA
