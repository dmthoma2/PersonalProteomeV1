package Peppy;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;

import Reports.HTMLReporter;
import Reports.TextReporter;
import Utilities.U;

/**
 * Peppy
 * A very stripped down Java version of the MS/MS to genome mapping of Morgan Gidding's GFS.
 * Designed with the following goals:
 * 1) More simple code to promote open source development
 * 2) Takes advantage of Java's popularity over Objective-C
 * 3) better multi-threading
 * @author Brian Risk
 *
 */
public class Peppy {
	
	//So that we may report the total amount of peptides found
	static int peptideTally = 0;
	
	public static void main(String [] args) {
		printGreeting();
		init(args);
		runJobs(args);
//		exportPeptideList();
//		checkSequences();
		U.p("done");
	}
	
	public static void checkSequences() {
		U.p("making sure we have the correct start and stop locations");
		Properties.isSequenceFileDNA = true;
		Properties.numberOfMissedCleavages = 0;
		ArrayList<Sequence> sequences = Sequence.loadSequences(new File("/Users/risk2/PeppyOverflow/sequences HG19/chr12.fa"));
		Sequence sequence = sequences.get(0);
		ArrayList<Peptide> peptides;

		peptides = sequence.extractMorePeptides(false);
		//continually extract peptides from the sequence until there aren't anymore
		while (peptides != null) {
			for (Peptide peptide: peptides) {
				if (peptide.equals("STDYGIFQINSR") || peptide.equals("NMQDMVEDYR")) {
					U.p (peptide.getAcidSequenceString() + ": " + peptide.getStartIndex() + "\t" + peptide.getStopIndex());
				}
			}
			peptides.clear();
			System.gc();
			peptides = sequence.extractMorePeptides(false);
		}
		U.p("done");
		
	}
	
	public static void runPeppy(String [] args) {
		U.startStopwatch();
		peptideTally = 0;
		
		//create new report directory
		File reportDir = new File(Properties.reportDirectory, Properties.spectraDirectoryOrFile.getName() + "_" + System.currentTimeMillis());
		
		//save our properties
		Properties.generatePropertiesFile(reportDir);
		
		//Load our spectra
		U.p("loading spectra...");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectra();
		U.p("loaded " +spectra.size() + " spectra.");
		
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		//initialize our ArrayList of matches
		ArrayList<Match> matches = new ArrayList<Match>();
		
		if (Properties.useSpliceVariants) {
			//gets the first nucleotide sequence in the first sequence file
			Sequence sequenceFile = sequences.get(0);
			RNA_Sequence rna = new RNA_Sequence(sequenceFile.getNucleotideSequences().get(0), Properties.sequenceRegionStart, Properties.sequenceRegionStop);
			
			U.p("digesting...");
			RNA_Digestor rnaDigestor = new RNA_Digestor(rna);
			ArrayList<Peptide> peptides  = rnaDigestor.getPeptides();
			U.p("peptide tally: " + peptides.size());
			
			U.p("getting matches...");
			//TODO get rid of this getMatches function when this is overhauled
			matches = getMatchesWithPeptides(peptides, spectra);
			
		} else {	
			if (Properties.useSequenceRegion) {
				U.p("digesting on part of sequence");
				ArrayList<Sequence> oneSequenceList = new ArrayList<Sequence>();
				oneSequenceList.add(sequences.get(0));
				sequences = oneSequenceList;
			}
			matches = getMatches(sequences, spectra);
			U.p("peptide tally: " + peptideTally);
		}
		
		
		
		U.p("creating text reports");
		TextReporter textReport = new TextReporter(matches, spectra, sequences, reportDir);
		textReport.generateFullReport();
		
		
		
		if (Properties.createHTMLReport) {
			U.p("creating HTML reports");
			HTMLReporter report = new HTMLReporter(matches, spectra, sequences, reportDir);
			report.generateFullReport();
		}	
		
		U.p();
		U.stopStopwatch();
	}
	
	
	public static void runJobs(String [] args) {
		File jobsDir = new File("jobs");
		File[] potentialJobsFiles = jobsDir.listFiles();
		ArrayList<File> jobFiles = new ArrayList<File>();
		if (potentialJobsFiles != null) {
			for (int i = 0; i < potentialJobsFiles.length; i++) {
				if (potentialJobsFiles[i].getName().toLowerCase().endsWith(".txt")) {
					jobFiles.add(potentialJobsFiles[i]);
				}	
			}
		}
		if (jobFiles.size() == 0) {
			U.p("no jobs in jobs folder.  running according to main properties file");
			runPeppy(null);
		} else {
			U.p("running " + jobFiles.size() + " jobs");
			for (int i = 0; i < jobFiles.size(); i++) {
				U.p("running job " + (i + 1) + "; " + jobFiles.get(i).getName());
				init(args);
				init(jobFiles.get(i).getAbsolutePath());
				runPeppy(null);
			}
		}
	}


	public static void init(String propertiesFile) {
		System.setProperty("java.awt.headless", "true"); 
		Properties.loadProperties(propertiesFile);
		AminoAcids.init();
	}
	
	
	public static void init(String [] args) {
		if (args.length == 0) {
			init("properties.txt");
		} else {
			init(args[0]);
			U.p("Initializing with properties file: " + args[0]);
		}
	}
	
	
	/**
	 * Assumes that we are doing the normal, forwards digestion of our sequences
	 * @param sequences
	 * @param spectra
	 * @return
	 */
	public static ArrayList<Match> getMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra) {
		return getMatches(sequences, spectra, false);
	}
	
	public static ArrayList<Match> getReverseMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra) {
		return getMatches(sequences, spectra, true);
	}
	
	
	/**
	 * 
	 * @param sequences our list of sequences where we will be getting our peptides
	 * @param spectra our list of spectra
	 * @param isReverse if we are doing a normal, forwards search or if this is a null, reverse search
	 * @return
	 */
	private static ArrayList<Match> getMatches(ArrayList<Sequence> sequences, ArrayList<Spectrum> spectra, boolean isReverse) {
		//TODO: turn these into preferences
		int numberOfSpectraPerSegment = 40000;
		
		int spectraStart = 0;
		int spectraStop = numberOfSpectraPerSegment;
		if (spectraStop > spectra.size()) spectraStop = spectra.size();
		ArrayList<Match> matches = new ArrayList<Match>();
		while (true) {
			ArrayList<Match> segmentMatches = new ArrayList<Match>();
			
			//grab our segment of spectra
			ArrayList<Spectrum> spectraSegment = new ArrayList<Spectrum>();
			for (int i = spectraStart; i < spectraStop; i++) {
				spectraSegment.add(spectra.get(i));
			}
			
			if (Properties.sequenceFilesContainMultipleSequences) {
				U.p("Processing all sequences at once");
				ArrayList<Peptide> peptides = new ArrayList<Peptide>();
				Sequence sequence;
				double percentComplete;
				NumberFormat nfPercent = NumberFormat.getPercentInstance();
				nfPercent.setMaximumFractionDigits(2);
				for (int i = 0; i < sequences.size(); i++) {
					sequence = sequences.get(i);
					peptides.addAll(sequence.extractAllPeptides(isReverse));
					
					//there may be a lot of peptides, so this segments the task
					if (peptides.size() > 10000000) {
						percentComplete = (double) i / sequences.size();
						U.p(nfPercent.format(percentComplete) + " percent complete"); 
						peptideTally += peptides.size();
						ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectraSegment)).getMatches();
						//Possible to add only matches with a decent e value
						evaluateMatches(newMatches, segmentMatches);
						peptides.clear();
						System.gc();
						removeDuplicateMatches(segmentMatches);
					}
				}
				U.p("evaulating for " + peptides.size() + " peptide subset");
				peptideTally += peptides.size();
				ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectraSegment)).getMatches();
				//Possible to add only matches with a decent e value
				evaluateMatches(newMatches, segmentMatches);
				peptides.clear();
				System.gc();
				removeDuplicateMatches(segmentMatches);

			} else {
				for (Sequence sequence: sequences) {
					U.p("working on sequence: " +sequence.getSequenceFile().getName());
					
					ArrayList<Peptide> peptides;
					if (Properties.isSequenceFileDNA) {
						peptides = sequence.extractMorePeptides(isReverse);
						//continually extract peptides from the sequence until there aren't anymore
						while (peptides != null) {
							peptideTally += peptides.size();
							//This is where the bulk of the processing in long jobs takes
							ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectraSegment)).getMatches();
							//Possible to add only matches with a decent e value
							evaluateMatches(newMatches, segmentMatches);
							//free up the memory of the old peptide arraylist
							peptides.clear();
							System.gc();
							peptides = sequence.extractMorePeptides(isReverse);
						}
						sequence.reset();
						removeDuplicateMatches(segmentMatches);
					} else {
						peptides = ProteinDigestion.getPeptidesFromDatabase(sequence.getSequenceFile(), isReverse);
						peptideTally += peptides.size();
						//This is where the bulk of the processing in long jobs takes
						ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectraSegment)).getMatches();
						//Add only matches with a decent e value
						evaluateMatches(newMatches, segmentMatches);
					}		
				}
			}
			
			
			assignConfidenceValuesToMatches(segmentMatches);
			evaluateMatches(segmentMatches, matches);
			
			//increment our spectrum segment markers
			if (spectraStop == spectra.size()) break;
			spectraStart = spectraStop;
			spectraStop += numberOfSpectraPerSegment;
			if (spectraStop > spectra.size()) spectraStop = spectra.size();
			
		}
		assignRankToMatches(matches);
		assignRepeatedPeptideCount(matches);	
		
		return matches;
	}
	
	/**
	 * Gets matches where a list of peptides is already derived
	 * @param peptides
	 * @param spectra
	 * @param sequence
	 * @return
	 */
	public static ArrayList<Match> getMatchesWithPeptides(ArrayList<Peptide> peptides, ArrayList<Spectrum> spectra) {
		ArrayList<Match> matches = new ArrayList<Match>() ;
		peptideTally += peptides.size();
		
		//This is where the bulk of the processing in long jobs takes
		ArrayList<Match> newMatches = (new ScoringThreadServer(peptides, spectra)).getMatches();
		
		//Add only matches with a decent e value
		evaluateMatches(newMatches, matches);
		
		if (Properties.isSequenceFileDNA) {
			removeDuplicateMatches(matches);
		}
		assignRankToMatches(matches);
		assignRepeatedPeptideCount(matches);
		assignConfidenceValuesToMatches(matches);
		
		return matches;
	}
	
		
	public static void assignRankToMatches(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_SCORE);
		Collections.sort(matches);
		Match match = matches.get(0);
		Match previousMatch = matches.get(0);
		//set for the first
		match.rank = 1;
		int rank = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				rank = 1;
			}
			if (match.getScore() == previousMatch.getScore()) {
				rank = previousMatch.rank;
			}
			match.rank = rank;
			rank++;
			previousMatch = match;
		}
		//Setting Score ratios for those with rank 1
		int i = matches.size() - 1;
		double previousScore = match.getScore();
		for (; i >= 0; i--) {
			match = matches.get(i);
			if (match.rank == 1) {
				match.setScoreRatio(match.getScore() / previousScore);
			} else {
				previousScore = match.getScore();
			}
		}
	}
	
	/**
	 * finds the number of times a certain amino acid is found for each spectrum 
	 * @param matches
	 */
	public static void assignRepeatedPeptideCount(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
		Collections.sort(matches);
		Match match;
		Match previousMatch = matches.get(0);
		int rankCount = 1;
		for (int i = 1; i < matches.size(); i++) {
			//see if these are matches for a different spectrum
			match = matches.get(i);
			if (match.getSpectrum().getId() != previousMatch.getSpectrum().getId()) {
				for (int j = i - rankCount; j < i; j++) {
					matches.get(j).repeatCount = rankCount;
				}
				rankCount = 1;
			} else {
				if (match.getPeptide().equals(previousMatch.getPeptide())) {
					rankCount++;
				} else {
					for (int j = i - rankCount; j < i; j++) {
						matches.get(j).repeatCount = rankCount;
					}
					rankCount = 1;
				}
			}
			previousMatch = match;
		}
	}
	
	public static void removeDuplicateMatches(ArrayList<Match> matches) {
		//first error check
		if (matches.size() == 0) return;
		
		Match.setSortParameter(Match.SORT_BY_SPECTRUM_ID_THEN_PEPTIDE);
		Collections.sort(matches);
		int numberOfMatches = matches.size();
		Match match;
		Match previousMatch = matches.get(0);
		int spectrumID;
		int previousSpectrumID = previousMatch.getSpectrum().getId();
		boolean areEqual;
		for (int i = 1; i < numberOfMatches; i++) {
			match = matches.get(i);
			spectrumID = match.getSpectrum().getId();
			areEqual = false;
			if (match.getPeptide().isForward() && previousMatch.getPeptide().isForward()) {
				if (match.equals(previousMatch) && match.getPeptide().getStartIndex() == previousMatch.getPeptide().getStartIndex() && spectrumID == previousSpectrumID) {
					areEqual = true;
					matches.remove(i);
					i--;
					numberOfMatches--;
				}
			}
			if (!match.getPeptide().isForward() && !previousMatch.getPeptide().isForward()) {
				if (match.equals(previousMatch) && match.getPeptide().getStartIndex() == previousMatch.getPeptide().getStartIndex() && spectrumID == previousSpectrumID) {
					areEqual = true;
					matches.remove(i);
					i--;
					numberOfMatches--;
				}
			}
			if (!areEqual) {
				previousMatch = match;
				previousSpectrumID = spectrumID;
			}
		}
	}
	
	public static void assignConfidenceValuesToMatches(ArrayList<Match> matches) {
		for (Match match: matches) {
			if (match.calculateEValue() < match.calculateIMP()) {
				match.setEValue(Double.MAX_VALUE);
			}
		}
	}
	
	/**
	 * Evaluates newMatches, adds appropriate ones to matches
	 * @param newMatches
	 * @param matches
	 */
	public static void evaluateMatches(ArrayList<Match> newMatches, ArrayList<Match> matches) {	
		for (Match match: newMatches) {
			//The match E value should be less than our cutoff
			if (match.getEValue() <= Properties.eValueCutOff) {
				//the match IMP value should always be less than its E value
				if (match.calculateIMP() < match.getEValue()) {
					matches.add(match);
				}
			}
		}
	}

	/**
	 * saves a file of peptide information
	 * @param peptides
	 */
	@SuppressWarnings("unused")
	private static void exportPeptideList() {
		U.p("exporting peptide list");
		//Get references to our sequence files -- no nucleotide data is loaded at this point
		ArrayList<Sequence> sequences = Sequence.loadSequences(Properties.sequenceDirectoryOrFile);
		
		try {
			//loop through each sequence in the sequences ArrayList
			File peptidesFolder = new File ("peptides");
			peptidesFolder.mkdir();
			for (Sequence sequence: sequences) {
				U.p("working on sequence " + sequence.getSequenceFile().getName());
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(peptidesFolder, sequence.getSequenceFile().getName()))));
				
				ArrayList<Peptide> peptides = null;
				if (Properties.isSequenceFileDNA) {
					peptides = sequence.extractMorePeptides(false);
					while (peptides != null) {
						for (Peptide peptide: peptides) {
							pw.println(peptide);
						}
						peptides = sequence.extractMorePeptides(false);
					}
				} else {
					peptides = ProteinDigestion.getPeptidesFromDatabase(sequence.getSequenceFile());
					for (Peptide peptide: peptides) {
						pw.println(peptide);
					}
				}
				
				peptides = null;
				System.gc();
				pw.flush();
				pw.close();
				U.p(sequence.getSequenceFile().getName() + " digested.");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	protected static void printGreeting() {
		U.p("Welcome to Peppy");
		U.p("Proteogenomic mapping software.");
		U.p("Developed 2010 by the Giddings Lab");
		U.p();
		

//		U.p("Peppy Copyright (C) 2011 Brian Risk");
//		U.p("This program comes with ABSOLUTELY NO WARRANTY;");

	}
	

}
