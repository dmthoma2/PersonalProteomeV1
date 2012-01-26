package Peppy;
import java.util.ArrayList;

import Utilities.U;


/**
 * Manages a group of threads that score spectra against the peptide list.
 * each ScoringThread holds one spectrum and a reference to the peptide list.
 * When that thread is done executing it asks ScoringEngine for another
 * spectrum that it may search for.
 * @author Brian Risk
 *
 */
public class ScoringThreadServer {
	
	ArrayList<Spectrum> spectra;
	/*
	 * Having a ArrayList of ArrayLists may look overly complicated, but I am doing this for a reason.
	 * The easier way is to just have a ArrayList of SpectrumPeptideMatch objects and each time 
	 * getNextSpectrum we take the result we're given and "addAll".  The unfortunate thing is
	 * that since getNextSpectrum is synchronized it makes all the other possibly dozens of
	 * threads wait.  Let's not make them wait.  Let's let them do their work as fast as they
	 * can and then merge all the match ArrayLists together at the end.
	 */
	ArrayList<ArrayList<Match>> matches;
	
	ArrayList<Thread> threads;
	
	//this is how we keep track of which Spectrum to give out next
	private int spectrumIndex = 0;
	
	/* why this simply isn't Properties.numberOfThreads is because we may have less spectra than that number */
	private int numberOfThreads;
	
	/**
	 * 
	 * @param peptides
	 * @param spectra
	 * @param matches the ArrayList where we store the best matches
	 */
	public ScoringThreadServer(ArrayList<Peptide> peptides, ArrayList<Spectrum> spectra) {
		this.spectra = spectra;	
		matches = new ArrayList<ArrayList<Match>>(spectra.size());
		
		//here we make sure we don't use more threads than we have spectra
		numberOfThreads = Properties.numberOfThreads;
		if (numberOfThreads > spectra.size()) numberOfThreads = spectra.size();
		threads = new ArrayList<Thread>(numberOfThreads);
		
		//spawn new threads as needed
		for (int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++) {
			ScoringThread scorer = new ScoringThread(getNextSpectrum(), peptides, this);
			Thread thread = new Thread(scorer);
			thread.start();
			threads.add(thread);	
		}
	}
	
	/**
	 * 
	 * @param matchesForOneSpectrum
	 * @return the next Spectrum on the list.  Null if there are no more.
	 */
	public synchronized Spectrum getNextSpectrum(ArrayList<Match> matchesForOneSpectrum) {
		matches.add(matchesForOneSpectrum);
		return getNextSpectrum();
	}
	
	/**
	 * takes no parameters because it is used at the start when there are no matches to incorporate.
	 * @return
	 */
	public synchronized Spectrum getNextSpectrum() {
		Spectrum out =  null;
		if (spectrumIndex < spectra.size()) {
			out = spectra.get(spectrumIndex);
			spectrumIndex++;
		}
		return out;
	}
	
	/**
	 * First we wait for all of  our scoring threads to finish.
	 * "Yo ScoringThreads, I'm happy for you and I'ma let you finish"
	 * @return accumulated 
	 */
	public ArrayList<Match> getMatches() {
		boolean going = true;
		while (going) {
			for (Thread thread: threads) {
				going = thread.isAlive();
				//if at least one thread is going, break out of this for loop
				if (going) break;
			}
			//Sleep for a bit and wait for the threads to finish.
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				U.p("getMatches in ScoringEngine thread was interrupted!");
				e.printStackTrace();
			}
		}
		//calculate size of combined ArrayLists
		int size = 0;
		for (ArrayList<Match> matchCluster: matches) {
			size += matchCluster.size();
		}
		//combine matches into result and return
		ArrayList<Match> out = new ArrayList<Match>(size);
		for (ArrayList<Match> matchCluster: matches) {
			out.addAll(matchCluster);
		}
		return out;
	}

}
