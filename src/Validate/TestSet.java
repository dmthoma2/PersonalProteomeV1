package Validate;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import javax.imageio.ImageIO;

import Peppy.Match;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.ScoringEngine;
import Peppy.Spectrum;
import Utilities.U;

public class TestSet {
	
	private String testName;
	private ArrayList<Spectrum> spectra;
	private ArrayList<Match> topPositiveMatches = null;
	private ArrayList<Match> positiveMatches = null;
	private ArrayList<Match> falsePositiveMatches = null;
	private ArrayList<Match> correctMatches = null;
	private ArrayList<MatchContainer> testedMatches = null;
	private int setSize = -1;
	
	//statistics
	private int trueTally = 0;
	private int falseTally = 0;
	private int trueTallyAtOnePercentError = -1;
	private double percentAtOnePercentError = -1;
	private double eValueAtOnePercentError = -1;
	
	//PR curve
	private double areaUnderPRCurve = 0;
	private String fileNameForPRCurve = "";
	
	long timeElapsed;
	String timeToComplete = "";
	double milisecondsPerSpectrum;
	
	public TestSet(String testName) {
		this.testName = testName;
		
		//load spectra for this test
		spectra = Spectrum.loadSpectraFromFolder("/Users/risk2/PeppyOverflow/tests/" + testName + "/spectra");
		setSize = spectra.size();
		
//		//load "correct" matches as defined by the test set
//		correctMatches = loadCorrectMatches();
//		
//		//use only spectra that don't have a big difference between precursor and the predicted protein mass
//		ArrayList<Spectrum> reducedSpectra = new ArrayList<Spectrum>();
//		for (Match match: correctMatches) {
//			double difference = match.getSpectrum().getPrecursorMass() - match.getPeptide().getMass();
//			if (Math.abs(difference) < Properties.peptideMassThreshold) {
//				reducedSpectra.add(match.getSpectrum());
//			}
//		}
//		spectra = reducedSpectra;
//		setSize = spectra.size();
//		
//		//add all the correct peptides to the peptide database
//		peptides.addAll(loadCorrectPeptides());
//		Collections.sort(peptides);		
	}
	
	public void findPositiveMatches(ArrayList<Peptide> peptides) {
		//get the matches
		long startTimeMilliseconds = System.currentTimeMillis();
		positiveMatches = (new ScoringEngine(peptides, spectra, null)).getMatches();
		long stopTimeMilliseconds = System.currentTimeMillis();
		timeElapsed = stopTimeMilliseconds - startTimeMilliseconds;
		timeToComplete = U.millisecondsToString(timeElapsed);
		milisecondsPerSpectrum = (double) timeElapsed / setSize;
		
		//Sort matches by e value	
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(positiveMatches);
		
		//See which of the positive matches are true
		testedMatches = new ArrayList<MatchContainer>();
		for (Match match: positiveMatches) {
			if (match == null) U.p("null match?  What?");
			testedMatches.add(new MatchContainer(match));
		}
		//Sort MatchContainer by e value (default)	
		Collections.sort(testedMatches);
		
		//find the #1 ranked match for each spectrum
		topPositiveMatches = new ArrayList<Match>();
		for (Match match: positiveMatches) {
			if (match.getRank() == 0) {
				topPositiveMatches.add(match);
			}
		}
			
		
		
		calculateStastics();
	}


	/**
	 * Here the peptide list is probably from a reverse database.
	 * We are wanting to see what matches we find to a database that does not
	 * have the correct answer.
	 * @param peptides
	 */
	public void findFalsePositiveMatches(ArrayList<Peptide> peptides) {
		Properties.maximumNumberOfMatchesForASpectrum = 1;
		//get the matches
		falsePositiveMatches = (new ScoringEngine(peptides, spectra, null)).getMatches();
		
		//Sort matches by e value	
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(falsePositiveMatches);
		
	}
	
	/**
	 * This should be called only after our set of matches has been found
	 */
	private void calculateStastics() {
		//stats for tested matches
		boolean onePercenThresholdHasBeenReached = false;
		for (MatchContainer match: testedMatches) {
			if (match.isTrue()) {
				trueTally++;
			} else {
				falseTally++;
			}
			if (!onePercenThresholdHasBeenReached) {
				if ((double) falseTally / setSize >= 0.01) {
					onePercenThresholdHasBeenReached = true;
					trueTallyAtOnePercentError =  trueTally;
					percentAtOnePercentError = (double) trueTallyAtOnePercentError / setSize;
					eValueAtOnePercentError = match.getEValue();
				}
			}
		}
		
		generatePrecisionRecallCurve();
	}

	private void generatePrecisionRecallCurve() {
			int width = 500;
			int height = 500;
			//setting up Graphics context
			BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			Graphics2D g = bdest.createGraphics();
			g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
			g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
			g.setColor(Color.white);
			g.fillRect(0,0,width,height);
			
			//setting the line color
			g.setColor(Color.red);
			
			int x1 = 0;
			int x2 = 0;
			int x2_old = 0;
			int y1 = 0;
			int y2 = 0;
			int trueCount = 0;
			double precision = 0; 
			double recall = 0;
			double recallPrevious = 0;
			double area = 0;
			for(int i = 0; i < testedMatches.size(); i++) {
				MatchContainer match = testedMatches.get(i);
				if (match.isTrue()) {
					trueCount++;
				}
				
				
				precision = (double) trueCount / (i + 1);
				recallPrevious = recall;
				recall = (double) trueCount / getSetSize();	
				
				area += (recall - recallPrevious) * precision;
				
				
				x2 = (int) (recall * width);
				y2 = (int) ((1.0 - precision) * height);
				
				//in case we are moving so little we are not advancing
				if (x1 + 1 >= x2) {
					continue;
				} else {
					x1 = x2_old;
					
	//				U.p(trueCount + ", " + precision);
					g.setColor(Color.red);
	//				g.setStroke(new BasicStroke(2.0f));
					g.drawLine(x1, y1, x2, y2);
					//let's fill in the area under the line, yes?
					Polygon polygon = new Polygon();
					polygon.addPoint(x1, y1);
					polygon.addPoint(x2, y2);
					polygon.addPoint(x2, height);
					polygon.addPoint(x1, height);
					g.setColor(new Color(255,0,0,128));
					g.fillPolygon(polygon);
					
					//updating variables
					x2_old = x2;
					y1 = y2;
				}
				
			}
			//draw final line straight down (just for looks)
			g.setColor(Color.red);
	//		g.setStroke(new BasicStroke(2.0f));
			g.drawLine(x2, y1, x2, height);
			
			areaUnderPRCurve = area;
			
	//		g.drawLine(x2, y2, width, height);
			
			try {
				File testDirectory = new File(Properties.validationDirectory, testName);
				testDirectory.mkdirs();
				File imageFile = new File(testDirectory, "/precision-recall.jpg");
				fileNameForPRCurve = testName + "/" + imageFile.getName();
				ImageIO.write(bdest,"JPG",imageFile);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	private ArrayList<Peptide> loadCorrectPeptides() {
		ArrayList<Peptide> correctPeptides = new ArrayList<Peptide>();
		for(Spectrum spectrum: spectra) {
			//find the file for the correct peptide
			File spectrumFile = spectrum.getFile();
			File testFolder = spectrumFile.getParentFile().getParentFile();
			File peptideFolder = new File(testFolder, "peptides");
			File peptideFile = new File(peptideFolder, spectrumFile.getName());
			
			//load in the correct peptide string
			boolean validPeptideFile = true;
			String correctAcidSequence = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFile));
				//read the first line;
				correctAcidSequence = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				validPeptideFile = false;
				e.printStackTrace();
			} catch (IOException e) {
				validPeptideFile = false;
				e.printStackTrace();
			}
			
			//testing that we've got a valid peptide file
			if (correctAcidSequence == null) {
				validPeptideFile = false;
			}
			correctAcidSequence = correctAcidSequence.trim();
			if (correctAcidSequence.equals("")) {
				validPeptideFile = false;
			}
			//adding to the array list
			if (validPeptideFile) {
				correctPeptides.add(new Peptide(correctAcidSequence));
			}
		}
		return correctPeptides;
	}
	
	private ArrayList<Match> loadCorrectMatches() {
		ArrayList<Match> correctMatches = new ArrayList<Match>();
		for(Spectrum spectrum: spectra) {
			//find the file for the correct peptide
			File spectrumFile = spectrum.getFile();
			File testFolder = spectrumFile.getParentFile().getParentFile();
			File peptideFolder = new File(testFolder, "peptides");
			File peptideFile = new File(peptideFolder, spectrumFile.getName());
			
			//load in the correct peptide string
			boolean validPeptideFile = true;
			String correctAcidSequence = "";
			try {
				BufferedReader br = new BufferedReader(new FileReader(peptideFile));
				//read the first line;
				correctAcidSequence = br.readLine();
				//close;
				br.close();
			} catch (FileNotFoundException e) {
				validPeptideFile = false;
				e.printStackTrace();
			} catch (IOException e) {
				validPeptideFile = false;
				e.printStackTrace();
			}
			
			//testing that we've got a valid peptide file
			if (correctAcidSequence == null) {
				validPeptideFile = false;
			}
			correctAcidSequence = correctAcidSequence.trim();
			if (correctAcidSequence.equals("")) {
				validPeptideFile = false;
			}
			//adding to the array list
			if (validPeptideFile) {
				correctMatches.add(new Match(spectrum, new Peptide(correctAcidSequence), null));
			}
			
//			if (spectrum.getFile().getName().equals("T10707_Well_H13_1768.77_19185.mgf..pkl")) {
//				U.p ("Valid? " + validPeptideFile);
//			}
		}
		return correctMatches;
	}
	
	
	public String getName() {
		return testName;
	}

	public int getTrueTally() {
		return trueTally;
	}

	public int getFalseTally() {
		return falseTally;
	}

	public int getTrueTallyAtOnePercentError() {
		return trueTallyAtOnePercentError;
	}

	public double getEValueAtOnePercentError() {
		return eValueAtOnePercentError;
	}

	public double getMilisecondsPerSpectrum() {
		return milisecondsPerSpectrum;
	}

	public double getAreaUnderPRCurve() {
		return areaUnderPRCurve;
	}

	public String getFileNameForPRCurve() {
		return fileNameForPRCurve;
	}

	public long getTimeElapsed() {
		return timeElapsed;
	}

	public String getTimeToComplete() {
		return timeToComplete;
	}

	public int getSetSize() {
		return setSize;
	}

	public double getPercentAtOnePercentError() {
		return percentAtOnePercentError;
	}

	public ArrayList<MatchContainer> getTestedMatches() {
		return testedMatches;
	}
	
	public double getEValueAtPercentForwards(double percent) {
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(positiveMatches);
		int level = (int) (positiveMatches.size() * percent);
		return positiveMatches.get(level).getEValue();
	}
	
	public double getEValueAtPercentReverse(double percent) {
		Match.setSortParameter(Match.SORT_BY_E_VALUE);
		Collections.sort(falsePositiveMatches);
		int level = (int) (falsePositiveMatches.size() * percent);
		return falsePositiveMatches.get(level).getEValue();
	}


}
