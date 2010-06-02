package SpectralVisualizer;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import Peppy.Definitions;
import Peppy.Peak;
import Peppy.Peptide;
import Peppy.Properties;
import Peppy.Spectrum;
import Reports.HTMLReporter;
import Utilities.U;

public class SpectralVisualizer {
	
	static Color yIonColor = Color.green;
	static Color bIonColor = Color.red;
	
	public static void main(String args[]) {
		U.printUserDirectory();
		U.p("drawing images for spectra");
		ArrayList<Spectrum> spectra = Spectrum.loadSpectraFromFolder("spectra");
		ArrayList<Peptide> peptides = loadPeptides("peptides.txt");
		generateFullReport(spectra, peptides);
	}
	
	public static void generateFullReport(ArrayList<Spectrum> spectra, ArrayList<Peptide> peptides) {
		File reportFile = new File(Properties.reportDirectory, "index.html");
		try {
			//create our report directories
			Properties.reportDirectory.mkdirs();
			File imageFolder = new File(Properties.reportDirectory, "images");
			imageFolder.mkdirs();
			
			//set up our main index file
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(reportFile)));
			
			//print headers
			HTMLReporter.appendFile(pw, Reports.Properties.reportWebHeaderFile);
			
			StringBuffer sb;
			for (int i = 0 ; i <spectra.size(); i++) {
				//mark potential peptide peaks
				markMatchingIons(spectra.get(i), peptides.get(i));
				
				//make an image file of the spectrum
				File imageFile = new File(imageFolder, i + ".jpg");
				drawSpectrum(spectra.get(i), 1000, 300, imageFile);
				
				//include in report file
				pw.println("<img src=\"images/" + i + ".jpg\">");
				pw.println("<br>");
				pw.println(peptides.get(i).getAcidSequence());
				pw.println("<p>");
			}
			
			//print headers
			HTMLReporter.appendFile(pw, Reports.Properties.reportWebHeaderFile);

			pw.flush();
			pw.close();

			
		} catch (FileNotFoundException e) {
			U.p("could not find file: " + reportFile.getName());
			e.printStackTrace();
		} catch (IOException e) {
			U.p("could not read file: " + reportFile.getName());
			e.printStackTrace();
		}
	}
	
	public static ArrayList<Peptide> loadPeptides(String peptideFileName) {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(peptideFileName));
			String line = br.readLine();
			while (line != null) {
				line = line.trim();
				if (!line.equals("")) {
					peptides.add(new Peptide(line));
				}
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return peptides;
	}
	
	public static void drawSpectrum(Spectrum spectrum, int width, int height, File dest) throws IOException {
		
		int xLoc, yLoc;
		
		//setting up Graphics context
		BufferedImage bdest = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bdest.createGraphics();
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_SPEED));
		g.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF));
		g.setColor(Color.BLACK);
		g.fillRect(0,0,width,height);
		
		//getting maximum spectrum value and intensity
		double maxValue = spectrum.getMaxMass();
		double maxIntensity = spectrum.getCalculatedMaxIntensity();
		
		//draw the lines for our spectra
		ArrayList<Peak> peaks = spectrum.getPeaks();
		for (Peak peak: peaks) {
			g.setColor(peak.getColor());
			xLoc = (int) (peak.getMass()  * width / maxValue);
			yLoc = (int) (height - ((peak.getIntensity() / maxIntensity) * height));
			g.drawLine(xLoc, yLoc, xLoc, height);
			if (!peak.getColor().equals(Color.gray)) {
				if (yLoc < 10) yLoc += 10;
				g.drawString("" + peak.getMass(), xLoc, yLoc);
			}
		}
		ImageIO.write(bdest,"JPG",dest);
	}
	

	
	public static void markMatchingIons(Spectrum spectrum, Peptide peptide) {
		String peptideString = peptide.getAcidSequence();

		int i;
		double theoreticalPeakMass, peakMass;
		int peakIndex, seqIndex;
		
		//we want -1 because most of these spectra will have a match with 
		//the last theoretical peak
		int peptideLengthMinusOne = peptideString.length() - 1;
		
		double [] bIonMatchesWithHighestIntensity = new double[peptideString.length()];
		for (i = 0; i < peptideString.length(); i++) bIonMatchesWithHighestIntensity[i] = 0.0;
		double [] yIonMatchesWithHighestIntensity = new double[peptideString.length()];
		for (i = 0; i < peptideString.length(); i++) yIonMatchesWithHighestIntensity[i] = 0.0;

		
		//find the ranges around our theoretical peptides where we
		//count spectrum peaks
		double [] theoreticalPeaksLeft = new double[peptideLengthMinusOne];
		double [] theoreticalPeaksRight = new double[peptideLengthMinusOne];
		
		
		/* y-ion  */
		//computing the left and right boundaries for the ranges where our peaks should land
		theoreticalPeakMass = peptide.getMass() + Properties.yIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass -= Definitions.getAminoAcidWeightMono(peptideString.charAt(i));
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = spectrum.getPeakCount() - 1;
		seqIndex = 0;
		while (peakIndex >= 0) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass < theoreticalPeaksLeft[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				if (yIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					spectrum.getPeak(peakIndex).setColor(yIonColor);
				}
			}
			
			peakIndex--;
		}
			
		
		/* b-ion  */
		theoreticalPeakMass = Properties.bIonDifference;
		for (i = 0; i < peptideLengthMinusOne; i++) {
			theoreticalPeakMass += Definitions.getAminoAcidWeightMono(peptideString.charAt(i));
			theoreticalPeaksLeft[i] = theoreticalPeakMass - Properties.peakDifferenceThreshold;
			theoreticalPeaksRight[i] = theoreticalPeakMass + Properties.peakDifferenceThreshold;
		}
		
		peakIndex = 0;
		seqIndex = 0;
		while (peakIndex < spectrum.getPeakCount()) {
			peakMass = spectrum.getPeak(peakIndex).getMass();
			while (peakMass > theoreticalPeaksRight[seqIndex]) {
				seqIndex++;
				if (seqIndex == peptideLengthMinusOne) break;
			}
			if (seqIndex == peptideLengthMinusOne) break;
			if (peakMass >= theoreticalPeaksLeft[seqIndex] && peakMass <= theoreticalPeaksRight[seqIndex]) {
				if (bIonMatchesWithHighestIntensity[seqIndex] < spectrum.getPeak(peakIndex).getIntensity()) {
					spectrum.getPeak(peakIndex).setColor(bIonColor);
				}
			}
			
			peakIndex++;
		}
		
	}
	

}