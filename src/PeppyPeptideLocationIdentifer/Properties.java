package PeppyPeptideLocationIdentifer;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import PersonalProteome.U;

/**
 * Properties file for PeppyPeptideLocationIdentifer
 * @author David "Corvette" Thomas
 *
 */
public class Properties {

	
	//File/Directory Locations
	public static String annotationFile = "No annotation file set!";
	public static String peppyFileDirectory = "No peppy file set!";
	public static String outputDir = "/output/";

	
	/**
	 * 
	 * @param fileName the name of our properties file
	 */
	public static void loadProperties(File propertiesFile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(propertiesFile));
			String line = br.readLine();
			while (line != null) {
				setPropertyFromString(line);
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			U.p("Could not find the properties file: " + propertiesFile.getName());
			U.p("Using default properties...");
		} catch (IOException e) {
			U.p("Could not read the properties file: " + propertiesFile.getName());
			e.printStackTrace();
		}
		
	}
	
	private static void setPropertyFromString(String line) {
		line = line.trim();
		
		/* ignore blank lines */
		if (line.equals("")) return;
		
		/* ignore comments */
		if (line.startsWith("//")) return;
		if (line.startsWith("##")) return;
		
		/* ignore lines that do not have a space in them */
		if (line.indexOf(" ") == -1) return;
		
		/* getting the property name and the propert value */
		String propertyName = line.substring(0, line.indexOf(" "));
		String propertyValue = line.substring(line.indexOf(" ") + 1, line.length());
		
		//Files
		if (propertyName.equals("annotationFile")) 
			annotationFile = propertyValue.trim();
		if (propertyName.equals("peppyFileDirectory")) 
			peppyFileDirectory = propertyValue.trim();
		if (propertyName.equals("outputDir")) 
			outputDir = propertyValue.trim();
		
		
	}
}
