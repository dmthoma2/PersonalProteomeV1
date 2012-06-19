package PersonalProteome.PeptideAnalysisTool;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import PersonalProteome.U;

public class AnalysisProperties {

	
	//File/Directory Locations
	public static String annotationFile;
	public static String outputDir;
	public static String peptideListFile;
	
	
	//Extra information required for splice junction check
	public static String outputFileName;
	public static String chrmDir;
	
	public static boolean checkSplice = false;
	


	
	
	/**
	 * loadProperties Loads in the properties from the propertiesfile and stores them in variables in the AnalysisProperties class.
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
	/**
	 * Each from the properties file is brought in and put in its variable.
	 * @param line A line from the properties file.
	 */
	private static void setPropertyFromString(String line) {
		line = line.trim();
		
		/* ignore blank lines */
		if (line.equals("")) return;
		
		/* ignore comments */
		if (line.startsWith("//")) return;
		if (line.startsWith("#")) return;
		
		/* ignore lines that do not have a space in them */
		if (line.indexOf(" ") == -1) return;
		
		/* getting the property name and the propert value */
		String propertyName = line.substring(0, line.indexOf(" "));
		String propertyValue = line.substring(line.indexOf(" ") + 1, line.length());
		
		//Files
		if (propertyName.equals("annotationFile")) 
			annotationFile = propertyValue.trim();
		if (propertyName.equals("outputDir")) 
			outputDir = propertyValue.trim();
		if(propertyName.equals("peptideListFile"))
			peptideListFile = propertyValue.trim();
		
		//Splice Junctions Check
		if (propertyName.equals("outputFileName")) 
			outputFileName = propertyValue.trim();
		if(propertyName.equals("chrmDir"))
			chrmDir = propertyValue.trim();

		
		//Determine which operation to perform
		if (propertyName.equals("checkSplice"))
			checkSplice = Boolean.valueOf(propertyValue.trim());


		
	}
}
