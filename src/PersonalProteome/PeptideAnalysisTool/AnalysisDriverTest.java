package PersonalProteome.PeptideAnalysisTool;

import PersonalProteome.U;

public class AnalysisDriverTest {

	
	
	public static void main(String[] args){
		 //Deal with truncated numbers becuase comp does not recognize 69.999999999999 = 70
		
		
		U.p("******TESTING THAT THE VALUE OF SPLIT JUNCTION LOCATIONS IS RETURNED CORRECTLY**********");
		U.p("Testing 6.0.  Should return 6.0 " + valTest(6.0));
		U.p("Testing 6.00000001.  Should return 6.0 " + valTest(6.00000001));
		U.p("Testing 5.99999999.  Should return 6.0 " + valTest(5.99999999));
		
		U.p("Testing 5.66666666.  Should return 5.66666666 " + valTest(5.66666666));
		U.p("Testing 5.33333333.  Should return 5.33333333 " + valTest(5.33333333));
	}
		
		
	public static double valTest(double val){
		
		 int floor = (int)val;
		 int ceiling = floor + 1;
	     double  minValue = floor + .9;
	 
	     if(val < ceiling && val >= minValue ){
	    	 val = ceiling;
	     }else{
	    	 if(val < (int)val + .1){
	    		 val = (int)val;
	    	 }
//	    	 val = floor;
	     }
		return val;
	}
}
