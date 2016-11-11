package org.biojava.bio.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.constants.OrganismCode;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.tools.FileConverterTools;

public class RunExternalExe
{
	static final String helloWorldC = "..\\HelloWorldC\\debug\\MotifCutt.exe"; 
	static final String motifCut = "..\\MotifCut\\debug\\MotifCut.exe ..\\..\\generic\\dm\\dm01g.fasta 10 10";
	
	static final String weederEXE = "..\\weeder\\debug\\weeder.exe ";
	//static final String weederPar = " -e 2 -R 50 ";
	static final String weederPar = " -R 50 -O SC -W 8 -e 2 -N -N -T 10";
	static final String bioProspectorEXE = "..\\bioProspector\\debug\\bioProspector.exe ";
	static final String bioProspectorPar = " -d 1 -r ";
	
	public static void mainn(String[] args)
	{
		//RunExternalExe runExternalExe = new RunExternalExe();
		
		try 
		{
			//runExternalExe.runWeeder();
			//runExternalExe.runMotifCut();
			String winName = "C:\\Documents and Settings\\Doaa1\\Desktop\\dm\\dm01g.fasta";
			runBioProspector(winName, Dataset.dm01g, "C:\\");
		}
		catch (Exception err) {
		     err.printStackTrace();
		}
	}

	public static String runMotifFinder(MotifFinder finder, String fileName, int numOfMotifs, int motifWidth, int numOfTrials, OrganismCode organism, Dataset dataset, String outFileDir) throws IOException, BioException, InterruptedException
	{
		String output = "";
		
		switch (finder)
		{
			case BP:
				output = runBioProspector(fileName, numOfMotifs, dataset, outFileDir);
				break;
			case MotifCut:				
				break;
			case MDSCAN:
				break;
			case WEEDER: runWeeder(fileName, numOfMotifs, motifWidth, organism, dataset, outFileDir);
				break;
		}
		return output;
	}
	
	/**
	 * -f inputFasta
	 * -W motifWidth
	 * -e missmatchNumber
	 * -R sequencePercentage
	 * -O organism code
	 * -T numOfMotifs
	 * -S process both strands (default is single), just add "-S" to the parameters to enable it
	 * output files are saved in the same directory as input file
	 * with same seqname+.wee  / seqname+.html   / seqname+.mix
	 * 
	 * problems: It takes soo long to finish, specify parameters to save time
	 * eg e (mutation)
	 * @throws IOException 
	 * @throws BioException 
	 * @throws InterruptedException 
	 */
	private static String runWeeder(String fileName, int numOfMotifs, int motifWidth, OrganismCode organism, Dataset dataset, String outFileDir) throws IOException, BioException, InterruptedException 
	{
		
		String line;			
		File dir = new File ("..\\weeder\\debug");
		String outFileName = fileName.substring(0, fileName.lastIndexOf('.'));
		
		String command = weederEXE + weederPar 
					+ " -f " + getLinuxFileName(fileName) ;
				//	+ " -W " + motifWidth
			//		+ " -T " + numOfMotifs
			//		+ " -O " + organism;
		

		Process p = Runtime.getRuntime().exec(command, null, dir);
		
		BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
		BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		
		//while ((line = input.readLine()) != null){
		//	System.out.println(line);
		//}
		while ((line = error.readLine()) != null){
			System.out.println(line);
		}
		
		//int err = p.waitFor();
		input.close();
		error.close();	
		
		/*
		// convert output file to Tomap's format
		if (err==0)	    
		{
			File temp = new File(outFileName);
			FileConverterTools.covertBioProspectorToTompaFile(new File(outFileName+".wee"), new File(outFileName+"we"));
			temp.delete();    	
		}
		*/
		return outFileName; 		 
	}
	
	private static void runMotifCut(String fileName, int numOfMotifs) throws IOException
	{
	
		String line;
		File dir = new File ("..\\motifcut\\debug");
		
	    Process p = Runtime.getRuntime().exec(motifCut, null, dir);
	    
	 
	    BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
	    BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	   
	    while ((line = input.readLine()) != null){
	        System.out.println(line);
	    }
	    while ((line = error.readLine()) != null){
	        System.out.println(line);
	    }
	    
	    input.close();
	    error.close();
	   // System.in.read();
		 
	}
	
	/**
	 * using (-d 1) to examine forward strand only
	 * (-r numOfMotifs) reporting numOfMotifs motif 
	 * @throws IOException
	 * @throws BioException
	 * @return output Filename in Tompa's format
	 * @throws InterruptedException 
	 */
	private static String runBioProspector(String fileName, int numOfMotifs, Dataset dataset, String outFileDir) throws IOException, BioException, InterruptedException
	{
	
		String line;	
		File dir = new File ("..\\bioProspector\\debug");
		String outFileName = fileName.substring(fileName.lastIndexOf('\\'), fileName.lastIndexOf('.'));
		String tempFileName = outFileDir + outFileName + ".bp.temp";
		String tompaFileName = outFileDir + outFileName+".bp.tompa";
		String voterFileName = outFileDir + outFileName+".bp.voter";
			
		String command = bioProspectorEXE + bioProspectorPar + numOfMotifs +" -i " + getLinuxFileName(fileName) +
									" -o " + getLinuxFileName(tempFileName);
		
	    Process p = Runtime.getRuntime().exec(command, null, dir);
	    
	    BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
	    BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	   
	    while ((line = input.readLine()) != null){
	        System.out.println(line);
	    }
	    while ((line = error.readLine()) != null){
	        System.out.println(line);
	    }
	    	       
	    int err = p.waitFor();
	    input.close();
	    error.close();	
	    
	    // convert output file to Tomap's format
	    if (err==0)	    
	    {
	    	File temp = new File(tempFileName);	 
	    	Vector<Motif> motifs = FileConverterTools.covertBioProspectorToTompaFile(temp, new File(tompaFileName), dataset);
	    	FileConverterTools.writeMotifVoter(motifs, new File(voterFileName));
	    	temp.delete();	    	
	    }
		return tompaFileName; 
	}
	
	private static String runBioProspector(String fileName, Dataset dataset, String outFileDir) throws IOException, BioException, InterruptedException
	{		
		return runBioProspector(fileName, 1, dataset, outFileDir);
	}
	
	// surround names having spaces with quotes
	private static String getLinuxFileName(String winName)
	{
		String linuxName = winName;
		int end = 0;
		
		while (linuxName.indexOf(' ', end) != -1)
		{
			int loc = linuxName.indexOf(' ', end);
			int start = linuxName.substring(0, loc).lastIndexOf('\\') +1;
			end = loc + linuxName.substring(loc, linuxName.length()).indexOf('\\');
			linuxName = linuxName.substring(0,start)+"/'"+linuxName.substring(start,end)+"/'"+linuxName.substring(end,linuxName.length());			
		}		
		
		return linuxName;
	}
	
}
