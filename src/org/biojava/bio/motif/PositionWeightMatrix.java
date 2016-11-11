package org.biojava.bio.motif;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.MuscleAlignment;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.tools.FileConverterTools;
import org.biojavax.bio.seq.RichSequence.IOTools;


public class PositionWeightMatrix
{
	private double [][] columns; //  j =1,2,... , w, i = {A,C,G,T}, column[i][j]= freq
	private static Map<Symbol, Integer> symbols;
	  
	public double [] columnsInfoContent;
	
	public double avgScore = 0.0;
	public double maxPossibleScore = 0.0;
	public double totalInfoContent = 0;
	
	int numOfWords;
	boolean saveCounts = false;
	
	static
	{
		symbols = new HashMap<Symbol, Integer>();
		int i=0;
		for (Iterator<Symbol> itr=DNATools.getDNA().iterator(); itr.hasNext();)
			symbols.put(itr.next(), i++);
	}
	
	public PositionWeightMatrix(Motif motif, String count) throws Exception
	{						
		List<Site> sites = motif.getAllSites();
		if (sites.size()<1)
			return;
		
		Vector<SymbolList> list = new Vector<SymbolList>();		
		for (Iterator<Site> itr=sites.iterator(); itr.hasNext();) 
		{
			String s = itr.next().dna;			
			if (s.length()<6) continue; ////////////
			SymbolList subSeq = DNATools.createDNASequence(s, "");	     
	        list.add(subSeq); 
	    }
				
		getColumns(list, count);	
	}
	
	public PositionWeightMatrix(File file) throws Exception
	{						
		Vector<SymbolList> list;
		
		if (file.getName().endsWith(".voter"))
			list = readVoter(file);
		else
			list = readFASTA(file);
		saveCounts = true;
	    getColumns(list, "0");	  
	}
	
	private void getColumns(Vector<SymbolList> list, String n) throws Exception
	{		
		if (list.size()<1)
			return;
		numOfWords = list.size();
		Alignment al = new MuscleAlignment(list, n);
	   
		getPWM(al);
		   
	    columnsInfoContent = new double[columns.length];
	        
	   	for(int c = 0; c < columns.length; c++) 
		{
	   		for (int i=0; i<symbols.size(); i++)
		    {	   		    	
		    	double d = columns[c][i]; 	
		    	double d2 = Math.log10(4.0 * d)/ Math.log10(2.0);
		    	columnsInfoContent[c]+= d * d2;  		    	
		    }
	   		//columnsInfoContent[c]+=2.0;
	   		totalInfoContent+= columnsInfoContent[c];
	   			   		
		}
	   	
	   	totalInfoContent = totalInfoContent / (double)columns();
	   	//printMatrix();
	 //   maxPossibleScore = maxScore(); 
	    
	  //  for (Iterator<SymbolList> itr = al.symbolListIterator(); itr.hasNext();)
	  //  	avgScore+= getScoreOf(itr.next());
	  //  avgScore = avgScore / (double)list.size();
	}
	
	private void getPWM(Alignment a) throws FileNotFoundException, IOException
	{		
		List seqs = a.getLabels();
		columns = new double[a.length()][symbols.size()];
		double[] sum = new double[a.length()];
		
		// get counts
		for(int i = 0; i < a.length(); i++)
		{// For each position			  		   
		    for(Iterator j = seqs.iterator(); j.hasNext();)
		    {// of each sequence		       
		       Symbol s = a.symbolAt(j.next(),i + 1);		
		       if(s == null)
		         continue;
		
		       if (DNATools.getDNA().getGapSymbol().equals(s))
		    	   continue;
		       columns[i][symbols.get(s)]++;
		       sum[i]++;
		     }	    
		 }
		
		if (saveCounts)
			saveCountsToFile();
		
		
		
		// calculate relative entropy (Information content) for each position
		// uses pseudo counts
		for(int i = 0; i < a.length(); i++)
		{	
			double eps = Math.sqrt(sum[i])/ 4.0;
			double d1 = eps / 4.0;
			double d2 = (double)sum[i] + eps;
			
		    for(int j=0; j<symbols.size(); j++)
		    {		      	
		       columns[i][j] = (columns[i][j]+d1) / d2;		       
		    }	    
		 }
		
	}
	
	public double getMAPScore() throws IllegalSymbolException
	{
		double score = 0;
		if (columns==null)
			return Double.NaN;
		
		for(int c = 0; c < columns.length; c++) 
	 	{
			for (int i=0; i<symbols.size(); i++)
		    {	   		    	
		    	double d = columns[c][i]; 	
		    	double d2 = Math.log10(d) / Math.log10(2.0);
		    	score+= d * d2;  		    	
		    }			
	    }	    	
	    		
		score = score * Math.log10(numOfWords) / Math.log10(2.0)/ (double)columns();
		//System.out.println("MAP="+score);
		
		//if (Double.isNaN(score))
		//	return 0;
		return score;
	}
	
		
	private Vector<SymbolList> readFASTA(File file) throws NoSuchElementException, BioException, FileNotFoundException
	{
		BufferedReader br = new BufferedReader(new FileReader(file));
		
		SequenceIterator itr = IOTools.readFastaDNA(br, null);
			
		Vector<SymbolList> list = new Vector<SymbolList>();
		for (int i=0; itr.hasNext(); i++) 
		    list.add(itr.nextSequence()); 
	    	   
		return list;
	}
	
	private Vector<SymbolList> readVoter(File file) throws IOException, BioException 
	{				
		Vector<Motif> motifs = FileConverterTools.readMotifVoter(file, true);
		Iterator<Site> sites = motifs.firstElement().getAllSites().iterator();
		
		Vector<SymbolList> list = new Vector<SymbolList>();
		for (int i=0; sites.hasNext(); i++) 
		{
			String s = sites.next().dna;				
			SymbolList subSeq = DNATools.createDNASequence(s, "");	     
			list.add(subSeq); 
	    }
		
		return list;
	}
	public Double getScoreOf(SymbolList targetSeq) throws BioException
	{
		Double score = 0.0, totalScore = 0.0;
		
		for (int i = 0; i < this.columns(); i++)
		{
			if (targetSeq.length() < i+1) break;
			if (DNATools.getDNA().getGapSymbol().equals(targetSeq.symbolAt(i+1)))
				continue;
			score = columns[i][symbols.get(targetSeq.symbolAt(i+1))];
			totalScore += score;			
		}
		//System.out.println("Score("+targetSeq.seqString()+")= "+totalScore.doubleValue());
		
		return totalScore;
	}
	
	private Double maxScore() throws BioException
	{
		Double totalScore = 0.0;
		double max;		
		for (int i = 0; i < this.columns(); i++)
		{
			max = -1;
			for (Iterator<Symbol> itr = symbols.keySet().iterator(); itr.hasNext();)
			{
				double score = columns[i][symbols.get(itr.next())];
				if (score>max)
					max = score;
			}			
			totalScore += max;			
		}
		
		return totalScore;
	}
	
	public int columns()
	{		
		if (columns==null) return -1;
		else return this.columns.length;
	}
	
	public void saveCountsToFile() throws FileNotFoundException, IOException
	{
		FileWriter out = new FileWriter(File.createTempFile("count", ".PWM", new File(".")));
		
		out.write(">m1\n");
		for(int i = 0; i <columns(); i++)
		{			  		   
		    for(int j=0; j<symbols.size(); j++)
		    {		      	
		       out.write((int)columns[i][j]+"\t");		       
		    }	    
		    out.write('\n');
		 }
		
		out.close();
	}
	
	public void printMatrix() throws IllegalSymbolException
	{
		System.out.print("\t");
	    for (int i=0; i<columns(); i++)
	     	System.out.print(i+"\t");
	    System.out.println();
	   	   
	    for (Iterator itr = symbols.keySet().iterator(); itr.hasNext();)
	    {	    	
	    	Symbol s = (Symbol)itr.next();
	    	System.out.print(s.getName().toUpperCase().charAt(0)+"\t");	  
	    	for(int c = 0; c < columns.length; c++) 
		 	{
		    	System.out.print(String.format("%.2f", columns[c][symbols.get(s)])+"\t");			    
		    }
	    	System.out.println();
	    }
	    System.out.print("------------------------------------\n\t");
	    for(int c = 0; c < columns.length; c++) 
	 	{
	    	System.out.print(String.format("%.2f", columnsInfoContent[c])+"\t");			    
	    }
	    System.out.println("\nTotalInfoContent="+totalInfoContent);
	}
	
}
