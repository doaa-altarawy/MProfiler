package org.biojava.bio.motif;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.MuscleAlignment;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.UniformDistribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.tools.FileConverterTools;
import org.biojavax.bio.seq.RichSequence.IOTools;

import com.sun.jndi.url.dns.dnsURLContext;
import com.sun.org.apache.bcel.internal.generic.RETURN;

public class CopyOfPositionWeightMatrix implements WeightMatrix
{
	private Distribution [] columns; 
	private Alphabet alpha;
	  
	public double [] columnsInfoContent;
	
	public double avgScore = 0.0;
	public double maxPossibleScore = 0.0;
	public double totalInfoContent = 0;
	
	int numOfWords;
	
	public CopyOfPositionWeightMatrix(Motif motif, String count) throws Exception
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
	
	public CopyOfPositionWeightMatrix(File file) throws Exception
	{						
		Vector<SymbolList> list;
		
		if (file.getName().endsWith(".voter"))
			list = readVoter(file);
		else
			list = readFASTA(file);
		
	    getColumns(list, "0");	  
	}
	
	private void getColumns(Vector<SymbolList> list, String n) throws Exception
	{		
		if (list.size()<1)
			return;
		numOfWords = list.size();
		Alignment al = new MuscleAlignment(list, n);
	    columns = DistributionTools.distOverAlignment(al, true, 0.1);
		   
	    columnsInfoContent = new double[columns.length];
	    
	    this.alpha = DNATools.getDNA();//columns[0].getAlphabet();
    
	   	for(int c = 0; c < columns.length; c++) 
		{
	   		for (Iterator itr = DNATools.getDNA().iterator(); itr.hasNext();)
		    {	    	
		    	Symbol s = (Symbol)itr.next();
		    	double d = columns[c].getWeight(s);
		    	double d2 = Math.log(d);
		    	if (!Double.isInfinite(d2) && !Double.isNaN(d2))
		    		columnsInfoContent[c]+= d * d2;  		    	
		    }
	   		columnsInfoContent[c]+=2;
	   		totalInfoContent+= columnsInfoContent[c];
	   		
	   		if(columns[c].getAlphabet() != alpha)
	   		{
		        throw new IllegalAlphabetException(
		          "All columns must emit the same alphabet. Expecting " +
		          alpha.getName() + ", but found " + columns[c].getAlphabet().getName());
		     }	
		}
	   	
	   	printMatrix();
	    maxPossibleScore = maxScore(); 
	    
	    for (Iterator<SymbolList> itr = al.symbolListIterator(); itr.hasNext();)
	    	avgScore+= getScoreOf(itr.next());
	    avgScore = avgScore / (double)list.size();
	}
	
		
	public double getMAPScore() throws IllegalSymbolException
	{
		double score = 0;
		if (columns==null)
			return Double.NaN;
		
		for(int c = 0; c < columns.length; c++) 
	 	{
	    	score+= columnsInfoContent[c];
	    }	    	
	    		
		score = score * Math.log(numOfWords) / (double)columns();
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
			
			score = columns[i].getWeight(targetSeq.symbolAt(i+1));
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
			for (Iterator<Symbol> itr = DNATools.getDNA().iterator(); itr.hasNext();)
			{
				double score = columns[i].getWeight(itr.next());
				if (score>max)
					max = score;
			}			
			totalScore += max;			
		}
		
		return totalScore;
	}
	
	public int columns()
	{		
		return this.columns.length;
	}
	
	public Distribution getColumn(int column) 
	{
	    return columns[column];
	}
	
	public Alphabet getAlphabet() 
	{
	    return alpha;
	}
	
	public void printMatrix() throws IllegalSymbolException
	{
		System.out.print("\t");
	    for (int i=0; i<columns(); i++)
	     	System.out.print(i+"\t");
	    System.out.println();
	   	   
	    for (Iterator itr = DNATools.getDNA().iterator(); itr.hasNext();)
	    {	    	
	    	Symbol s = (Symbol)itr.next();
	    	System.out.print(s.getName().toUpperCase().charAt(0)+"\t");	  
	    	for(int c = 0; c < columns.length; c++) 
		 	{
		    	System.out.print(String.format("%.2f", columns[c].getWeight(s))+"\t");			    
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
