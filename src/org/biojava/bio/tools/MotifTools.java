package org.biojava.bio.tools;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.motif.PositionWeightMatrix;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.SymbolList;

/* @author Doaa Altarawy */


public class MotifTools 
{

/*-----------------------------------------------------------------------------*/
	
	/**
	 * Todo : need recoding, use removeRedundant method (compoundLocation) in the extraction step
	 */
	// get the set of finalSites then get the final motif
	public static Motif extractSitesMultiFinder(List<Motif> motifs) throws Exception
	{
		Motif finalMotif = getCenter(motifs);
		
		
		// add all sites of the max confidence motif
	//	Motif motifWithMaxConf = getConfidenceScore(motifs);		
		
		// remove sites that have only 2 contributers
		extraRefinment(finalMotif, motifs);
		

		//if (motifWithMaxConf!=null)
		//	finalMotif.merge(motifWithMaxConf);
				
		finalMotif.removeRedundant();
		
		return finalMotif;
	}
	
	public static Motif extractSitesMV(List<Motif> motifs) throws Exception
	{
		Motif finalMotif = getCenter(motifs);
		
		
		// add all sites of the max confidence motif
		Motif motifWithMaxConf = getConfidenceScore(motifs);		
		

		if (motifWithMaxConf!=null)
			finalMotif.merge(motifWithMaxConf);
				
		finalMotif.removeRedundant();
		
		return finalMotif;
	}
	
	public static Motif extractSitesOneFinder(List<Motif> motifs) throws Exception
	{
		Motif finalMotif = getPerwiseIntersection(motifs);
		
		finalMotif.removeRedundant();
		
		return finalMotif;
	}
	/*-----------------------------------------------------------------------------*/	
	public static Motif getCenter(List<Motif> motifs) throws Exception
	{
		Motif covered = new Motif(MotifFinder.MProfiler);
		if (motifs.size()>0)
			covered.setDataset(motifs.get(0).getDataset());//covered.setDataset(dataset);			
		List<Motif> temp = new ArrayList<Motif>(motifs);
		
		// get intersection covering
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			temp.remove(m);
			for (Iterator<Motif> itr2=temp.iterator(); itr2.hasNext();)
			{
				Motif m2 = itr2.next();	
				if (m.getFinder().equals(m2.getFinder()))	// intersection must be of different motifFinders 
					continue;
				for (Iterator<String> n=m.getSites().keySet().iterator(); n.hasNext();)
				{
					String seqName = n.next();
					List<Site> list1 = m.getSites().get(seqName);
					List<Site> coveredList = covered.getSites().get(seqName);
					if (coveredList==null) 
						coveredList = new ArrayList<Site>();
					for (Iterator<Site> s1=list1.iterator(); s1.hasNext();)
					{
						Site seq1 = s1.next();
						List<Site> list2 = m2.getSites().get(seqName);
						if (list2==null) continue;
						for (Iterator<Site> s2=list2.iterator(); s2.hasNext();)
						{
							Site seq2 = s2.next();
							Site intersectSeq = seq1.getLocationIntersection(seq2);  //////// check it
							
							if (intersectSeq!=null)
								coveredList.add(intersectSeq);
						}						
					}
					if (!coveredList.isEmpty())
						covered.getSites().put(seqName, coveredList);					
				}
			}
		}
		covered.removeRedundantFast();
		
		return covered;
	}
	/*-----------------------------------------------------------------------------*/
	
	/**
	 * Returns the motif with Max confidence score
	 * @throws Exception 
	 */
	/*
	private static Motif getMaxConfidenceScore(List<Motif> motifs) throws Exception
	{
		System.out.println("Start getMaxConfidenceScore for x.size()="+motifs.size());
			
		/*double x_compact = 0;
		Motif x_Center = getCenter(motifs);
		Motif best = null;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			
			double d = x_Center.sim(m);
			if (d >= x_compact)
			{
				x_compact = d;
				best = m;
			}
		}
		
		Sortable<Motif>[] list = getConfidenceScore(motifs);
		System.out.println("End getMaxConfidenceScore for x.size()="+motifs.size());
		
		if (list!= null && list.length>0)
			return list[0].obj;
		else
			return null;
		
	}
	*/
	/*-----------------------------------------------------------------------------*/
	
	private static Motif getConfidenceScore(List<Motif> motifs)
	{
		//Sortable<Motif>[] scores = new Sortable[motifs.size()];
		double maxScore = -1;
		Motif maxMotif = null;
		// for each finder, store union of locations of all other finders except this finder
		Map<MotifFinder, Motif> locations = new HashMap<MotifFinder, Motif>(MotifFinder.values().length);
		Map<MotifFinder, Motif> locationsComp = new HashMap<MotifFinder, Motif>(MotifFinder.values().length);
						
		for (Iterator<Motif> itr= motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			if (locations.get(m.getFinder())==null)
			{
				locations.put(m.getFinder(), new Motif(m));
				locationsComp.put(m.getFinder(), new Motif());
			}
			else
				locations.get(m.getFinder()).merge(m);				
		}
		
		// refine Motifs
		for (Iterator<Motif> itr=locations.values().iterator(); itr.hasNext();)
			itr.next().removeRedundantFast();
		
		for (Iterator<MotifFinder> itr= locations.keySet().iterator(); itr.hasNext();)
		{
			MotifFinder finder = itr.next();
			for (Iterator<MotifFinder> itr2= locations.keySet().iterator(); itr2.hasNext();)
			{
				MotifFinder finder2 = itr2.next();
				if (finder2.equals(finder)) continue;
				locationsComp.get(finder2).merge(locations.get(finder));	
			}
		}
		
		// refine Motifs
		for (Iterator<Motif> itr=locationsComp.values().iterator(); itr.hasNext();)
			itr.next().removeRedundantFast();
		
		
		for (Iterator<Motif> itr= motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			double intersect = 0;
			double size = m.getAllLocationsCoverageCount();
			for (Iterator<Motif> itr2=locationsComp.values().iterator(); itr2.hasNext();)
			{
				Motif m2 = itr2.next();
				if (m.getFinder().equals(m2.getFinder())) 
					intersect += m.countIntersect(m2);				
			}			
			double score = intersect / size;
			if (score>maxScore)
			{
				maxScore = score;
				maxMotif = m;
			}
			//scores[i++] = new Sortable<Motif>(m, intersect/size);
		}
		
		//Arrays.sort(scores);
		//for (i=0; i<scores.length; i++)
		//	System.out.println("Motif "+scores[i].obj.getFinder().name()+" score: "+scores[i].value+", correct="+scores[i].obj.getCorrect());
		return maxMotif==null?null:new Motif(maxMotif);
	}
	/*-----------------------------------------------------------------------------*/
	
	
	public static Motif getPerwiseIntersection(List<Motif> motifs) throws Exception
	{		
		Motif covered = new Motif(MotifFinder.MProfiler);
		if (motifs.size()>0)
			covered.setDataset(motifs.get(0).getDataset());//covered.setDataset(dataset);		
		List<Motif> temp = new ArrayList<Motif>(motifs);
		
		// get intersection covering
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			temp.remove(m);
			for (Iterator<Motif> itr2=temp.iterator(); itr2.hasNext();)
			{
				Motif m2 = itr2.next();					
				for (Iterator<String> n=m.getSites().keySet().iterator(); n.hasNext();)
				{
					String seqName = n.next();
					List<Site> list1 = m.getSites().get(seqName);
					List<Site> coveredList = covered.getSites().get(seqName);
					if (coveredList==null) 
						coveredList = new ArrayList<Site>();
					for (Iterator<Site> s1=list1.iterator(); s1.hasNext();)
					{
						Site seq1 = s1.next();
						List<Site> list2 = m2.getSites().get(seqName);
						if (list2==null) continue;
						for (Iterator<Site> s2=list2.iterator(); s2.hasNext();)
						{
							Site seq2 = s2.next();
							Site intersectSeq = seq1.getLocationIntersection(seq2);
							
							if (intersectSeq!=null)
								coveredList.add(intersectSeq);
						}						
					}
					if (!coveredList.isEmpty())
						covered.getSites().put(seqName, coveredList);					
				}
			}
		}
		
		covered.removeRedundantFast();
		
		return covered;
	}
	/*-----------------------------------------------------------------------------*/
	public static void extraRefinment(Motif motif, List<Motif> x)
	{
		Motif finalMotif = new Motif(motif.getFinder(), motif.getDataset(), "", 0);
		
		for (Iterator<String> itr=motif.getSites().keySet().iterator(); itr.hasNext();)
		{
			String seqName = itr.next();
			List<Site> sites = new ArrayList<Site>();
			for (Iterator<Site> itr2=motif.getSites().get(seqName).iterator(); itr2.hasNext();)
			{
				Site site = itr2.next();
				Location loc1 = site.location;
				Map<MotifFinder, Integer> finderExist = getEmptyFindersMap();
				
				for (Iterator<Motif> itr3=x.iterator(); itr3.hasNext();)
				{
					Motif m = itr3.next();
					Location loc2 = m.getLocations(seqName);
					if (loc2!=null && LocationTools.overlaps(loc1, loc2))
						finderExist.put(m.getFinder(), finderExist.get(m.getFinder())+1);
				}
				
				// If 3 or more finders contribute to this site accept it
				int count = 0;
				for (Iterator<MotifFinder> itr4=finderExist.keySet().iterator(); itr4.hasNext();)
					if (finderExist.get(itr4.next())>0)
						count ++;
				
				if (count > 2)
					sites.add(site);
				
			}
			finalMotif.getSites().put(seqName, sites);
		}	
		
		motif.setSites(finalMotif.getSites());		
	}
	/*-----------------------------------------------------------------------------*/
	
	public static double sim(List<Motif> motifs)
	{
		double sum = 0;
		List<Motif> list = new ArrayList<Motif>(motifs);
		int count = 0;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			list.remove(m);
			//Map<Motif, Double> mSim = sim.get(m);
			for (Iterator<Motif> itr2=list.iterator(); itr2.hasNext();)
			{
				Motif temp = itr2.next();
				sum += m.sim(temp);	
				count++;
			}
		}
		
		sum = 2 * sum + motifs.size();	// symmetric matrix ^ 2, plus diagonal sim = 1 * motifs.size()
		
		double total = (sum / Math.pow(motifs.size(), 2));
		return total;
	}

	 
	/*-----------------------------------------------------------------------------*/
	
	public static double simFinders(List<Motif> motifs)
	{
		double sum = 0;
		Map<MotifFinder, Motif> map = new HashMap<MotifFinder, Motif>(10);
		
		int count = 0;
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			if (map.containsKey(m.getFinder()))
				map.get(m.getFinder()).merge(m);
			else
				map.put(m.getFinder(), m);
		}
			
		for (Iterator<MotifFinder> itr=map.keySet().iterator(); itr.hasNext();)
		{
			Motif m = map.get(itr.next());
			for (Iterator<MotifFinder> itr2=map.keySet().iterator(); itr2.hasNext();)
			{
				Motif m2 = map.get(itr2.next());
				if (m == m2) continue;
				sum += m.sim(m2);
				count ++;
			}
			
		}
		
		sum = 2 * sum + map.keySet().size();
		
		double total = sum / (double)map.keySet().size(); //(sum / Math.pow(motifs.size(), 2));
		return total;
	}

	
	/*-----------------------------------------------------------------------------*/
	/*-----------------------------------------------------------------------------*/
	
	public static Double[] getsequenceScore(File file, List<String> seq)
	{
		PositionWeightMatrix pwm;
		Double[] score = new Double[seq.size()+1];
		score[seq.size()] = 0.0;		
		try
		{			
		    pwm = new PositionWeightMatrix(file);
		    pwm.printMatrix();
		    System.out.println("MaxPossibleScore="+pwm.maxPossibleScore);
			System.out.println("AvgScore="+pwm.avgScore);
		    int i = 0;
		    for (Iterator<String> itr = seq.iterator(); itr.hasNext();)
		    {
		    	SymbolList target = DNATools.createDNA(itr.next());
		    	score[i] = pwm.getScoreOf(target);
		    	score[seq.size()]+= score[i];
				i++;
		    }
			
		}
		catch (Exception e) {
			e.printStackTrace();
			return score;
		}
	    
		System.out.println("Sum of scores="+score[seq.size()]);
		
		
	    return score;
	}
	
	public static Map<MotifFinder, Integer> getEmptyFindersMap()
	{
		Map<MotifFinder, Integer> finderExistsFalse = new HashMap<MotifFinder, Integer>();
		for (int i=0; i<MotifFinder.values().length; i++)
			finderExistsFalse.put(MotifFinder.values()[i], 0);
		
		return finderExistsFalse;
	}
	
	public static int getNumOfFinders(List<Motif> motifs)
	{		
		Set<MotifFinder> finders = new HashSet<MotifFinder>();
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			finders.add(m.getFinder());			// This is set, elements doesn't repeat				
		}
		
		return finders.size();
	}
}
