package org.biojava.bio.motif;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.tools.TompaDataset;
import org.biojava.utils.ChangeVetoException;




public class FuzzyMotif extends  Motif
{		

	public FuzzyMotif()
	{
		super();
	}
	
	public FuzzyMotif(MotifFinder finder)
	{
		super(finder);
	}
	
	public FuzzyMotif(Dataset dataset)
	{
		super(dataset);
	}
	
	public FuzzyMotif(Dataset dataset, String comment, long totalLen)
	{
		super(dataset, comment, totalLen);
	}
	
	public FuzzyMotif(MotifFinder finder, Dataset dataset, String comment, long totalLen)
	{
		super(finder, dataset, comment, totalLen);
	}
	
	public FuzzyMotif(Motif m)
	{
		super(m);
	}
	
	public FuzzyMotif(Motif m, double w)
	{
		this.finder = m.finder;
		this.dataset = m.dataset;
		this.comment = m.comment;
		this.totalLength = m.totalLength;
		sites = new HashMap<String,List<Site>>(m.sites.size());
		
		for (Iterator<String> itr=m.getSites().keySet().iterator(); itr.hasNext();)
		{
			String seq = itr.next();
			List<Site> list = new ArrayList<Site>(m.getSites().get(seq).size());
			for (Iterator<Site> itr2=m.getSites().get(seq).iterator(); itr2.hasNext();)
			{
				Site site = itr2.next();
				list.add(new FuzzySite(site, w));
			}
			sites.put(seq, list);
		}
	}
	
	/*----------------------------------------------------------------------------------*/
	
	// rewrite
	public FuzzyMotif intersect(Motif m, double weight) throws NoSuchElementException, BioException, IOException
	{
	
		FuzzyMotif finalMotif = new FuzzyMotif(finder, dataset, "", 0);				
		for (Iterator<String> itr=getSites().keySet().iterator(); itr.hasNext();)
		{
			String seqName = itr.next();			
			List<Site> sites1 = getSites().get(seqName);
			List<Site> sites2 = m.getSites().get(seqName);
			
			if (sites2 != null)
			{
				List<Site> finalSites = new ArrayList<Site>();
				for (Iterator<Site> site1Itr=sites1.iterator(); site1Itr.hasNext();)
				{
					FuzzySite site1 = (FuzzySite)site1Itr.next();
					for (Iterator<Site> site2Itr=sites2.iterator(); site2Itr.hasNext();)
					{
						Site site2 = site2Itr.next();
						Site temp = site1.getIntersection(site2, weight);
						if (temp!=null)
							finalSites.add(temp);
					}
				}
				finalMotif.getSites().put(seqName, finalSites);	
			}
		}
		
		return finalMotif;
	}
		
	
	// this motif is fuzzy, but m is not.. weight is the weight associated with m
	public double sim(Motif m, double weight) throws IndexOutOfBoundsException, ChangeVetoException, BioException
	{
		double total = 0;
		double union = 0;
		double intersection = 0;
		Set<String> seqNames = new HashSet<String>(getSites().keySet());
		seqNames.addAll(m.getSites().keySet());
		
		for (Iterator<String> itr=getSites().keySet().iterator(); itr.hasNext();)
		{
			String seqName = itr.next();			
			List<Site> sites1 = getSites().get(seqName);
			List<Site> sites2 = m.getSites().get(seqName);
			
			if (sites1!=null && sites2 == null)
			{
				for (Iterator<Site> site1Itr=sites1.iterator(); site1Itr.hasNext();)
				{
					FuzzySite site1 = (FuzzySite)site1Itr.next();				
					union+=	site1.getWeightedCoverage();
				}
			}
			else if (sites1==null && sites2 != null)
			{
				union+= (double)LocationTools.coverage(m.getLocations(seqName)) * weight;
			}
			else if (sites1!=null && sites2 != null)
			{				
				for (Iterator<Site> site1Itr=sites1.iterator(); site1Itr.hasNext();)
				{
					FuzzySite site1 = (FuzzySite)site1Itr.next();
					for (Iterator<Site> site2Itr=sites2.iterator(); site2Itr.hasNext();)
					{
						Site site2 = site2Itr.next();
						intersection += site1.getIntersectionCount(site2, weight);
						union += site1.getUnionCount(site2, weight);
						
					}
				}					
			}
		}
		if (union != 0)
			total += ((double)intersection/(double)union);
		
		return total;
	}
	
	@SuppressWarnings("serial")
	public static class FuzzySite extends Motif.Site
	{
		double weight;
		
		public FuzzySite(Site site, double weight)
		{
			dna = site.dna;
			strand = site.strand;
			source = site.source;
			location = site.location;
			annotation = site.annotation;
			
			this.weight = weight;
		}
		
		public Site getIntersection(Site s2, double weight2) throws IndexOutOfBoundsException, ChangeVetoException, BioException
		{
			Site intersect = null;
			
			if (LocationTools.overlaps(location, s2.location))
			{
				Location loc = LocationTools.intersection(location, s2.location);		
				intersect = new FuzzySite(getSubSite(loc.getMin(), loc.getMax()), (weight+weight2)/2.0);
			}
			
			return intersect;
		}
		
		public double getIntersectionCount(Site s2, double weight2) throws IndexOutOfBoundsException, ChangeVetoException, BioException
		{
			if (LocationTools.overlaps(location, s2.location))
			{
				Location loc = LocationTools.intersection(location, s2.location);		
				return ((double)LocationTools.coverage(loc))*(weight+weight2)/2.0;
			}
			
			return 0.0;
		}
		
		public double getUnionCount(Site s2, double weight2) throws IndexOutOfBoundsException, ChangeVetoException, BioException
		{
			
			if (s2==null)
				return getWeightedCoverage();			
			else if (LocationTools.overlaps(location, s2.location))//////// revise: give non intersecting bits different weight?
			{
				Location loc = LocationTools.union(location, s2.location);		
				return ((double)LocationTools.coverage(loc))*(weight+weight2)/2.0;
			}
			else
				return getWeightedCoverage()+ s2.getWeightedCoverage(weight2);
			
		}
		
		public double getWeightedCoverage()
		{
			return getWeightedCoverage(weight);
		}
	}
}
