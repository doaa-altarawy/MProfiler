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
import java.util.Vector;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.tools.TompaDataset;
import org.biojava.utils.ChangeVetoException;


public class Motif
{
	MotifFinder finder;			// motif finder name
	
	Dataset dataset;				// dataset for this motif
	String comment;
	
	long totalLength;			// the length of the seq from which the motif was extracted
	
	int correct;				// number of correct sites (for testing);
	
	Map<String, List<Site>> sites;		// <sequence name, its sites> locations are in negative, sequence starts from 0
	
	public Motif()
	{
		sites = new HashMap<String,List<Site>>();
	}
	
	public Motif(MotifFinder finder)
	{
		this.finder = finder;		
		sites = new HashMap<String,List<Site>>();
	}
	
	public Motif(Dataset dataset)
	{
		this.dataset = dataset;		
		sites = new HashMap<String,List<Site>>();
	}
	
	public Motif(Dataset dataset, String comment, long totalLen)
	{
		this.dataset = dataset;
		this.comment = comment;
		this.totalLength = totalLen;
		sites = new HashMap<String,List<Site>>();
	}
	
	public Motif(MotifFinder finder, Dataset dataset, String comment, long totalLen)
	{
		this.finder = finder;
		this.dataset = dataset;
		this.comment = comment;
		this.totalLength = totalLen;
		sites = new HashMap<String,List<Site>>();
	}
	
	public Motif(Motif m)
	{
		this.finder = m.finder;
		this.dataset = m.dataset;
		this.comment = m.comment;
		this.totalLength = m.totalLength;
		sites = new HashMap<String,List<Site>>();
		
		// deep copy
		for (Iterator<String> itr=m.sites.keySet().iterator(); itr.hasNext();)
		{
			String seq = itr.next();
			sites.put(seq, new ArrayList<Site>(m.sites.get(seq)));		
		}
	}
	
	public void addSite(Site seq)
	{
		String source = seq.source;
		List<Site> list = sites.get(source);
		if (list==null)
			list = new ArrayList<Site>();
		list.add(seq);
		sites.put(source, list);
	}
	
	public void addSite(String dna, String source, int start, int end, String direction) throws ChangeVetoException, BioException
	{
		Site seq = new Site(dna, source, start, end, direction);
		List<Site> list = sites.get(source);
		if (list==null)
			list = new ArrayList<Site>();
		list.add(seq);
		sites.put(source, list);		
	}

	public List<Site> getAllSites()
	{
		List<Site> list = new ArrayList<Site>();
		
		for (Iterator<List<Site>> itr=sites.values().iterator(); itr.hasNext();)
		{
			list.addAll(itr.next());
		}
		return list;
	}
	
	public void removeAllSites()
	{
		sites = null;
	}
	
	public Location getLocations(String seq)
	{
		List<Location> total = new ArrayList<Location>();
		if (sites.get(seq) == null)
			return null;
		
		for (Iterator<Site> itr=sites.get(seq).iterator(); itr.hasNext();)
		{
			total.add(itr.next().location);	
		}
		return LocationTools.union(total);
	}
	
	public long getTotalLength()
	{
		return totalLength;
	}

	public void setTotalLength(long totalLength)
	{
		this.totalLength = totalLength;
	}
	
	/**
	 * Merge a given Motif to this motif (from same dataset)
	 * @param m : motif to merge
	 * @throws Exception 
	 */
	public void merge(Motif m)
	{
		for (Iterator<String> itr=m.getSites().keySet().iterator(); itr.hasNext();)
		{
			String seq = itr.next();
			List<Site> list = sites.get(seq);
			List<Site> list2 = m.getSites().get(seq);
			if (list != null)
			{
				List<Site> temp = new ArrayList<Site>(list2);
				list.addAll(temp);
			}
			else
				list = new Vector<Site>(list2);
			sites.put(seq, list);	
		}
		
		totalLength += m.getTotalLength();	
	}
	
	/**
	 * intersect a given Motif with this motif (from same dataset)
	 * @param m : motif to intersect with
	 * @throws BioException 
	 * @throws NoSuchElementException 
	 * @throws IOException 
	 */
	public Motif intersect(Motif m) throws Exception
	{	
		Motif finalMotif = new Motif(finder, dataset, "", 0);
		Map<String, Site> sequences = TompaDataset.getSequences(dataset);
				
		for (Iterator<String> itr=getSites().keySet().iterator(); itr.hasNext();)
		{
			String seqName = itr.next();			
			Site seq = sequences.get(seqName);
			Location loc1 = getLocations(seqName);
			Location loc2 = m.getLocations(seqName);
			
			if (loc2 != null)
			{
				Location intersectLoc = LocationTools.intersection(loc1, loc2);
				List<Site> finalSites = new ArrayList<Site>();
				for (Iterator<Location> itr2=intersectLoc.blockIterator(); itr2.hasNext();)
				{
					Location loc = itr2.next();
					String dna = seq.getDNA(loc.getMin(), loc.getMax());
					finalSites.add(new Site(dna, seqName, loc.getMin(), loc.getMax(), "+")); ///////////////
				}
				finalMotif.getSites().put(seqName, finalSites);	
			}
		}
		
		return finalMotif;
	}
	
	public MotifFinder getFinder()
	{
		return finder;
	}

	public void setFinder(MotifFinder finder)
	{
		this.finder = finder;
	}
	public void setDataset(Dataset dataset)
	{
		this.dataset = dataset;
	}

	public int getAllLocationsCoverageCount()
	{
		int count = 0;
		
		for (Iterator<String> itr=sites.keySet().iterator(); itr.hasNext();)
			count+= LocationTools.coverage(getLocations(itr.next()));
		
		return count;
	}

	public Dataset getDataset()
	{
		return dataset;
	}
	public void setComment(String comment)
	{
		this.comment = comment;
	}
	public String getComment()
	{
		return comment;
	}
	
	public Map<String, List<Site>> getSites()
	{
		return sites;
	}

	public void setSites(Map<String, List<Site>> sites)
	{
		this.sites = sites;
	}
	
	public int getCorrect()
	{
		return correct;
	}
	
	public void setCorrect(int i)
	{
		correct = i;
	}
	/**
	 * Returns a copy of this Motif but with sites in the forward strand only
	 * @return
	 */
	public Motif getForwardOnly()
	{
		Motif m = new Motif(finder, dataset, comment, totalLength);
		
		for (Iterator<String> itr=getSites().keySet().iterator(); itr.hasNext();)
		{
			String seqName = itr.next();
			List<Site> list = new ArrayList<Site>(getSites().get(seqName).size());
			
			for (Iterator<Site> itr2=getSites().get(seqName).iterator(); itr2.hasNext();)
			{
				Site seq = itr2.next();
				if (seq.type.equals("+"))
					list.add(seq);
			}
			m.getSites().put(seqName, list);
		}
		return m;
	}
	
	@SuppressWarnings("serial")
	public static class Site extends StrandedFeature.Template
	{
		public String dna;
		
		public Site(){}
		
		public Site(String dna, String source, int start, int end, String direction) throws IllegalSymbolException
		{			
			this.dna = dna;
			strand = direction.equals("-")?StrandedFeature.NEGATIVE:StrandedFeature.POSITIVE;
			////template.type;
			this.source = source;
			location = new RangeLocation(start, end);
			annotation = Annotation.EMPTY_ANNOTATION;
			
		}		
		
		// light site
		public Site(String source, Location loc, Strand strand) throws IllegalSymbolException
		{			
			this.source = source;
			this.location = loc;
			this.strand = strand;
			
		}
		
		public Site getSubSite(int start, int end)
		{
			Site site = new Site();
			site.dna = getDNA(start, end);
			site.strand = strand;
			site.source = source;
			site.location = new RangeLocation(start, end);
			site.annotation = annotation; 
			return site;
		}
		
		// start and end are locations on the original DNA sequence in -ve's
		public String getDNA(int start, int end)
		{
			int a = start-location.getMin();
			int b = end-location.getMin();			
			return dna.substring(a, b+1);
		}
		
		public String getStrandString()
		{
			return strand.equals(StrandedFeature.NEGATIVE)?"-":"+";
		}
		
		public Site getIntersection(Site s2) throws IndexOutOfBoundsException, ChangeVetoException, BioException
		{
			Site intersect = null;
			
			if (LocationTools.overlaps(location, s2.location))
			{
				Location loc = LocationTools.intersection(location, s2.location);		
				intersect = getSubSite(loc.getMin(), loc.getMax());
			}			
			return intersect;
		}
		
		public Site getLocationIntersection(Site s2) throws IndexOutOfBoundsException, ChangeVetoException, BioException
		{
			Site intersect = null;
			
			if (LocationTools.overlaps(location, s2.location))
			{
				Location loc = LocationTools.intersection(location, s2.location);		
				intersect = new Site();
				intersect.location = loc;
				intersect.source = source;
				intersect.strand = strand;
			}			
			return intersect;
		}
		
		public double getWeightedCoverage(double weight)
		{
			return ((double)LocationTools.coverage(location))*weight;
		}
	}
	
	public int countIntersect(Motif m)
	{
		int count = 0;
		
		for (Iterator<String> itr=m.getSites().keySet().iterator(); itr.hasNext();)
		{
			String seq = itr.next();
			Location loc1 = getLocations(seq);
			Location loc2 = m.getLocations(seq);
			if (loc1 != null)
				count+= LocationTools.coverage(LocationTools.intersection(loc1, loc2));			
		}
				
		return count;
	}
	
	public double sim(Motif m)
	{
		double total = 0;
		int union = 0;
		int intersection = 0;
		Set<String> seqNames = new HashSet<String>(getSites().keySet());
		seqNames.addAll(m.getSites().keySet());
		
		for (Iterator<String> itr=seqNames.iterator(); itr.hasNext();)
		{
			String seqName = itr.next();
			Location loc1 = getLocations(seqName);
			Location loc2 = m.getLocations(seqName);
			
			if (loc1!=null && loc2!=null)
			{
				union += LocationTools.coverage(LocationTools.union(loc1, loc2));
				intersection += LocationTools.coverage(LocationTools.intersection(loc1, loc2));			
			}
			else if (loc1==null)
				union += LocationTools.coverage(loc2);
			else if (loc2==null)
				union += LocationTools.coverage(loc1);
		}
	
		if (union != 0)
			total += ((double)intersection/(double)union);
		
		return total;
	}
	
	@SuppressWarnings("unchecked")
	public void removeRedundant()
	{
		try
		{
			Map<String, Site> sequences = TompaDataset.getSequences(dataset); 
						
			for (Iterator<String> sitr=getSites().keySet().iterator(); sitr.hasNext();)
			{
				String seqName = sitr.next();						
				Location totalLocations = getLocations(seqName);//new CompoundRichLocation(sites);
				Site seq = sequences.get(seqName);
					
				List<Site> finalSites = new ArrayList<Site>();
				for (Iterator<Location> itr=totalLocations.blockIterator(); itr.hasNext();)
				{
					Location loc = itr.next();
					String dna = seq.getDNA(loc.getMin(), loc.getMax());
					finalSites.add(new Site(dna, seqName, loc.getMin(), loc.getMax(), "+")); ///////////////
				}
				getSites().put(seqName, finalSites);
			}
		}
		catch (Exception e)
		{e.printStackTrace();}
	}
	
	@SuppressWarnings("unchecked")
	public void removeRedundantFast()
	{
		try
		{			
			for (Iterator<String> sitr=getSites().keySet().iterator(); sitr.hasNext();)
			{
				String seqName = sitr.next();						
				Location totalLocations = getLocations(seqName);//new CompoundRichLocation(sites);
					
				List<Site> finalSites = new ArrayList<Site>();
				for (Iterator<Location> itr=totalLocations.blockIterator(); itr.hasNext();)
				{
					Location loc = itr.next();					
					finalSites.add(new Site(seqName, loc, StrandedFeature.POSITIVE)); ///////////////
				}
				getSites().put(seqName, finalSites);
			}
		}
		catch (Exception e)
		{e.printStackTrace();}
	}
	
	public String toString()
	{
		StringBuffer temp = new StringBuffer();
		temp.append(finder.name()+", ");
		int count = 0;
		
		for (Iterator<String> itr=sites.keySet().iterator(); itr.hasNext();)
		{
			String seq = itr.next();
			temp.append("site"+seq+"="+sites.get(seq).size()+", ");		
			count+= sites.get(seq).size();
		}
		temp.append("Total="+count);
		return temp.toString();
	}
}
