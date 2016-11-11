package org.biojava.bio.motif;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.Vector;



import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.tools.TompaDataset;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.CompoundRichLocation;

/**
 * Implementation of MotifVoter: an Ensemble method for motif finding
 * @author Doaa Altarawy
 *
 */
public class MotifVoter
{
	// Input
	Dataset dataset;		
	Vector<Motif> inputMotifs;				// vector of all given motifs
	
	// Processing data structure	
	Map<Motif, Map<Motif, Double>> sim;
	Map<Motif, Vector<Motif>> orderedSim;
	
	// Results
	List<Motif> accepted;
	List<Motif> finalSites;
	Motif finalMotif;
		
	
	/* Constructors */	
	public MotifVoter(Dataset dataset)
	{
		this.dataset = dataset;		
	}
	
	public MotifVoter(Dataset dataset, Vector<Motif> motifs)
	{
		this.dataset = dataset;
		this.inputMotifs = motifs;
	}	
	/*----------------------------------------------------------------------------*/
	
	public Motif findMotif() throws IndexOutOfBoundsException, ChangeVetoException, BioException, NoSuchElementException, IOException
	{
		computeSortedSimilarity();
		
		List<Motif> maxSet = null;
		int maxNumOfFinders = 0;
		double maxWeight = -1.0;
		
		List<Motif> x, p_x;
		int numOfFinders = 0;
		double weight = -1.0;
		int MAX_NUM_OF_FINDERS = getNumOfFinders(inputMotifs);
		
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
		{
			x = new ArrayList<Motif>(inputMotifs.size());
			p_x = new ArrayList<Motif>(inputMotifs);
			Motif z = itr.next(); 
			x.add(z);
			p_x.remove(z);
			System.out.println("Finder z:"+z.getFinder().name());
			
			Map<Motif, Double> map = sim.get(z);
						
			for (Iterator<Motif> itr2 = orderedSim.get(z).iterator(); itr2.hasNext();)	
			{
				Motif p = itr2.next();
				if (z.equals(p)) //skip same motif
					continue;
				x.add(p);
				p_x.remove(p);				
				
				weight = weightRatio(x, p_x);
				numOfFinders = getNumOfFinders(x);
				
				//System.out.println("p="+map.get(p)+", weight="+weight+ ", Finders="+numOfFinders);
					
				//if ((weight>=maxWeight) && (numOfFinders>2))
				//if ((weight>=maxWeight) &&  numOfFinders==MAX_NUM_OF_FINDERS) 
				//if ((weight>=maxWeight) && (numOfFinders>maxNumOfFinders) && (numOfFinders>3))
				if ((weight>=maxWeight) && (numOfFinders>maxNumOfFinders))
				{
					System.out.println("\n---------\n Max changed to:"+weight);
					System.out.println("Max Finders to:"+numOfFinders);	
					maxWeight = weight;
					maxNumOfFinders = numOfFinders;
					maxSet = new Vector<Motif>(x);
				}
			}
		}
			
		accepted = maxSet;
		finalMotif = extractSites(accepted);
		
		return finalMotif;
	}
	
	/*-----------------------------------------------------------------------------*/
	
	private void computeSortedSimilarity()
	{
		Sortable<Motif>[][] similarities = new Sortable[inputMotifs.size()][inputMotifs.size()];
		sim = new HashMap<Motif, Map<Motif,Double>>(inputMotifs.size());
		orderedSim = new HashMap<Motif, Vector<Motif>>(inputMotifs.size());
		
		for (int i=0; i<inputMotifs.size(); i++)
		{			
			for (int j=0; j<inputMotifs.size(); j++)
			{				
				double d = sim(inputMotifs.get(i), inputMotifs.get(j));
				similarities[i][j] = new Sortable<Motif>(inputMotifs.get(j), d);
			}
			Arrays.sort(similarities[i]);	
			Map<Motif, Double> map = new HashMap<Motif, Double>(inputMotifs.size()-1);
			Vector<Motif> v = new Vector<Motif>(inputMotifs.size()-1);
			for (int k=0; k<similarities[i].length; k++)
			{
				map.put(similarities[i][k].obj, similarities[i][k].value);
				v.add(similarities[i][k].obj);
			}
			sim.put(inputMotifs.get(i), map);	
			orderedSim.put(inputMotifs.get(i), v);
		}		
		
	}
	
	/*-----------------------------------------------------------------------------*/

	int NANcount = 0;
	int infCount = 0;
	private double weightRatio(List<Motif> x, List<Motif> p_x)
	{
		double a = weight(x);
		double b = weight(p_x);
		
		if (Double.isInfinite(b))	// doesn't happen
			{b = 0; infCount++;}
		if (b==0)
			{NANcount++; return Double.NaN;}
		//System.out.println("w(X)="+a+" , w(P-X)="+b);
		return a / b;
	}
	
	/*-----------------------------------------------------------------------------*/
	
	private int getNumOfFinders(List<Motif> motifs)
	{		
		Set<MotifFinder> finders = new HashSet<MotifFinder>();
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			finders.add(m.finder);			// This is set, elements doesn't repeat				
		}
		
		return finders.size();
	}
	/*-----------------------------------------------------------------------------*/
	public Motif getMergedSites()
	{
		Motif finalMotif = new Motif(MotifFinder.MotifVoter);
		
		for (Iterator<Motif> itr=accepted.iterator(); itr.hasNext();)
		{
			finalMotif.merge(itr.next());
		}
		return finalMotif;
	}
	/*-----------------------------------------------------------------------------*/
	/**
	 * Todo : need recoding, use removeRedundant method (compoundLocation) in the extraxtion step
	 */
	// get the set of finalSites then get the final motif
	public Motif extractSites(List<Motif> motifs) throws IndexOutOfBoundsException, ChangeVetoException, BioException, NoSuchElementException, IOException
	{
		Motif finalMotif = new Motif(MotifFinder.MProfiler);
		finalMotif.setDataset(dataset);
		Motif covered = new Motif(MotifFinder.MProfiler);
		covered.setDataset(dataset);
		List<Motif> temp = new ArrayList<Motif>(motifs);
		Motif motifWithMaxConf = getMaxConfidenceScore(motifs);
		
		// get intersection covering
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			temp.remove(m);
			for (Iterator<Motif> itr2=temp.iterator(); itr2.hasNext();)
			{
				Motif m2 = itr2.next();				
				for (Iterator<String> n=m.sites.keySet().iterator(); n.hasNext();)
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
							Site intersectSeq = getIntersection(seq1, seq2);
							
							if (intersectSeq!=null)
								coveredList.add(intersectSeq);
						}						
					}
					if (!coveredList.isEmpty())
						covered.getSites().put(seqName, coveredList);					
				}
			}
		}
		
		finalMotif.merge(motifWithMaxConf);
		finalMotif.merge(covered);
		
		removeRedundant(finalMotif);
		
		return finalMotif;
	}
	/*-----------------------------------------------------------------------------*/
	
	@SuppressWarnings("unchecked")
	private Motif removeRedundant(Motif motif) throws NoSuchElementException, BioException, IOException
	{
		SequenceIterator seqItr = TompaDataset.getDatasetSequences(dataset);
		Map<String, Site> sequences = new HashMap<String, Site>();
		while (seqItr.hasNext())
		{
			Sequence seq = seqItr.nextSequence();
			String source = seq.getName().substring(seq.getName().indexOf('_')+1);
			sequences.put(source, new Site(seq.seqString(), source, -seq.length(), -1, "+"));
		}
		
		for (Iterator<String> sitr=motif.getSites().keySet().iterator(); sitr.hasNext();)
		{
			String seqName = sitr.next();						
			Location totalLocations = motif.getLocations(seqName);//new CompoundRichLocation(sites);
			Site seq = sequences.get(seqName);
				
			List<Site> finalSites = new ArrayList<Site>();
			for (Iterator<Location> itr=totalLocations.blockIterator(); itr.hasNext();)
			{
				Location loc = itr.next();
				String dna = seq.getDNA(loc.getMin(), loc.getMax());
				finalSites.add(new Site(dna, seqName, loc.getMin(), loc.getMax(), "+")); ///////////////
			}
			motif.getSites().put(seqName, finalSites);
		}
		return motif;
	}
	/*-----------------------------------------------------------------------------*/

	private Site getIntersection(Site s1, Site s2) throws IndexOutOfBoundsException, ChangeVetoException, BioException
	{
		Site intersect = null;
		
		if (LocationTools.overlaps(s1.location, s2.location))
		{
			Location loc = LocationTools.intersection(s1.location, s2.location);		
			intersect = s1.getSubSite(loc.getMin(), loc.getMax());
		}
		
		return intersect;
	}
	
	/*-----------------------------------------------------------------------------*/
	/**
	 * Returns the motif with Max confidence score
	 */
	private Motif getMaxConfidenceScore(List<Motif> motifs)
	{
		Map<Motif, Double> scores = getConfidenceScore(motifs);
		Motif motif = null;
		double max = -1;
		for (Iterator<Motif> itr=scores.keySet().iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			double value = scores.get(m);
			System.out.println("Motif "+m.getFinder().name()+" score: "+value);
			if (value>max)
			{
				max = value;
				motif = m;
			}
		}
		
		return motif;
	}
	/*-----------------------------------------------------------------------------*/
	private Map<Motif, Double> getConfidenceScore(List<Motif> motifs)
	{
		Map<Motif, Double> scores = new HashMap<Motif, Double>(motifs.size());
				
		for (Iterator<Motif> itr= motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			double intersect = 0;
			double size = 0;
			for (Iterator<Motif> itr2=motifs.iterator(); itr2.hasNext();)
			{
				Motif m2 = itr2.next();
				if (m == m2) continue;
				for (Iterator<String> seqItr=m.getSites().keySet().iterator(); seqItr.hasNext();)
				{
					String seq = seqItr.next();
					Location loc1 = m.getLocations(seq);
					Location loc2 = m2.getLocations(seq);
					size += LocationTools.coverage(loc1);
					if (loc2!=null && LocationTools.overlaps(loc1, loc2))
					{
						intersect += LocationTools.coverage(LocationTools.intersection(loc1, loc2));
					}
				}
				
			}
			scores.put(m, intersect/size);
		}
		
		return scores;
	}
	/*-----------------------------------------------------------------------------*/
	
	private double weight(List<Motif> motifs)
	{
		double sum = 0;		
		double total = 0;
		double simm = sim(motifs);
		//List<Motif> list = new ArrayList<Motif>(motifs);
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			//list.remove(m);
			Map<Motif, Double> mSim = sim.get(m);
			for (Iterator<Motif> itr2=motifs.iterator(); itr2.hasNext();)
			{
				Motif temp = itr2.next();
				//if (m==temp) continue;
				sum += Math.pow(mSim.get(temp) - simm, 2);
			}
		}
		
		if (sum != 0)
			total = (simm / Math.sqrt(sum));//(Math.pow(motifs.size(),2) *simm / Math.sqrt(sum));
		return total;
	}
	
	/*-----------------------------------------------------------------------------*/
	
	private double sim(List<Motif> motifs)
	{
		double sum = 0;
		//List<Motif> list = new ArrayList<Motif>(motifs);
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			//list.remove(m);
			Map<Motif, Double> mSim = sim.get(m);
			for (Iterator<Motif> itr2=motifs.iterator(); itr2.hasNext();)
			{
				Motif temp = itr2.next();
				//if (m==temp) continue;
				sum += mSim.get(temp);
			}
		}
			
		double total = (sum / Math.pow(motifs.size(), 2));
		return total;
	}
	/*-----------------------------------------------------------------------------*/

	private double sim(Motif m1, Motif m2)
	{
		double total = 0;
		int union = 0;
		int intersection = 0;
		Set<String> seqNames = new HashSet<String>(m1.getSites().keySet());
		seqNames.addAll(m2.getSites().keySet());
		
		for (Iterator<String> itr=seqNames.iterator(); itr.hasNext();)
		{
			String seqName = itr.next();
			Location loc1 = m1.getLocations(seqName);
			Location loc2 = m2.getLocations(seqName);
			
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
	/*-----------------------------------------------------------------------------*/

	private class Sortable<T> implements Comparable<Sortable<T>>
	{
		T obj;
		Double value = 0.0;
		
		public Sortable(T m)
		{
			obj = m;			
		}
		
		public Sortable(T m, double val)
		{
			obj = m;
			this.value = val;
		}

		// reversed to sort descending
		public int compareTo(Sortable<T> o)
		{
			if (value>o.value) return -1;		// returns -ve if greater to sort descending
			else if (value<o.value) return 1;	// returns +ve if smaller to sort descending
			else return 0;
		}
	}	
	/*-----------------------------------------------------------------------------*/

	/* Setters and getters */
	public Dataset getDataset()
	{
		return dataset;
	}

	public void setDataset(Dataset dataset)
	{
		this.dataset = dataset;
	}

	public List<Motif> getAccepted()
	{
		return accepted;
	}

	public void setAccepted(List<Motif> accepted)
	{
		this.accepted = accepted;
	}

	public Motif getFinalMotif()
	{
		return finalMotif;
	}

		
	public List<Motif> getFinalSites()
	{
		return finalSites;
	}

}
