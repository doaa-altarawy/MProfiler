
package org.biojava.bio.tools;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.constants.Dataset;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;


/**
 * @author Doaa Altarawy
 *
 */
public class MotifEvaluation
{
	public static enum MEAURES {
		nTP, nTN, nFP, nFN, nSn, 
		nPPV, nPC, nCC,
		sTP, sFP, sFN, 
		sSn, sPPV, sASP};
	
	String name = "";
	
	long seqLength[];
	long totalLen;
	
	Motif target;
	Motif predicted;
										
	public int nTP, nTN, nFP, nFN;
	public double nSn;					// Sensitivity
	public double nPPV;					// positive predictive value(specificity)
	public double nPC;					// performance coefficient
	public double nCC;					// Correlation coefficient
	
	public int sTP, sFP, sFN;
	public double sSn;					// Sensitivity
	public double sPPV;					// sequence positive predictive value(specificity)
	public double sASP;					// (site level) average site performance
	
	
	public MotifEvaluation(String name)
	{
		this.name = name;		
		target = new Motif(Dataset.unknown);
	}
	
	public MotifEvaluation(String name, Motif target, Motif predicted)
	{
		this.name = name;
		this.totalLen = target.getTotalLength();
		this.target = target;
		this.predicted = predicted;	
		
		calculate();
	}
		
	
	/**
	 * Calculate the performance measures
	 */
	private void calculate()
	{
		Map<String,Location> intersection = new HashMap<String, Location>(target.getSites().keySet().size());
		
		Location loc1 = null, loc2 = null;
		
		// for each sequence in target		
		for (Iterator<String> strItr=target.getSites().keySet().iterator();strItr.hasNext();)
		{
			String n = strItr.next();
			List<Site> targets = target.getSites().get(n);
			List<Site> predicts = predicted.getSites().get(n);
			Location intersect = null;
			
			for (Iterator<Site> seqItr=targets.iterator(); seqItr.hasNext()&& predicts!=null;)
			{
				loc1 = seqItr.next().location;
				
				for (Iterator<Site> seqItr2=predicts.iterator();seqItr2.hasNext();)
				{
					loc2 = seqItr2.next().location;
					if (LocationTools.overlaps(loc1, loc2))
					{						
						if (intersect==null)
							intersect = LocationTools.intersection(loc1, loc2);
						else
							intersect = LocationTools.union(intersect,(LocationTools.intersection(loc1, loc2)));
					}				
				}
			}
			intersection.put(n, intersect);						
		}
		
		// for each sequence in target		
		for (Iterator<String> strItr=target.getSites().keySet().iterator();strItr.hasNext();)
		{
			String n = strItr.next();
			Location found = null;
			List<Site> targets = target.getSites().get(n);
			
						
			for (Iterator<Site> seqItr=targets.iterator(); seqItr.hasNext();)
			{
				loc1 = seqItr.next().location;
				if ((intersection.get(n)!=null) && LocationTools.overlaps(intersection.get(n), loc1))
				{
					Location temp = LocationTools.intersection(intersection.get(n), loc1);
					if (found ==null)
						found = temp;
					else 
						found = LocationTools.union(found, temp);
				
					// consider the overlap as found location if overlapping is more than 25%
					if ((double)LocationTools.coverage(temp)>= (0.25 * (double)LocationTools.coverage(loc1)))
						sTP++;
					else
						sFN++;
				}
				else
					sFN++;
			}
			if (found!=null)
			{
				nTP += LocationTools.coverage(found);
				nFN += LocationTools.coverage(LocationTools.subtract(target.getLocations(n), found));
			}
			else
				nFN += LocationTools.coverage(target.getLocations(n));
		}
		
		// for each sequence in predicted		
		for (Iterator<String> strItr=predicted.getSites().keySet().iterator();strItr.hasNext();)
		{
			String n = strItr.next();
			Location found = null;
			List<Site> predicts = predicted.getSites().get(n);
			
						
			for (Iterator<Site> seqItr=predicts.iterator(); seqItr.hasNext();)
			{
				loc1 = seqItr.next().location;
				Location temp = null;
				if ((intersection.get(n)!=null) &&  LocationTools.overlaps(intersection.get(n), loc1))
				{
					temp = LocationTools.intersection(intersection.get(n), loc1);
					if (found ==null)
						found = temp;
					else 
						found = LocationTools.union(found, temp);				
				}
				if ((temp==null) || ((double)LocationTools.coverage(temp)< (0.25 * (double)LocationTools.coverage(loc1))))
					sFP++;
			}
			if (found!=null)
				nFP += LocationTools.coverage(LocationTools.subtract(predicted.getLocations(n),found));
			else
				nFP += LocationTools.coverage(predicted.getLocations(n));
				
		}
		
		
		nTN += totalLen - nTP - nFP - nFN;
		
		calculateStatistics();
		
	}	
	
	private void calculateStatistics()
	{
		double x, y, z;
		if ((x = nTP + nFN)!= 0 )
			nSn = (double)nTP / x;
		else nSn += 0;
		
		if ((y = nTP + nFP)!= 0)
			nPPV = (double)nTP / y;
		else nPPV += 0;
		
		if ((x + nFP)!=0)
			nPC = (double)nTP / (double)(x + nFP);
		else nPC += 0;
		
		z = x * y * (nTN + nFP) * (nTN +nFN);
		
		if (z==0)
			nCC += 0; ///////////////////
		else
			nCC = (double)(nTP*nTN-nFN*nFP)/ Math.sqrt(z);
		
		
		if ((x = sTP + sFN)!= 0 )
			sSn = (double)sTP / x;
		else sSn += 0;
		
		if ((y = sTP + sFP)!= 0)
			sPPV = (double)sTP / y;
		else sPPV += 0;
		
		sASP = (sSn + sPPV) / 2.0;
	}
	
	public String toString()
	{
		StringBuffer buffer = new StringBuffer("Data set: "+ target.getDataset()+"\n");
		
		buffer.append("nTP= "+nTP+"\n");
		buffer.append("nFP= "+nFP+"\n");
		buffer.append("nFN= "+nFN+"\n");
		buffer.append("nTN= "+nTN+"\n");
		buffer.append("nSn= "+nSn+"\n");
		buffer.append("nPPV= "+nPPV+"\n");
		buffer.append("nPC= "+nPC+"\n");
		buffer.append("nCC= "+nCC+"\n");
		buffer.append("sTP= "+sTP+"\n");
		buffer.append("sFP= "+sFP+"\n");
		buffer.append("sFN= "+sFN+"\n");
		buffer.append("sSn= "+sSn+"\n");
		buffer.append("sPPV= "+sPPV+"\n");
		buffer.append("sASP= "+sASP+"\n");
		buffer.append("---------------"+"\n");
		buffer.append("---------------"+"\n");
		return buffer.toString();
	}
	
	/**
	 * Used to merge results from different Datasets
	 * @param m
	 */
	public void mergeDatasetEvaluation(MotifEvaluation m)
	{
		nTP += m.nTP;
		nTN += m.nTN;
		nFP += m.nFP;
		nFN += m.nFN;
		
		sTP += m.sTP;
		sFP += m.sFP;
		sFN += m.sFN;
		
		calculateStatistics();
	}
	
	public String getName()
	{
		return name;
	}
	
	public double getMeasure(MEAURES m)
	{
		switch (m)	
		{
			case nTP: return nTP;
			case nFN: return nFN;
			case nFP: return nFP;
			case nTN: return nTN;
			case nCC: return nCC;
			case nPC: return nPC;
			case nSn: return nSn;
			case nPPV: return nPPV;
			case sTP: return sTP;
			case sFN: return sFN;
			case sFP: return sFP;
			case sSn: return sSn;
			case sPPV: return sPPV;
			case sASP: return sASP;			
			default: return 0;
		}
	}	
}
