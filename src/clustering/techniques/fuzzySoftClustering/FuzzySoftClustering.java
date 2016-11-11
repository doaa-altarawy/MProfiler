package clustering.techniques.fuzzySoftClustering;

import java.util.*;

import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.FuzzyMotif;
import org.biojava.bio.motif.Motif;

import clustering.util.Util;
import clustering.dataTypes.FuzzySoftClusteringResult;
import clustering.dataTypes.Pattern;
import clustering.distances.Distance;

public abstract class FuzzySoftClustering
{
	Distance distance;				// distance measure
	public FuzzySoftClusteringResult result;		// clustering results
	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	protected abstract void calculatePseudoMembershipMatrix();	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public void partition()
	{
		boolean flag = false;
		
		// step 1
		//assignInitialClustersCenters();
		getInitialClustersOfMotifs();
		
		do
		{
			// step 2
			calculateSquaredDistanceMatrix();
			
			// step 3
			calculatePseudoMembershipMatrix();
			
			// step 4
			normalizeMembershipMatrix();
			
			// step 5
			flag = distance.calculateClustersCenters(result);	////////////	
				
			result.numOfIterations++;
			if (result.numOfIterations>Util.maxNumOfIterations)
			{
				System.out.println("Max # of itrations exceeded");
				flag=false;
			}
				
		}
		while (flag);
	}
	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
		
	
	
	public FuzzySoftClusteringResult partition(Vector<Pattern> data, Distance distance,int m, int k, Object T)
	{
		this.distance = distance;		
				
		result = new FuzzySoftClusteringResult(data, m, k, T);				
		
		if (m<=1)
		{
			System.out.println("Error: m must be > 1");
			return null;
		}
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		partition();
		long endTime = Calendar.getInstance().getTimeInMillis();
		
		result.runningTime = endTime - startTime ;
		
		return result;
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public void assignInitialClustersCenters()
	{
		// taken as first k patterns
		result.clustersCenters = new HashMap<Integer,Vector>(result.numOfClusters);
		for (int i=0; i<result.numOfClusters; i++)
		{
			FuzzyMotif m = new FuzzyMotif((Motif)result.data.get(i).properties.firstElement(), 1.0);
			Vector v = new Vector(1);
			v.add(m);
			result.clustersCenters.put(i, v);
		}
			
	}
	
	public void getInitialClustersOfMotifs()
	{
		Map<MotifFinder, Boolean> finder = getFinders(result.data);
		
		int j=0;
		for (int i=0; i<result.numOfClusters; i++)
		{
			while(j < result.numOfPatterns)	// while more patterns
			{
				MotifFinder f = ((Motif)result.data.get(j).properties.firstElement()).getFinder();
				if (!finder.get(f))
				{
					result.clustersCenters.put(i, new Vector(result.data.get(j).properties));
					result.membership[i][j] = 1.0;
					j++;
					finder.put(f, true);
					break;
				}
				j++;
			}
			if (j==result.numOfPatterns) // if num of clusters more than num of finders. get any center
			{
				result.clustersCenters.put(i, new Vector(result.data.get(0).properties));
				result.membership[i][0] = 1.0;
				
			}
		}
		
		// not imp for k-means, other methods need a pattern to belong to some cluster
		for (int i=0; i<result.numOfClusters; i++)
			for (j=0; j<result.numOfPatterns; j++)
			{
				if (result.membership[i][j] == null)
					result.membership[i][j] = 0.0;
			}
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/

	private Map<MotifFinder, Boolean> getFinders(Vector<Pattern> motifs)
	{		
		Map<MotifFinder, Boolean> finders = new HashMap<MotifFinder, Boolean>();
		
		for (Iterator<Pattern> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = (Motif)itr.next().properties.firstElement();
			finders.put(m.getFinder(), false);						
		}
		
		return finders;
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public boolean calculateClustersCenters()
	{
		boolean flag = false;		
		Vector<Double> sum;
		Double sumWi;
		
		for (int i=0; i<result.numOfClusters; i++) // for every cluster
		{
						
			sum = new Vector<Double>(result.data.get(i).properties.size());
			sumWi = new Double(0);
			
			for (int j=0; j<result.data.get(i).properties.size(); j++)
				sum.add(j, new Double(0));				
					
						
			// sum member vectors
			for (int j=0; j<result.numOfPatterns; j++) // for all patterns
			{	
				double d = Math.pow(result.membership[i][j],result.m);
				sumWi += d;
				
				Vector<Double> x =  new Vector<Double>(result.data.get(j).properties);
				for (int k=0; k<x.size(); k++) // for all properties
					sum.set(k, sum.get(k)+ x.get(k)* d);				
			}
									
			for (int k=0; k<sum.size(); k++)
			{
				if (!Util.equal(sumWi, Util.ZERO))
					sum.set(k, sum.get(k)/sumWi);
				else
				{
					System.out.println("sumWiiiiiiiiiiiiiiiiiii="+sumWi);					
					System.out.println("cluster="+i);
					sum.set(k, sum.get(k)/sumWi);
				}
				
				//if (!Util.equal(sum.get(k), result.clustersCenters.get(i).get(k)))
				//	flag = true;					
			}			
			result.clustersCenters.put(i, new Vector<Double>(sum));
		}		
		
		return flag;
	}

	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public void calculateSquaredDistanceMatrix()
	{		
		for (int i=0; i<result.numOfClusters; i++)
			for (int j=0; j<result.numOfPatterns; j++)
				result.distanceMatrix[i][j] = ((Double)distance.getDistanceSquare(result.clustersCenters.get(i), result.data.get(j).properties)).doubleValue();
		
	}	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public void normalizeMembershipMatrix()
	{
		Double sum;
		
		for (int i=0; i<result.numOfPatterns; i++) // for every pattern
		{
			sum = new Double(0);
			for (int j=0; j<result.numOfClusters; j++)
				sum += result.pseudoMembership[j][i];
			
						
			for (int j=0; j<result.numOfClusters; j++)
				result.membership[j][i] = result.pseudoMembership[j][i] / sum;
				
		}
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/	
}
