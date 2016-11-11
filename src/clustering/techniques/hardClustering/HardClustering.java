package clustering.techniques.hardClustering;

import java.util.*;

import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.tools.MotifTools;

import clustering.util.Util;

import clustering.dataTypes.*;
import clustering.distances.Distance;


public abstract class HardClustering 
{
	Distance distance;				// distance measure
	int k; 							// num of clusters
	int m;							// num of patterns in the data set
	public HardClusteringResult result;		// clustering results
	
	
		
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
		
	public abstract void partition();
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public HardClusteringResult partition(Vector<Pattern> data, Distance distance,int numOfClusters)
	{
		return partition(data, distance, numOfClusters, null);
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public HardClusteringResult partition(Vector<Pattern> data, Distance distance,int numOfClusters, int clusterSizes[])
	{		
		k = numOfClusters;			
		m = data.size();
		result = new HardClusteringResult(data, m, k);				
		this.distance = distance;		

		getInitialClustersOfMotifs(); //partitionAtRandom(clusterSizes); 
		
		long startTime = Calendar.getInstance().getTimeInMillis();		
		
		partition();
		long endTime = Calendar.getInstance().getTimeInMillis();
		
				
		result.runningTime = endTime - startTime ;
			
		
		return result;
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	/*
	public boolean calculateClustersCenters()
	{
		boolean flag = false;
		double n;
		Vector<Double> sum;
		
				
		for (int i=0; i<result.membership.length; i++) // for every row (cluster)
		{
			//n = 0;
			n = result.getNumOfPatternsInCluster(i);
			sum = new Vector<Double>(result.data.get(i).properties.size());
			
			for (int j=0; j<result.data.get(i).properties.size(); j++)
				sum.add(j, new Double(0));
			
			// sum member vectors
			for (int j=0; j<result.membership[i].length; j++) // for all member patterns
				if (!Util.equal(result.membership[i][j], Util.ZERO))
				{
					//n++;
					Vector<Double> temp =  new Vector<Double>(result.data.get(j).properties);
					for (int k=0; k<temp.size(); k++ )
					{
						Double d = new Double(sum.get(k)+temp.get(k));
						sum.set(k, d);
					}
				}
			
			
			for (int k=0; k<sum.size(); k++ )
			{
				if (n!=0)
					sum.set(k,sum.get(k)/n);
				else
					sum.set(k,new Double(0));
				if (!Util.equal(sum.get(k),result.clustersCenters.get(i).get(k)))
					flag = true;				
			}
			
			result.clustersCenters.put(i, new Vector<Double>(sum));
		}		
		
		return flag;
	}
*/
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
/*	
	public boolean calculateClustersCenters(int i, int j, double sign)
	{
		boolean flag = false;
		double n;
		Vector<Double> m, x;
		
		n = result.getNumOfPatternsInCluster(i);
		m = new Vector<Double>(result.clustersCenters.get(i));
		x = result.data.get(j).properties;
	
		
		for (int k=0; k<m.size(); k++ )
		{
			if (n!=0)
				m.set(k,m.get(k)+sign*(x.get(k)-m.get(k))/n);
			else
				m.set(k,new Double(0));
			
			
			if (!Util.equal((Double)m.get(k),(Double)result.clustersCenters.get(i).get(k)))
				flag = true;				
			
		}
		
		result.clustersCenters.put(i, new Vector<Double>(m));
		
		return flag;
	}
*/

	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public void partitionAtRandom(int n[])
	{		
		int count = 0;
			
		if (n==null)
			n = getRandomSizes(k,m);	
		
		for (int i=0; i<k; i++)
		{			
			for (int j=0; j<n[i]; j++)
			{
				result.clustersCenters.put(i, new Vector<Double>(result.data.get(count).properties));
				result.membership[i][count++] = 1.0;
			}
		}		
	}	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static int[] getRandomSizes(int k, int m)
	{
		int n[] = new int[k];
		int count = m-k;
			
		
		for (int i=0; i<k-1; i++)
		{
			if (count==0) break;
			
			n[i]= Util.generateRandomNumbers(1,count)[0];			
			count -= n[i];
		}
		
		n[k-1] = count;
		
		for (int i=0; i<k; i++)
		{
			n[i]++;
			System.out.println("\t cluster num"+i+" has: " + n[i]+" elements.");
		}
		
		return n;
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public Double calculateObjectiveFunction()
	{
		Double objFunc = new Double(0);
		
		for (int i=0; i<k; i++)
			for (int j=0; j<m; j++)
				objFunc+= (double) result.membership[i][j] * ((Double)distance.getDistanceSquare(result.data.get(j).properties, result.clustersCenters.get(i))).doubleValue();
		
		return objFunc;
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public Double calculateObjectiveFunctionWeight()
	{
		Double objFunc = new Double(0);
		
		for (int i=0; i<k; i++)
		{
			Vector<Motif> cluster = getCluster(i);
			Double d = MotifTools.sim(cluster);  // or use any weight function
			System.out.println("Weight(Cluster "+i+", numOfMotifs "+cluster.size()+")="+d);
			if (!d.equals(Double.NaN))
				objFunc+= d;
			
		}
		System.out.println("Obj Fun="+objFunc);
		return objFunc;
	}
	
/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public Double calculateObjectiveFunctionWeightRatio()
	{
		Double objFunc = new Double(0);
		
		for (int i=0; i<k; i++)
		{
			Vector<Motif> cluster = getCluster(i);
			Vector<Motif> rest = new Vector<Motif>();
			for (int j=0; j<k; j++)
				if (j!=i)	rest.addAll(getCluster(j));
			
			Double a = MotifTools.sim(cluster);  	// or use any weight function
			Double b = MotifTools.sim(rest);		// or use any weight function
			Double d = a / b;
			System.out.println("Ratio(Cluster "+i+", numOfMotifs "+cluster.size()+")="+d);
			if (!d.equals(Double.NaN))
				objFunc+= d;
			
		}
		System.out.println("Obj Fun="+objFunc);
		return objFunc;
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/		
	
	public void getInitialClustersOfMotifs()
	{
		Map<MotifFinder, Boolean> finder = getFinders(result.data);
		
		int j=0;
		for (int i=0; i<k; i++)
		{
			while(j < m)	// while more patterns
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
			if (j==m) // if num of clusters more than num of finders. get any center
			{
				result.clustersCenters.put(i, new Vector(result.data.get(0).properties));
				result.membership[i][0] = 1.0;
				
			}
		}
		
		// not imp for k-means, other methods need a pattern to belong to some cluster
		for (int i=0; i<k; i++)
			for (j=0; j<m; j++)
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
	
	public Vector<Motif> getCluster(int i)
	{
		Vector<Motif> cluster = new Vector<Motif>();
		for (int j=0; j<m; j++)
		{
			if (!Util.equalZero(result.membership[i][j]))
				cluster.add((Motif)result.data.get(j).properties.firstElement());		
		}
		
		return cluster;
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
}
