package clustering.dataTypes;

import java.util.*;

import org.biojava.bio.motif.Motif;


public class HierarchialClusteringResult extends ClusteringResult 
{
	public int minNumOfClusters;
	public Map<Integer,Vector<Double>> dissimilarityMatrix;		/* 	Symmetric matrix thus lower part is only stored, 
																	diagonal is 0, null entry means deleted	*/
	public Vector<Vector<Integer>> patternAssigHistory;		// History of pattern assigned to clusters
	
	public Tree tree;
		
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	
	/**	 
	 * where num of clusters in unknown
	 */
	public HierarchialClusteringResult( Vector<Pattern> data, int numOfPatterns)
	{
		super(data, numOfPatterns);
		
		numOfClusters = numOfPatterns;		// initial num of clusters
		
		this.minNumOfClusters = 1;
		tree = new Tree();
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	/**	 
	 * where num of clusters in known
	 */
	public HierarchialClusteringResult(Vector<Pattern> data,  int numOfPatterns,int minNumOfClusters )
	{
		super(data, numOfPatterns);
		
		this.numOfClusters = numOfPatterns;		// initial num of clusters
				
		this.minNumOfClusters = minNumOfClusters;
		tree = new Tree();
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public int getNumOfPatternsInCluster(int cluster)
	{
		int sum = 0;
		
		for (int i=0; i<patternAssigHistory.lastElement().size(); i++)
			if (patternAssigHistory.lastElement().get(i)==cluster)
				sum++;
		
		return sum;
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	// get last level cluster
	public Collection<Vector<Motif>> getClusters()
	{
		Map<Integer, Vector<Motif>> clusters = new HashMap<Integer,Vector<Motif>>();
		for (int i=0; i<patternAssigHistory.lastElement().size(); i++)
		{
			int cluster = patternAssigHistory.lastElement().get(i);
			if (clusters.containsKey(cluster))
				clusters.get(cluster).add((Motif)data.get(i).properties.firstElement());
			else
			{
				Vector<Motif> v = new Vector<Motif>();
				v.add((Motif)data.get(i).properties.firstElement());
				clusters.put(cluster, v);
			}
		}
		return clusters.values();
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	public Vector<Motif> getCluster(int n)
	{
		Vector<Motif> cluster = new Vector<Motif>();
		for (int i=0; i<patternAssigHistory.lastElement().size(); i++)
		{
			int clusterIndex = patternAssigHistory.lastElement().get(i);
			if (clusterIndex==n)
				cluster.add((Motif)data.get(i).properties.firstElement());
			
		}
		return cluster;
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public Collection<Vector<Motif>> getClusters(int level)
	{
		Map<Integer, Vector<Motif>> clusters = new HashMap<Integer,Vector<Motif>>();
		for (int i=0; i<patternAssigHistory.get(level).size(); i++)
		{
			int cluster = patternAssigHistory.get(level).get(i);
			if (clusters.containsKey(cluster))
				clusters.get(cluster).add((Motif)data.get(i).properties.firstElement());
			else
			{
				Vector<Motif> v = new Vector<Motif>();
				v.add((Motif)data.get(i).properties.firstElement());
				clusters.put(cluster, v);
			}
		}
		return clusters.values();
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
		
	
}
