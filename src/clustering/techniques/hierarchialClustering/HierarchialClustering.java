package clustering.techniques.hierarchialClustering;

import java.util.*;

import clustering.util.Util;
import clustering.dataTypes.*;
import clustering.distances.Distance;

public abstract class HierarchialClustering
{
	Distance distance;								// distance measure
	public HierarchialClusteringResult result;		// clustering results
	
	
		
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	protected abstract double getai(int i, int j, int k);
	protected abstract double getaj(int i, int j, int k);
	protected abstract double getb(int i, int j, int k);
	protected abstract double getc(int i, int j, int k);	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public ClusteringResult partition( Vector<Pattern> data, Distance distance, int numOfPatterns, int minNumOfClusters) throws Exception
	{
		this.distance = distance;		
		
		result = new HierarchialClusteringResult(data, numOfPatterns, minNumOfClusters);
				
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		partition();
		long endTime = Calendar.getInstance().getTimeInMillis();
		
		result.runningTime = endTime - startTime ;
		
		return result;
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	
	public ClusteringResult partition( Vector<Pattern> data, Distance distance, int numOfPatterns) throws Exception
	{				
		return partition(data, distance, numOfPatterns,1);
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public ClusteringResult partition( Vector<Pattern> data, Distance distance) throws Exception
	{
		return partition(data, distance, data.size(),1);
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public void partition() throws Exception
	{		
		Tree treeList[] = new Tree[result.data.size()];
		Double min = null;
		int c1=0, c2=0;
		
		for (int i=0; i<result.data.size(); i++)
			treeList[i] = new Tree();
		
		// step 1: each pattern is assigned to a unique cluster
		initialize();
		
		if (result.numOfClusters<= result.minNumOfClusters) return;
		
		do
		{
			min = null; c1=0; c2=0;
			// step 2: find smallest entry in dissimilarityMatrix 
			for (int i=0; i<result.dissimilarityMatrix.size(); i++)
			{
				for (int j=0 ;j<result.dissimilarityMatrix.get(i).size(); j++)
				{
					Double temp = result.dissimilarityMatrix.get(i).get(j);
					if (temp==null || Util.equal(temp, Util.ZERO) ) // means a deleted cluster (merged)
						continue;		
									
					if(min== null || temp<=min)
					{
						min = temp;
						c1 = i;
						c2 = j;
					}					
				}				
			}// find smallest entry
			
			// store data in the tree			
			Tree temp;
			if (c1!=c2)
			{
				temp = new Tree(min, treeList[c2], treeList[c1]);				
				treeList[c1] = null;
				treeList[c2] = temp;
			}
			
			
			// merge the 2 clusters c1, c2 into c2 (c2 is always < c1)			
			Vector<Integer> oldClusters = new Vector<Integer>(result.patternAssigHistory.lastElement());
			for (int i=0; i<oldClusters.size(); i++)
				if(oldClusters.get(i)==c1)	oldClusters.set(i,c2);
			result.patternAssigHistory.add(oldClusters);
			
			// step 3: update dissimilarityMatrix, update c2, delete c1 distances			
			// a- update c2
			
			for (int k=0; k<result.dissimilarityMatrix.size(); k++)
			{				
				if (k == c2 || k == c1) 
					continue;
				if (result.dissimilarityMatrix.get(c2>k?c2:k).get(c2>k?k:c2)==null)
					continue;
				if (result.dissimilarityMatrix.get(c1>k?c1:k).get(c1>k?k:c1)==null)
					continue;
				
				double d1,d2,d3, newD;
								
				
				d1 = result.dissimilarityMatrix.get(c1>k?c1:k).get(c1>k?k:c1);					
				d2 = result.dissimilarityMatrix.get(c2>k?c2:k).get(c2>k?k:c2); 
				d3 = result.dissimilarityMatrix.get(c1).get(c2); 	
				
				newD = getai(c1,c2,k)*d1 + getaj(c1,c2,k)*d2 + getb(c1,c2,k)*d3 + getc(c1,c2,k)*Math.abs(d1-d2);
				
				
				result.dissimilarityMatrix.get(c2>k?c2:k).set(c2>k?k:c2,newD);
			}
			
			
			// b- delete c1
			for (int i=0; i<result.dissimilarityMatrix.size(); i++)
				if (i!=c1)
					result.dissimilarityMatrix.get(c1>i?c1:i).set(c1>i?i:c1,null);
			
			result.dissimilarityMatrix.get(c1).set(c2,new Double(0));	// distance bet c1 & c2 = 0 (merged clusters)
			
			
			result.numOfClusters--;
			result.numOfIterations++;
			
		}
		while (result.numOfClusters>result.minNumOfClusters);
		
		result.tree = treeList[0];
				
		
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	
	public void initialize()
	{
		result.patternAssigHistory = new Vector<Vector<Integer>>();
		int unsorted[] = Util.generateUnsortedList(result.numOfPatterns); // same effect in the algorithm if it sorted
				
		
		Vector<Integer> temp = new Vector<Integer>(result.numOfPatterns); 
		for (int i=0; i<result.numOfPatterns; i++)		// assign every pattern to a different cluster
			temp.add(unsorted[i]);					
		
		result.patternAssigHistory.add(temp);			// first trace (history) of clusters assignment
		
		// fill dissimilarityMatrix 
		result.dissimilarityMatrix = new HashMap<Integer,Vector<Double>>(result.numOfClusters);
		for (int i=0; i< result.numOfClusters; i++)
		{	
			result.dissimilarityMatrix.put(i,new Vector<Double>());
			for (int j=0; j<i; j++)
			{
				Double d = (Double)distance.getDistanceSquare(result.data.get(i).properties, result.data.get(j).properties);
				
				result.dissimilarityMatrix.get(i).add(d);
				
			}
		}	
		
	}	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
}
