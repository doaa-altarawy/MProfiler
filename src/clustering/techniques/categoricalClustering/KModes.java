package clustering.techniques.categoricalClustering;

import java.util.*;
import java.util.Map.Entry;

import clustering.util.Util;

import clustering.dataTypes.Pattern;
import clustering.dataTypes.categoricalData.CategoricalClusteringResult;

public class KModes extends CategoricalClusteringKMeans<Vector<Integer>, Integer> 
{

	@Override
	public void calculateClustersCenters() 
	{	
		int numOfAtrr = result.data.get(0).properties.size(); 
		Vector<Integer> freq;
		
		result.clustersCenters = new HashMap<Integer, Vector<Integer>>();
		result.frequences = new Vector<HashMap<Integer, Vector<Integer>>>();
				
		for (int i=0; i<result.numOfClusters; i++) // for every col (cluster)
		{
			result.clustersCenters.put(i,new Vector<Integer>());
			result.frequences.add(i,new HashMap<Integer, Vector<Integer>>());
			
			for (int j=0; j<numOfAtrr; j++)
			{				
				int numOfCatOfAttrJ = result.properties.get(j).categories.size();
				freq = new Vector<Integer>();
				for (int cat=0; cat<numOfCatOfAttrJ; cat++) // for every category			
					freq.add(cat,0);
			
				for (int p=0; p<result.data.size(); p++) // for every row (pattern)			
				{
					if (result.membership[i][p]==1)
					{
						int cat = ((Integer)result.data.get(p).properties.get(j)).intValue(); 
						freq.set(cat,freq.get(cat).intValue()+1);
					}
				}
				result.frequences.get(i).put(j,freq);
				result.clustersCenters.get(i).add(j,getMax(freq));
			}
		}
		
	}

	@Override
	public int getMinDistance(int pattern) 
	{
		int index=getClusterNum(pattern);
		Integer minDist=null;
		if (index!=-1)
			minDist = getDistance(result.data.get(pattern).properties,result.clustersCenters.get(index));
		
		for (int i=0; i<result.clustersCenters.size(); i++)
		{
			int distance = getDistance(result.data.get(pattern).properties,result.clustersCenters.get(i));
			if (minDist==null || distance<minDist.intValue())
			{
				minDist = distance;
				index = i;
			}
		}
		
		return index;
	}

	public int getDistance(Vector<Integer> v1, Vector<Integer> v2)	
	{
		int count=0;
		for (int i=0; i<v1.size(); i++)
			if (v1.get(i).intValue()!= v2.get(i).intValue())
				count++;
		
		return count;
	}
	
	@Override
	public void assignInitialClusterCenters(int n[])
	{
		if (n==null) n = generateRandomClusterCenters(result);	
		
		result.clustersCenters = new HashMap<Integer,Vector<Integer>>();
		
		for (int i=0; i<result.numOfClusters; i++)
			for (int j=0; j<result.data.size(); j++)
				result.membership[i][j] = 0;
		
		for (int i=0; i<result.numOfClusters; i++) // for every col (cluster)
		{
			result.clustersCenters.put(i,result.data.get(n[i]).properties);
			//result.membership[i][n[i]] = 1;
			//result.patternCountPerCluster[i]++;
		}	
		
	}
	
	public int getMax(Vector<Integer> v)
	{
		Integer max = new Integer(-1);
		int index = -1;
		for (int i=0; i<v.size(); i++)
			if (v.get(i).intValue()>max.intValue())
			{
				max = v.get(i).intValue();
				index = i;
			}
		
		return index;
	}

	@Override
	public void calculateClustersCenters(int pattern, int from, int to) 
	{
		int numOfAtrr = result.data.get(0).properties.size();		
		Vector<Integer> prop = result.data.get(pattern).properties;
				
		
		for (int j=0; j<numOfAtrr; j++)
		{				
			int i = result.frequences.get(from).get(j).get(prop.get(j));
			
			result.frequences.get(from).get(j).set(prop.get(j),i-1);			
			result.clustersCenters.get(from).set(j,getMax(result.frequences.get(from).get(j)));
			
			i = result.frequences.get(to).get(j).get(prop.get(j));
			
			result.frequences.get(to).get(j).set(prop.get(j),i+1);			
			result.clustersCenters.get(to).set(j,getMax(result.frequences.get(from).get(j)));
			
		}
		
	}
}
