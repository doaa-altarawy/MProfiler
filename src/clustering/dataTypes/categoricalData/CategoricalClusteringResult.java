package clustering.dataTypes.categoricalData;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;
import java.util.Map.Entry;

import clustering.util.Util;

import clustering.dataTypes.ClusteringResult;
import clustering.dataTypes.Pattern;

public class CategoricalClusteringResult<T, T2> extends ClusteringResult
{	
	public Vector<CategoricalPropertyDefinition> properties;	// labels for categorical properties
	
	/**
	 * @param T is either a vector<Integer> : for regular cluser center (mode)
	 * 						Vector<Vector<Double>> : for <value,RelativeFreq> pair 
	 */
	public Map <Integer,T> clustersCenters;	// Map of (centers) <clusterIndex, properties>		
	public double membership[][]; 			// 2D array of [clusterIndex][pattenIndex]
	public Map<String,Integer> clusterNames;	// <clustername,numOfPatternsInThatCluster>
	public Vector<HashMap<Integer, Vector<T2>>> frequences;
	public int patternCountPerCluster[];
	
	public CategoricalClusteringResult(Vector<CategoricalPropertyDefinition> prop, Vector<Pattern> data, int numOfClusters) 
	{
		super(data, data.size(), numOfClusters);
		properties = prop;
		membership = new double[numOfClusters][data.size()];	
		patternCountPerCluster = new int[numOfClusters];
		for (int i=0; i<numOfClusters; i++)
			patternCountPerCluster[i] = 0;
		buildClusterNamesFromData();
	}
	
	public int getNumOfPatternsInCluster(int cluster)
	{
		int n = 0;
		
		for (int i=0; i<membership[cluster].length; i++)
		{
			if (!Util.equal(membership[cluster][i],Util.ZERO))
				n++;
		}		
		
		return n;
	}
	
	private void buildClusterNamesFromData()
	{
		clusterNames = new HashMap<String, Integer>();
		
		for (int i=0; i<data.size(); i++)			
			{
				String correctCluster = data.get(i).correctCluster;
				if (clusterNames.containsKey(correctCluster))
					clusterNames.put(correctCluster, clusterNames.get(correctCluster)+1);
				else
					clusterNames.put(correctCluster, 1);
			}
		
		
	}
	
	public HashMap<String, Integer> getPatternDistributionInCluster(int clusterNum)
	{
		HashMap<String, Integer> temp = new HashMap<String, Integer>();
		
		for (Iterator i=clusterNames.entrySet().iterator(); i.hasNext();)
		{
			Entry<String, Integer> entry = (Entry<String, Integer>)i.next();
			temp.put(entry.getKey(),0);
		}
		
		for (int i=0; i<data.size(); i++)
			if (membership[clusterNum][i]==1)
			{
				String correctCluster = data.get(i).correctCluster;
				temp.put(correctCluster, temp.get(correctCluster)+1);
				
			}
		
		
		int max = -1;
		String cluster=null;		
		for (Iterator i=temp.entrySet().iterator(); i.hasNext();)
		{
			Entry<String, Integer> entry = (Entry<String, Integer>)i.next();
			if (entry.getValue()>max)
			{
				max = entry.getValue();
				cluster = entry.getKey();
			}
		}
		
		//return cluster;
		
		return temp;
	}
}
