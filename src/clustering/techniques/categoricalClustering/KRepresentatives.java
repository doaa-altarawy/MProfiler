package clustering.techniques.categoricalClustering;

import java.util.HashMap;
import java.util.Vector;

import clustering.util.Util;


public class KRepresentatives extends CategoricalClusteringKMeans< Vector<Vector<Double>>, Double> 
{

	@Override
	public void calculateClustersCenters() 
	{	
		int numOfAtrr = result.properties.size(); 
		Vector<Double> freq;
				
		result.clustersCenters = new HashMap<Integer, Vector<Vector<Double>>>();
		result.frequences = new Vector<HashMap<Integer, Vector<Double>>>();
		
				
		for (int i=0; i<result.numOfClusters; i++) // for every col (cluster)
		{
			result.clustersCenters.put(i,new Vector<Vector<Double>>());
			result.frequences.add(i,new HashMap<Integer, Vector<Double>>());
			
			for (int j=0; j<numOfAtrr; j++)
			{				
				int numOfCatOfAttrJ = result.properties.get(j).categories.size();
								
				freq = new Vector<Double>();							
				for (int cat=0; cat<numOfCatOfAttrJ; cat++) // for every category			
					{ freq.add(cat,0.0); }
			
				int numOfPatInCluster=0;
				for (int p=0; p<result.data.size(); p++) // for every row (cluster)			
				{
					if (result.membership[i][p]==1)
					{
						numOfPatInCluster++;
						int cat = ((Integer)result.data.get(p).properties.get(j)).intValue(); 
						freq.set(cat,freq.get(cat)+1);												
					}
				}
				result.frequences.get(i).put(j,(Vector<Double>)(freq.clone()));
				for (int cat=0; cat<freq.size(); cat++)
					freq.set(cat, freq.get(cat)/(double)numOfPatInCluster);
				
				result.clustersCenters.get(i).add(j,freq);
			}
		}
		
	}

	@Override
	public void calculateClustersCenters(int pattern, int from, int to) 
	{					
		Vector<Double> temp, freq;
		int numOfAtrr = result.data.get(0).properties.size();		
		Vector<Integer> prop = result.data.get(pattern).properties;
		double i;		
		
		for (int j=0; j<numOfAtrr; j++)
		{				
			temp = new Vector<Double>();
			
			i = result.frequences.get(from).get(j).get(prop.get(j));
			
			result.frequences.get(from).get(j).set(prop.get(j),i-1);			
			freq = result.frequences.get(from).get(j);
					
			
			for (int cat=0; cat<freq.size(); cat++)
				temp.add(cat, (double)freq.get(cat)/(double)result.patternCountPerCluster[from]);
			
			result.clustersCenters.get(from).set(j,temp);
			
			
			///////
			
			
			temp = new Vector<Double>();
			
			i = result.frequences.get(to).get(j).get(prop.get(j));
			
			result.frequences.get(to).get(j).set(prop.get(j),i+1);			
			freq = result.frequences.get(to).get(j);
				
			
			for (int cat=0; cat<freq.size(); cat++)
				temp.add(cat, (double)freq.get(cat)/(double)result.patternCountPerCluster[to]);
			
			result.clustersCenters.get(to).set(j,temp);
			
			
			
		}
		
		
	}
	@Override
	public int getMinDistance(int pattern) 
	{
		int index=getClusterNum(pattern);
		double distance;
		Double minDist=null;
		if (index!=-1)
			minDist = getDistance(result.data.get(pattern).properties,result.clustersCenters.get(index));
		
		for (int i=0; i<result.clustersCenters.size(); i++)
		{
			distance = getDistance(result.data.get(pattern).properties,result.clustersCenters.get(i));
			if (minDist==null || distance<minDist.intValue())
			{
				minDist = distance;
				index = i;
			}
		}
		
		return index;
	}

	public double getDistance(Vector<Integer> v1, Vector<Vector<Double>> v2)	
	{
		double count=0.0;
		for (int i=0; i<v1.size(); i++)
			count+= (1 - v2.get(i).get(v1.get(i)));
		
		return count;
	}
	
	@Override
	public void assignInitialClusterCenters(int n[])
	{
		int numOfAttr = result.properties.size();
		
		if (n==null) n = generateRandomClusterCenters(result);	
		
		result.clustersCenters = new HashMap<Integer,Vector<Vector<Double>>>();
		
		for (int i=0; i<result.numOfClusters; i++)
			for (int j=0; j<result.data.size(); j++)
				result.membership[i][j] = 0;
		
		for (int i=0; i<result.numOfClusters; i++) // for every col (cluster)
		{
			Vector<Vector<Double>> temp = new Vector<Vector<Double>>();
			
			for (int attr=0; attr<numOfAttr; attr++)
			{
				Vector<Double> attribute = new Vector<Double>();
				int numOfCatOfAttr = result.properties.get(attr).categories.size();
				for (int cat=0; cat<numOfCatOfAttr; cat++)
				{
					attribute.add(new Double(0));
				}
				attribute.set(((Integer)result.data.get(n[i]).properties.get(attr)).intValue(), new Double(1.0));
				temp.add(attr, attribute);
				
			}
			result.clustersCenters.put(i,temp);
			//result.membership[i][n[i]] = 1;
			//result.patternCountPerCluster[i]++;
		}
		
		
	}

	
}
