package clustering.techniques.categoricalClustering;

import java.util.Calendar;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;

import clustering.util.Util;

import clustering.dataTypes.FuzzySoftClusteringResult;
import clustering.dataTypes.Pattern;
import clustering.dataTypes.categoricalData.CategoricalClusteringResult;
import clustering.dataTypes.categoricalData.CategoricalPropertyDefinition;

public abstract class CategoricalClusteringKMeans<T, T2> 
{	
	public CategoricalClusteringResult<T, T2> result;
	
	public abstract void assignInitialClusterCenters(int n[]);
	public abstract void calculateClustersCenters();
	public abstract void calculateClustersCenters(int pattern, int from, int to);
	public abstract int getMinDistance(int pattern);
	
	
	public CategoricalClusteringResult<T, T2> partition(Vector<CategoricalPropertyDefinition> prop, Vector<Pattern> data, int numOfClusters, int clusterSizes[])
	{
				
		result = new CategoricalClusteringResult<T, T2>(prop, data, numOfClusters);				
				
		//partitionAtRandom(clusterSizes);
		assignInitialClusterCenters(clusterSizes);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		partition();
		long endTime = Calendar.getInstance().getTimeInMillis();
		
		result.runningTime = endTime - startTime ;
		
		return result;
	}

	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	private void partition()
	{
		boolean flag = false;
		int m = result.data.size();
		int k = result.numOfClusters;
		int min;
		int c=0;
		
				
		do
		{
			System.out.println("Iteration# "+result.numOfIterations);
			flag = false;
			
						
			// step 3: distribute patterns among clusters			
			for (int i=0; i<m; i++)
			{					
				min = getMinDistance(i);
				c=-1;			
			
				for (int j=0; j<k; j++)
					if (result.membership[j][i]==1) 
					{
						if (j!=min)
						{
							flag = true; // pattern will change its cluster j to cluster min
							c = j;
						}
						
						result.patternCountPerCluster[j]--;					
						result.membership[j][i] = 0;
					}
				result.membership[min][i] = 1;
				result.patternCountPerCluster[min]++;
				
//				step 2: get clusters centers
				if (result.numOfIterations>0 && c!=-1)
					calculateClustersCenters(i,c,min);					
				
			}
			
			if (result.numOfIterations==0)
				calculateClustersCenters();
			
			result.numOfIterations++;
			
			if (result.numOfIterations>Util.maxNumOfIterations)
			{
				System.out.println("Max # of itrations exceeded");
				flag=false;
			}
		}
		while (flag || result.numOfIterations==1);
	
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public CategoricalClusteringResult<T, T2> partition2(Vector<CategoricalPropertyDefinition> prop, Vector<Pattern> data, int numOfClusters, int clusterSizes[])
	{
				
		result = new CategoricalClusteringResult<T, T2>(prop, data, numOfClusters);				
				
		//partitionAtRandom(clusterSizes);
		assignInitialClusterCenters(clusterSizes);
		
		long startTime = Calendar.getInstance().getTimeInMillis();
		partition2();
		long endTime = Calendar.getInstance().getTimeInMillis();
		
		result.runningTime = endTime - startTime ;
		
		return result;
	}

	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	private void partition2()
	{
		boolean flag = false;
		int m = result.data.size();
		int k = result.numOfClusters;
		int min;
		int c=0;
		
				
		do
		{
			System.out.println("Iteration# "+result.numOfIterations);
			flag = false;
			
						
			// step 3: distribute patterns among clusters			
			for (int i=0; i<m; i++)
			{					
				min = getMinDistance(i);
				c=-1;			
			
				for (int j=0; j<k; j++)
					if (result.membership[j][i]==1) 
					{
						if (j!=min)
						{
							flag = true; // pattern will change its cluster j to cluster min
							c = j;
						}
						
						result.patternCountPerCluster[j]--;					
						result.membership[j][i] = 0;
					}
				result.membership[min][i] = 1;
				result.patternCountPerCluster[min]++;
				
//				
				calculateClustersCenters();
								
				
			}
			
			//if (result.numOfIterations==0)
							
			
			result.numOfIterations++;
			
			if (result.numOfIterations>Util.maxNumOfIterations)
			{
				System.out.println("Max # of itrations exceeded");
				flag=false;
			}
		}
		while (flag || result.numOfIterations==1);
	
	}

	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static int [] generateRandomClusterCenters(CategoricalClusteringResult result)
	{			
			
		int n[] = new int[result.numOfClusters];
		Object names[] = result.clusterNames.keySet().toArray();
	
		for (int i=0; i<result.numOfClusters; i++)
		{			
			int x = Util.generateRandomNumbers(1,result.data.size())[0];
			
			n[i] = x;
			
			if (names.length<=i) continue;
						
			if (names[i].equals("recommend"))
					{n[i] = 0; continue;}
			
			while (!names[i].equals(result.data.elementAt(x).correctCluster))
				x = Util.generateRandomNumbers(1,result.data.size())[0];
			
			n[i] = x;
				
		}	
		
		return n;
	}	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	/*public static int[] getRandomSizes(int k, int m)
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
	*/
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public int getClusterNum(int pattern)
	{
		for (int i=0; i<result.membership.length; i++)
			if (result.membership[i][pattern]==1)
				return i;
		
		return -1;
	}
	
}
