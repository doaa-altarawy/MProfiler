package clustering.techniques.hardClustering;

import clustering.util.Util;

public class DHF extends HardClustering 
{

	@Override
	public void partition() 
	{
		boolean flag = false;
		int cluster, n, nj;
		Double delta, deltaJ;
		
		//	step 1: choose k initial clusters
		//partitionAtRandom();
		
		// step 2: calculate clusters centers and objective function
		distance.calculateClustersCenters(result);
		result.objectiveFunctionTrace.add(calculateObjectiveFunction());
		
		// step 3: core iteration block for DHF
		do
		{
					
			flag = false;
			for (int i=0; i<m; i++)
			{
				cluster = result.getClusterOfPattern(i);				
				n = result.getNumOfPatternsInCluster(cluster);
				if (n!=1)
				{
					delta = new Double(((double)n /(double)(n-1)) * ((Double)distance.getDistanceSquare(result.data.get(i).properties, result.clustersCenters.get(cluster))).doubleValue());
					
					for (int j=0; j<k; j++)
					{
						if (i!=j)
						{
							nj = result.getNumOfPatternsInCluster(j);
							deltaJ = new Double(((double)nj /(double)(nj+1)) * ((Double)distance.getDistanceSquare(result.data.get(i).properties, result.clustersCenters.get(j))).doubleValue());
							if (deltaJ < delta  && cluster!= j)
							{
								//	move pattern i to newCluster
								result.membership[cluster][i] = 0.0;
								result.membership[j][i] = 1.0;
								distance.calculateClustersCenters(result, cluster, i, -1);
								distance.calculateClustersCenters(result, j, i, +1);
								Double d = calculateObjectiveFunction();
								if (!Util.equal(result.objectiveFunctionTrace.lastElement(),d))
									flag = true;
								result.objectiveFunctionTrace.add(d);							
								j=k;								
							}
						}// end (i!=j)						
					}// end for j					
				}// enf (n!=1)
			}// end for i
			
						
			result.numOfIterations++;
			if (result.numOfIterations>Util.maxNumOfIterations)
			{
				System.out.println("Max # of itrations exceeded");
				flag=false;
			}
		}
		while (flag);			

	}

}
