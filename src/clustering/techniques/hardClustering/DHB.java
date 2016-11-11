package clustering.techniques.hardClustering;

import clustering.util.Util;

public class DHB extends HardClustering 
{

	@Override
	public void partition() 
	{
		boolean flag = false;
		int cluster, newCluster, n, nj;
		Double delta, deltaMin, deltaJ;
		
		//	step 1: choose k initial clusters
		//partitionAtRandom();
		
		// step 2: calculate cluseters centers and objective function
		distance.calculateClustersCenters(result);
		result.objectiveFunctionTrace.add(calculateObjectiveFunction());
		
		// step 3: core iteration block for DHB
		do
		{		
			
			flag = false;
			for (int i=0; i<m; i++)
			{
				cluster = result.getClusterOfPattern(i);
				newCluster = cluster;
				n = result.getNumOfPatternsInCluster(cluster);
				if (n!=1)
				{
					delta = new Double(((double)n /(double)(n-1)) * ((Double)distance.getDistanceSquare(result.data.get(i).properties, result.clustersCenters.get(cluster))).doubleValue());
					deltaMin = new Double(delta);	///////
					
					for (int j=0; j<k; j++)
					{
						if (i!=j)
						{
							nj = result.getNumOfPatternsInCluster(j);
							deltaJ = new Double(((double)nj /(double)(nj+1)) * ((Double)distance.getDistanceSquare(result.data.get(i).properties, result.clustersCenters.get(j))).doubleValue());
							if (deltaJ < deltaMin )
							{
								deltaMin = deltaJ;
								newCluster = j;
							}
						}// end (i!=j)						
					}// end for j	
					
					if (deltaMin < delta && cluster!= newCluster)
					{
						// move pattern i to newCluster
						result.membership[cluster][i] = 0.0;
						result.membership[newCluster][i] = 1.0;
						distance.calculateClustersCenters(result, cluster, i, -1);
						distance.calculateClustersCenters(result, newCluster, i, +1);
						Double d = calculateObjectiveFunction();
						if (!Util.equal(result.objectiveFunctionTrace.lastElement(),d))
							flag = true;
						result.objectiveFunctionTrace.add(d);						
					}
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
