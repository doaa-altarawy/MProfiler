package clustering.techniques.hardClustering;

import clustering.util.*;


public class KMeans extends HardClustering
{
		
	@Override
	/*
	public void partition() 
	{
		boolean flag = true;
		
				
		do
		{
			System.out.println("Iteration num:"+result.numOfIterations);
			// step 2: distribute patterns among clusters			
			for (int i=0; i<m; i++)
			{				
				int min = ((Integer)distance.getMinDistance(result.data.get(i).properties, result.clustersCenters)).intValue();
				int max = ((Integer)distance.getMaxDistance(result.data.get(i).properties, result.clustersCenters)).intValue();
				
				if (min==0)
				{	
					for (int j=0; j<k; j++)
						if (j!=min)
							result.membership[j][i] = 0.0;
					result.membership[min][i] = 1.0;
				}
				else
				{
					for (int j=0; j<k; j++)
						if (j!=max)
							result.membership[j][i] = 0.0;
					result.membership[max][i] = 1.0;
				}
			}
			
			// step 3: get clusters centers
			flag = distance.calculateClustersCenters(result); 			
			
			System.out.println("Calc Objective Function....");
			result.objectiveFunctionTrace.add(calculateObjectiveFunctionWeightRatio());
			
			int maxIndex = result.objectiveFunctionTrace.size()-1;
			if (maxIndex>1 && Util.equal(result.objectiveFunctionTrace.lastElement(),result.objectiveFunctionTrace.get(maxIndex-1)))
				flag = false;
				
			result.numOfIterations++;
			if (result.numOfIterations>Util.maxNumOfIterations)
			{
				System.out.println("Max # of itrations exceeded");
				flag=false;
			}
		}
		while (flag);
		System.out.println("Number of Iterations: "+result.numOfIterations);
	}

	*/
	 public void partition() 
	{
		boolean flag = true;
		
				
		do
		{
			System.out.println("Iteration num:"+result.numOfIterations);
			// step 2: distribute patterns among clusters			
			for (int i=0; i<m; i++)
			{				
				int min = ((Integer)distance.getMinDistance(result.data.get(i).properties, result.clustersCenters)).intValue();
									
				for (int j=0; j<k; j++)
					if (j!=min)
						result.membership[j][i] = 0.0;
				result.membership[min][i] = 1.0;
			}
			
			// step 3: get clusters centers
			flag = distance.calculateClustersCenters(result); 			
			
			//System.out.println("Calc Objective Function....");
			result.objectiveFunctionTrace.add(calculateObjectiveFunction());
			System.out.println("Objective Function="+result.objectiveFunctionTrace.lastElement());
			
			int maxIndex = result.objectiveFunctionTrace.size()-1;
			if (maxIndex>1 && Util.equal(result.objectiveFunctionTrace.lastElement(),result.objectiveFunctionTrace.get(maxIndex-1)))
				flag = false;
				
			result.numOfIterations++;
			if (result.numOfIterations>Util.maxNumOfIterations)
			{
				System.out.println("Max # of itrations exceeded");
				flag=false;
			}
		}
		while (flag);
		System.out.println("Number of Iterations: "+result.numOfIterations);
	}

		
}
