package clustering.techniques.fuzzySoftClustering;

import clustering.util.Util;

public class SoftTW extends FuzzySoftClustering
{
		
	@Override
	protected void calculatePseudoMembershipMatrix() 
	{
		int r;
		
		// compute pseudo membership
		for (int i=0; i<result.numOfPatterns; i++ ) // for each pattern
		{
			r = -1;			
			for (int j=0; j<result.numOfClusters; j++)	// search for a zero distance
				if (Util.equal(result.distanceMatrix[j][i], Util.ZERO))
					r = j;
			
			if (r != -1)
			{
				for (int j=0; j<result.numOfClusters; j++)
					result.pseudoMembership[j][i] = new Double(0);
				result.pseudoMembership[r][i] = new Double(1);				
			}
			else
			{
				double segma = 0;
				for (int k=0; k<result.numOfClusters; k++)
					segma += Math.pow(result.distanceMatrix[k][i], 1.0/(1.0-result.m));
				
				for (int j=0; j<result.numOfClusters; j++)
				{
					double lamda = Math.pow(result.distanceMatrix[j][i], 1.0/(1.0-result.m)) / segma;
					if(lamda < result.TW)
						result.pseudoMembership[j][i] = 0.0;
					else result.pseudoMembership[j][i] = lamda;			
				}
			} // else	
		}
		
		// update membership (soft clustering part)
		
		
	}

		

}
