package clustering.techniques.fuzzySoftClustering;

import clustering.util.Util;
import java.util.*;

public class SoftTN extends FuzzySoftClustering
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
				for (int j=0; j<result.numOfClusters; j++)
					result.pseudoMembership[j][i] = Math.pow(result.distanceMatrix[j][i], 1.0/(1.0-result.m));			
		}
		
		// update membership (soft clustering part)
		for (int i=0; i<result.numOfPatterns; i++)
		{			
			Double alfa = getLargestTNElementOf(i);	// smallest value
			
			for (int j=0; j<result.numOfClusters; j++)
				if (result.pseudoMembership[j][i]<alfa)
					result.pseudoMembership[j][i] = 0.0;
			
		}// for every pattern		
		
	}

	Double getLargestTNElementOf(int pattern)
	{
		Vector<Double> list = new Vector<Double>();
		
		for (int i=0; i<result.numOfClusters; i++)
			list.add(result.pseudoMembership[i][pattern]);
		
		Collections.sort(list);
		
		return list.get(list.size()-result.TN);
	}

}
