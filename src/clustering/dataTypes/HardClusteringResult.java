package clustering.dataTypes;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;

import org.biojava.bio.motif.Motif;

import clustering.util.Util;

public class HardClusteringResult extends ClusteringResult
{
		
	//public Map <Integer,Vector> clustersCenters;// Map of (centers) <clusterIndex, vector of properties>	
	public Vector<Double> objectiveFunctionTrace;		// value of objectiveFunc in every iteration	
	
		
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	/**	 
	 * when num of clusters in known or estimated
	 */
	public HardClusteringResult(Vector<Pattern> data, int numOfPatterns,int numOfClusters )
	{
		super(data, numOfPatterns, numOfClusters);
		
		clustersCenters = new HashMap<Integer,Vector>();
		objectiveFunctionTrace = new Vector<Double>();		
		membership = new Double[numOfClusters][numOfPatterns];		
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
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
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public int getClusterOfPattern(int pattern)
	{
		for (int i=0; i<numOfClusters; i++)
		{
			if (membership[i][pattern]==1)
				return i;
		}
		
		return -1;		// shouldn't get here
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public Vector<Motif> getCluster(int i)
	{
		Vector<Motif> cluster = new Vector<Motif>();
		for (int j=0; j<numOfPatterns; j++)
		{
			if (!Util.equalZero(membership[i][j]))
				cluster.add((Motif)data.get(j).properties.firstElement());		
		}
		
		return cluster;
	}
}
