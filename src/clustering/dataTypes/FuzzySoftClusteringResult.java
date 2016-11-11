package clustering.dataTypes;

import java.util.Map;
import java.util.Vector;

public class FuzzySoftClusteringResult extends ClusteringResult 
{	
	public int m;		// power of membership matrix
	public int TN;		// max num of clusters per pattern
	public double TD;	// min accepted distance
	public double TW;	// min accepted weight
	
		
	public Double pseudoMembership[][]; 		// 2D array of [clusterIndex][pattenIndex], values are not restricted 
	public Double distanceMatrix[][];
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	
	public FuzzySoftClusteringResult(Vector<Pattern> data, int m, int k, Object T) 
	{
		super(data, data.size(), k);
		
		this.m = m;
		
		try {TN = (Integer)T; if (TN>k) TN=k;}
		catch (Exception e){TN = k;}
		
		try {TD = (Double)T;}
		catch (Exception e){TD=0;}
		
		try {TW = (Double)T / (double)k ; if (TW>1) TW=1.0;}
		catch (Exception e){TW=1.0;}
		
		membership = new Double[numOfClusters][numOfPatterns];
		pseudoMembership = new Double[numOfClusters][numOfPatterns];
		distanceMatrix = new Double[numOfClusters][numOfPatterns];
		
	
	}
}
