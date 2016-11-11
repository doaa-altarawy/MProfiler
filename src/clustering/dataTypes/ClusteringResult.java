package clustering.dataTypes;

import java.util.*;


public abstract class ClusteringResult 
{
	public Vector<Pattern> data;		// data set
	public int numOfClusters;	
	public int numOfPatterns;
	
	public int numOfIterations;		
	public long runningTime;
	
	public Double membership[][]; 		// 2D array of [clusterIndex][pattenIndex]
	// fuzzy public Double membership[][]; 		// 2D array of [clusterIndex][pattenIndex], 0<= w <=1
	public Map <Integer,Vector> clustersCenters;// Map of (centers) <clusterIndex, vector of properties>
	
	public ClusteringResult( Vector<Pattern> data, int numOfPatterns)
	{		
		numOfIterations = 0;	
		this.numOfPatterns = numOfPatterns;
		this.data = new Vector<Pattern>(data.subList(0,numOfPatterns));
	}
	
	/**	 
	 * when num of clusters in known or estimated
	 */
	public ClusteringResult(Vector<Pattern> data,  int numOfPatterns, int numOfClusters)
	{		
		this.numOfClusters = numOfClusters;
		this.numOfPatterns = numOfPatterns;
		
		numOfIterations = 0;				
		this.data = new Vector<Pattern>(data.subList(0,numOfPatterns));
	}
}
