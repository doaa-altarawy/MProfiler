package clustering.distances;


import java.util.Vector;

import clustering.dataTypes.HardClusteringResult;

public class Euclidian extends Distance 
{
	@Override
	public Object getDistanceSquare(Vector v1, Vector v2) 
	{
		Double distance = new Double(0);
		
		
		if (v1.size()!= v2.size()) return null;
		
		// find v1'.v1
		double d1=0;
		for (int i=0; i<v1.size(); i++)
			d1 = d1 + (Double)v1.get(i) * (Double)v1.get(i);
		
		// find v1'.v2
		double d2=0;
		for (int i=0; i<v1.size(); i++)
			d2 = d2 + (Double)v1.get(i) * (Double)v2.get(i);
		
		
		// find v2'.v2
		double d3=0;
		for (int i=0; i<v2.size(); i++)
			d3 = d3 + (Double)v2.get(i) * (Double)v2.get(i);
		
		
		// find Euclidian Distance Square
		distance = d1 - 2.0 *d2 + d3;
		
		return distance;
	}

	@Override
	public boolean calculateClustersCenters(HardClusteringResult result, int i)
	{
		// TODO Auto-generated method stub
		return false;
	}

	
}
