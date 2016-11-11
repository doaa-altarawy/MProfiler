package clustering.distances;

import java.util.Vector;

import clustering.dataTypes.HardClusteringResult;
import clustering.util.Util;

public class BinaryEuclidian extends Distance {

	
	@Override
	public Object getDistanceSquare(Vector v1, Vector v2) 
	{
		Double distance = 0.0;
		Double l = new Double(v1.size());
		Double b = 0.0;
		
		if (v1.size()!= v2.size()) return null;
		
		for (int i=0; i<v1.size(); i++)
			if (Util.equal(((Double)v1.get(i)).doubleValue(),((Double)v2.get(i)).doubleValue()))
				b++;
		
		distance = (l-b) / l;
		
		return distance;
	}

	@Override
	public boolean calculateClustersCenters(HardClusteringResult result, int i)
	{
		// TODO Auto-generated method stub
		return false;
	}

}
