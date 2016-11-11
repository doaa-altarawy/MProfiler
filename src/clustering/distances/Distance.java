package clustering.distances;

import java.util.Map;
import java.util.Vector;

import clustering.dataTypes.ClusteringResult;
import clustering.dataTypes.HardClusteringResult;
import clustering.util.Util;

public abstract class Distance
{
	
	public abstract Object getDistanceSquare(Vector v1, Vector v2);
	
	public Object getDistance(Vector v1, Vector v2)
	{
		Double distance = (Double) getDistanceSquare(v1,v2);
		distance = Math.abs(distance);
		distance = Math.sqrt(distance);
		
		return distance;
	}
	 	
	public Object getMinDistance(Vector v, Map <Integer,Vector> z)
	{
		Integer min = new Integer(-1);
		Double minDist = null;
		
		for (int i=0; i<z.size(); i++)
		{
			Double d =  (Double) getDistance(v,z.get(i));
			if (minDist==null || minDist.doubleValue()> d.doubleValue())
			{
				min = i;
				minDist = d.doubleValue();
			}
		}
		
		return min;
	}
	
	public Object getMaxDistance(Vector v, Map <Integer,Vector> z)
	{
		Integer max = new Integer(-1);
		Double maxDist = null;
		
		for (int i=0; i<z.size(); i++)
		{
			Double d =  (Double) getDistance(v,z.get(i));
			if (maxDist==null || maxDist.doubleValue()< d.doubleValue())
			{
				max = i;
				maxDist = d.doubleValue();
			}
		}
		
		return max;
	}
	
	// assumes clusters centers are vectors of double (override it for other data types)
	public boolean calculateClustersCenters(ClusteringResult result)
	{
		boolean flag = false;
		double n;
		Vector<Double> sum;
		
				
		for (int i=0; i<result.membership.length; i++) // for every row (cluster)
		{
			//n = 0;
			n = ((HardClusteringResult)result).getNumOfPatternsInCluster(i);
			sum = new Vector<Double>(result.data.get(i).properties.size());
			
			for (int j=0; j<result.data.get(i).properties.size(); j++)
				sum.add(j, new Double(0));
			
			// sum member vectors
			for (int j=0; j<result.membership[i].length; j++) // for all member patterns
				if (!Util.equal(result.membership[i][j], Util.ZERO))
				{
					//n++;
					Vector<Double> temp =  new Vector<Double>(result.data.get(j).properties);
					for (int k=0; k<temp.size(); k++ )
					{
						Double d = new Double(sum.get(k)+temp.get(k));
						sum.set(k, d);
					}
				}
			
			
			for (int k=0; k<sum.size(); k++ )
			{
				if (n!=0)
					sum.set(k,sum.get(k)/n);
				else
					sum.set(k,new Double(0));
				if (!Util.equal((Double)sum.get(k),(Double)result.clustersCenters.get(i).get(k)))
					flag = true;				
			}
			
			result.clustersCenters.put(i, new Vector<Double>(sum));
		}		
		
		return flag;
	}
	public abstract boolean calculateClustersCenters(HardClusteringResult result, int i);
	
	/**
	 * Todo: wrong, always add the new pattern , sign = -1 is disabled
	 * 
	 */ 
	public boolean calculateClustersCenters(HardClusteringResult result, int i, int j, double sign)
	{
		boolean flag = false;
		double n;
		Vector<Double> m, x;
		
		n = result.getNumOfPatternsInCluster(i);
		m = new Vector<Double>(result.clustersCenters.get(i));
		x = result.data.get(j).properties;
	
		
		for (int k=0; k<m.size(); k++ )
		{
			if (n!=0)
				m.set(k,m.get(k)+sign*(x.get(k)-m.get(k))/n);
			else
				m.set(k,new Double(0));
			
			
			if (!Util.equal((Double)m.get(k),(Double)result.clustersCenters.get(i).get(k)))
				flag = true;				
			
		}
		
		result.clustersCenters.put(i, new Vector<Double>(m));
		
		return flag;
	}
}
