package clustering.distances;


import java.util.Vector;
import org.biojava.bio.motif.FuzzyMotif;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.tools.MotifTools;

import clustering.dataTypes.ClusteringResult;
import clustering.dataTypes.FuzzySoftClusteringResult;
import clustering.dataTypes.HardClusteringResult;
import clustering.util.Util;

public class MotifDistance extends Distance 
{
	public static double MEMBER = 0.7;
	
	@Override
	public Object getDistanceSquare(Vector v1, Vector v2)
	{
		Double distance = (Double) getDistance(v1,v2);
		distance = Math.pow(distance, 2);
		
		return distance;
	}
	
	//////////////
	@Override
	public Object getDistance(Vector v1, Vector v2)
	{
		return 1.0 - ((Motif)v1.firstElement()).sim((Motif)v2.firstElement());
	}
	
	public boolean calculateClustersCenters(ClusteringResult result)
	{		
		if (result instanceof HardClusteringResult)
			return calculateClustersCenters((HardClusteringResult)result);
		else
			return calculateClustersCenters((FuzzySoftClusteringResult)result);
	}
	
	public boolean calculateClustersCenters(HardClusteringResult result)
	{
		boolean flag = false;		
		Motif intersection = null;
				
		for (int i=0; i<result.membership.length; i++) // for every row (cluster)
		{					
			// get Intersections 
			Vector<Motif> motifs = result.getCluster(i);
			try{
			intersection = MotifTools.getPerwiseIntersection(motifs);
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			// if no change in all cluster centers return false
			if (!Util.equalZero(intersection.sim((Motif)result.clustersCenters.get(i).firstElement())))
				flag = true;				
			
			Vector<Motif> temp = new Vector<Motif>(1);
			temp.add(intersection);
			result.clustersCenters.put(i, temp);
		}		
		
		if (!flag)
			System.out.println("No change in Clutsers Centers..");
		
		return flag;
	}
	
	public boolean calculateClustersCenters(HardClusteringResult result, int i)
	{
		boolean flag = false;		
		Motif intersection = null;
				
		// get Intersections 
		Vector<Motif> motifs = result.getCluster(i);
		try{
		intersection = MotifTools.getPerwiseIntersection(motifs);
		}
		catch (Exception e) {e.printStackTrace();}
		
		// if no change in all cluster centers return false
		if (!Util.equalZero(intersection.sim((Motif)result.clustersCenters.get(i).firstElement())))
			flag = true;				
		
		Vector<Motif> temp = new Vector<Motif>(1);
		temp.add(intersection);
		result.clustersCenters.put(i, temp);
		
		if (!flag)
			System.out.println("No change in Clutsers Centers..");
		
		return flag;
	}
	/*
	public boolean calculateClustersCenters(HardClusteringResult result)
	{
		boolean flag = false;		
		Motif intersection;
				
		for (int i=0; i<result.membership.length; i++) // for every row (cluster)
		{			
			intersection = new Motif(((Motif)result.data.get(0).properties.firstElement()).getDataset());
			
			// get Intersections 
			for (int j=0; j<result.membership[i].length; j++) // for all member patterns
			{
				if (!Util.equalZero(result.membership[i][j]))
				{	
					Motif m = (Motif)result.data.get(j).properties.firstElement();
					for (int k=j+1; k<result.membership[i].length; k++)
					{
						if (!Util.equalZero(result.membership[i][k]))
						{
							Motif m2 = (Motif)result.data.get(k).properties.firstElement();
							if (m.getFinder().equals(m2.getFinder())) continue;
							try{intersection.merge(m.intersect(m2));}
							catch (Exception e){e.printStackTrace();}
						}
					}
				}
			}
			intersection.removeRedundantFast();
						
			// if no change in all cluster centers return false
			if (!Util.equalZero(intersection.sim((Motif)result.clustersCenters.get(i).firstElement())))
				flag = true;				
			
			Vector<Motif> temp = new Vector<Motif>(1);
			temp.add(intersection);
			result.clustersCenters.put(i, temp);
		}		
		
		if (!flag)
			System.out.println("No change in Clutsers Centers..");
		
		return flag;
	}
	*/
	// revise
	public boolean calculateClustersCenters(FuzzySoftClusteringResult result)
	{
		boolean flag = false;		
		FuzzyMotif intersection;
				
		for (int i=0; i<result.membership.length; i++) // for every row (cluster)
		{			
			intersection = new FuzzyMotif(((Motif)result.data.get(0).properties.firstElement()).getDataset());
			
			// get Intersections 
			for (int j=0; j<result.membership[i].length; j++) // for all member patterns
			{
				if (!Util.equalZero(result.membership[i][j]))
				{	
					FuzzyMotif m = new FuzzyMotif((Motif)result.data.get(j).properties.firstElement(), result.membership[i][j]);
					for (int k=j+1; k<result.membership[i].length; k++)
					{
						if (!Util.equalZero(result.membership[i][k]))
						{
							Motif m2 = (Motif)result.data.get(k).properties.firstElement();
							if (m.getFinder().equals(m2.getFinder())) continue;
							try{intersection.merge(m.intersect(m2, result.membership[i][k]));}
							catch (Exception e){e.printStackTrace();}
						}
					}
				}
			}
			intersection.removeRedundant();
						
			// if no change in all cluster centers return false
			if (!Util.equalZero(intersection.sim((FuzzyMotif)result.clustersCenters.get(i).firstElement())))
				flag = true;				
			
			Vector<FuzzyMotif> temp = new Vector<FuzzyMotif>(1);
			temp.add(intersection);
			result.clustersCenters.put(i, temp);
		}		
		
		if (!flag)
			System.out.println("No change in Clutsers Centers..");
		
		return flag;
	}
	/*
	public boolean calculateClustersCenters(FuzzySoftClusteringResult result)
	{
		boolean flag = false;		
		Vector<Double> sum;
		Double sumWi;
		
		for (int i=0; i<result.numOfClusters; i++) // for every cluster
		{
						
			sum = new Vector<Double>(result.data.get(i).properties.size());
			sumWi = new Double(0);
			
			for (int j=0; j<result.data.get(i).properties.size(); j++)
				sum.add(j, new Double(0));				
					
						
			// sum member vectors
			for (int j=0; j<result.numOfPatterns; j++) // for all patterns
			{	
				double d = Math.pow(result.membership[i][j],result.m);
				sumWi += d;
				
				Vector<Double> x =  new Vector<Double>(result.data.get(j).properties);
				for (int k=0; k<x.size(); k++) // for all properties
					sum.set(k, sum.get(k)+ x.get(k)* d);				
			}
									
			for (int k=0; k<sum.size(); k++)
			{
				if (!Util.equal(sumWi, Util.ZERO))
					sum.set(k, sum.get(k)/sumWi);
				else
				{
					System.out.println("sumWiiiiiiiiiiiiiiiiiii="+sumWi);					
					System.out.println("cluster="+i);
					sum.set(k, sum.get(k)/sumWi);
				}
				
				//if (!Util.equal(sum.get(k), result.clustersCenters.get(i).get(k)))
				//	flag = true;					
			}			
			result.clustersCenters.put(i, new Vector<Double>(sum));
		}		
		
		return flag;
	}

	*/
	public boolean calculateClustersCenters(HardClusteringResult result, int i, int j, double sign)
	{
		boolean flag = false;
		double n;
		Vector<Motif> m, x;
		
		n = result.getNumOfPatternsInCluster(i);
		m = new Vector<Motif>(result.clustersCenters.get(i));
		x = result.data.get(j).properties;
	
		
		for (int k=0; k<m.size(); k++ )
		{
			if (n!=0)
				m.get(k).merge(x.get(k));
			else
				m.set(k,new Motif());
			
			
			if (!Util.equalZero(m.get(k).sim((Motif)result.clustersCenters.get(i).get(k))))
				flag = true;			
			
		}
		
		result.clustersCenters.put(i, new Vector<Motif>(m));
		
		return flag;
	}
	
	
}
