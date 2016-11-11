package clustering.techniques.hierarchialClustering;



public class Centroid extends HierarchialClustering
{
	@Override
	protected double getai(int i, int j, int k) 
	{				
		double n1 = result.getNumOfPatternsInCluster(i);
		double n2 = result.getNumOfPatternsInCluster(j);	
				
		return n1/(n1+n2);
	}

	@Override
	protected double getaj(int i, int j, int k) 
	{		
		double n1 = result.getNumOfPatternsInCluster(i);
		double n2 = result.getNumOfPatternsInCluster(j);		
				
		return n2/(n1+n2);
				
	}

	@Override
	protected double getb(int i, int j, int k) 
	{		
		double n1 = result.getNumOfPatternsInCluster(i);
		double n2 = result.getNumOfPatternsInCluster(j);		
				
		return - n1*n2/Math.pow(n1+n2,2);
	}

	@Override
	protected double getc(int i, int j, int k) 
	{		
		return 0;
	}

}
