package clustering.techniques.hierarchialClustering;

public class WardMinVariance extends HierarchialClustering 
{
	@Override
	protected double getai(int i, int j, int k) 
	{				
		double n1 = result.getNumOfPatternsInCluster(i);
		double n2 = result.getNumOfPatternsInCluster(j);
		double n3 = result.getNumOfPatternsInCluster(k);	
				
		return (n1+n3)/(n1+n2+n3);
	}

	@Override
	protected double getaj(int i, int j, int k) 
	{		
		double n1 = result.getNumOfPatternsInCluster(i);
		double n2 = result.getNumOfPatternsInCluster(j);
		double n3 = result.getNumOfPatternsInCluster(k);	
				
		return (n2+n3)/(n1+n2+n3);
				
	}

	@Override
	protected double getb(int i, int j, int k) 
	{		
		double n1 = result.getNumOfPatternsInCluster(i);
		double n2 = result.getNumOfPatternsInCluster(j);	
		double n3 = result.getNumOfPatternsInCluster(k);	
				
		return - n3/(n1+n2+n3);
	}

	@Override
	protected double getc(int i, int j, int k) 
	{		
		return 0;
	}
}
