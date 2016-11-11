package clustering.techniques.hierarchialClustering;

public class Median extends HierarchialClustering 
{
	@Override
	protected double getai(int i, int j, int k) 
	{		
		return 0.5;
	}

	@Override
	protected double getaj(int i, int j, int k) 
	{		
		return 0.5;
	}

	@Override
	protected double getb(int i, int j, int k) 
	{		
		return -025;
	}

	@Override
	protected double getc(int i, int j, int k) 
	{		
		return 0;
	}

}
