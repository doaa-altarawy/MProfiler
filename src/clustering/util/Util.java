package clustering.util;

import java.util.Calendar;
import java.util.Random;

public class Util 
{
	public static Double RANGE = 0.001;
	public static Double ZERO = RANGE / 1000.0;
	public static int maxNumOfIterations=25;
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	public static boolean equal(double d1, double d2)
	{
		if (Math.abs(d1-d2) < RANGE) return true;
		else return false;
		
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	public static boolean equalZero(double d)
	{
		if (Math.abs(d-ZERO) < RANGE) return true;
		else return false;
		
	}
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static int[] generateRandomNumbers(int num, int to)
	{
		Random random = new Random(Calendar.getInstance().getTimeInMillis());
		int init[] = new int[num];
		
		
		// generate random unique points		
		for (int i=0; i<num; i++)
			init[i] = random.nextInt(to);
						
		return init;
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static int[] generateRandomUniqueNumbers(int num, int to)
	{
		Random random = new Random(Calendar.getInstance().getTimeInMillis());
		int list[] = new int[num];
		
		
		// generate random unique points		
		for (int i=0; i<num; i++)
		{
			boolean unique = true;			
			do
			{
				unique = true;
				list[i] = random.nextInt(to);
				for (int j=0; j<num ; j++)
					if (list[i] == list[j] && i!=j)
						unique = false;
			}
			while (!unique);
		}
		
		return list;
	}
	
/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static int[] generateUnsortedList(int num)
	{
		Random random = new Random(Calendar.getInstance().getTimeInMillis());
		int list[] = new int[num];
		
		for (int i=0; i<num; i++)
			list[i] = i;
		
		// unsort list		
		for (int i=0; i<num; i++)
		{
			int n1 = random.nextInt(num);
			int n2 = random.nextInt(num);
			
			int temp = list[n1];
			list[n1] = list[n2];
			list[n2] = temp;			
			
		}
		
		return list;
	}
}
