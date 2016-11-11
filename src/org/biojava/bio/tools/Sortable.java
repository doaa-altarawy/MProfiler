package org.biojava.bio.tools;


public class Sortable<T> implements Comparable<Sortable<T>>
{
	public T obj;
	public Double value = 0.0;
	boolean descending = true;
	
	public Sortable(T m)
	{
		obj = m;			
	}
	
	public Sortable(T m, double val)
	{
		obj = m;
		this.value = val;	
	}
	
	public Sortable(T m, double val, boolean descending)
	{
		obj = m;
		this.value = val;
		this.descending = descending;
	}

	// reversed to sort descending
	public int compareTo(Sortable<T> o)
	{
		if (value>o.value) return descending?-1:1;		// returns -ve if greater to sort descending
		else if (value<o.value) return descending?1:-1;	// returns +ve if smaller to sort descending
		else return 0;
	}
	
	public void setDescending(boolean b)
	{
		descending = b;
	}
}	
/*-----------------------------------------------------------------------------*/

