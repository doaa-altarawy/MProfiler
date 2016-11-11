package clustering.dataTypes;

import java.util.Vector;

public class Pattern 
{
	public String label;		// label for the pattern	
	public Vector properties;	// properties of the pattern
	public Vector<Character> labels;
	public String correctCluster;
	
	public Pattern(){label = new String();}
	public Pattern(String label){this.label = label;}
	public Pattern(String label, Vector prop)
	{
		this.label = label;
		properties = new Vector(prop);
	}
}
