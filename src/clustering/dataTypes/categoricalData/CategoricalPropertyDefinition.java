package clustering.dataTypes.categoricalData;

import java.util.HashMap;

/**
 * 
 * @author Doaa Altarawy
 *
 */

public class CategoricalPropertyDefinition
{
	public String propertyName;
	public HashMap<String,Integer> categories;
	
	public CategoricalPropertyDefinition(){}
	
	public CategoricalPropertyDefinition(String name) 
	{
		propertyName = name;
	}
	
	public CategoricalPropertyDefinition(String name, HashMap<String,Integer> cat) 
	{
		propertyName = name;
		categories = cat;
	}
	
	
	
}
