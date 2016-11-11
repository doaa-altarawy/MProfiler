package clustering.util;

import java.io.*;
import java.util.StringTokenizer;
import java.util.*;

import clustering.dataTypes.Pattern;
import clustering.dataTypes.categoricalData.CategoricalPropertyDefinition;

public class DataReader 
{
	/**
	 * Reads data from a file..
	 * for each data, it reads:
	 * 	1- Pattern id
	 * 	2- a vector of properties
	 * @param fileName : source of the data
	 * @return a vector of patterns in the data file
	 * @author Doaa Altarwy
	 */
	public static Vector<Pattern> getNumericalData(String fileName) throws Exception
	{
		String space = "\t";		// default separator is "\t" (tab)
				
		return getNumericalData(fileName, space);
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static Vector<Pattern> getNumericalData(String fileName, String separator) throws Exception
	{
		Vector<Pattern> temp = new Vector<Pattern>();	
		StringTokenizer tokenizer;
		int numOfProp=0;
		String str;	
		Pattern pattern;
		int i;
		FileReader fileReader = new FileReader(fileName);
        BufferedReader bufReader = new BufferedReader(fileReader);
        
       
        
        while ((str = bufReader.readLine()) != null)
		{
        	if (str==null || str.equals("")) continue;
        	
        	tokenizer = new StringTokenizer(str,separator);
        	
        	if (tokenizer.countTokens()<0) 
        		throw new ClusteringException("invalid number of properties");
        	if (numOfProp==0)
        		numOfProp =  tokenizer.countTokens();
        	else 
        		{
        			if (numOfProp != tokenizer.countTokens())        		
        				throw new ClusteringException("all patterns must have same number of properties");        			
        		}
        			
        	pattern = new Pattern();        	
        	pattern.properties = new Vector<Double>(numOfProp);
        	pattern.labels = new Vector<Character>(numOfProp);
        	i=0;
        	while (tokenizer.hasMoreTokens())
        	{
        		String token = tokenizer.nextToken();
        		try{pattern.properties.add(i,new Double(token));}
        		catch(Exception e){pattern.properties.add(i,null);}
        		pattern.labels.add(i, token.charAt(0));
        		i++;
        	}
        	
        	temp.add(pattern);        
		}
		
		return temp;
	}

	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	public static Vector<CategoricalPropertyDefinition> getCategoricalDataHeader(BufferedReader bufReader, String separator1, String separator2) throws Exception
	{
		Vector<CategoricalPropertyDefinition> temp = new Vector<CategoricalPropertyDefinition>();	
		StringTokenizer tokenizer;
		int i, numOfProp=0;
		String str;		
	 
        str = bufReader.readLine();
        if ((numOfProp = Integer.parseInt(str.trim()))<0) 
        	throw new ClusteringException("invalid number of properties");
        
        while (numOfProp>0)
		{        	
        	if ((str=bufReader.readLine())==null) break;
        	
        	if (str.equals("")) continue;
        	
        	tokenizer = new StringTokenizer(str,separator1);
        	CategoricalPropertyDefinition cat = new CategoricalPropertyDefinition();
        	        	
        	if (tokenizer.countTokens()<0) 
        		throw new ClusteringException("invalid number of properties");
        	
        	cat.propertyName = tokenizer.nextToken();
        	tokenizer = new StringTokenizer(tokenizer.nextToken(),separator2);
        	
        	        			
        	      	
        	cat.categories = new HashMap<String,Integer>();
        	        	
        	i=0;
        	while (tokenizer.hasMoreTokens())
        	{
        		String token = tokenizer.nextToken();
        		try{cat.categories.put(token,i);}
        		catch(Exception e){cat.categories.put(null,i);}        		
        		i++;
        	}
        	
        	temp.add(cat);
        	numOfProp--;
        
		}		
		
		
		return temp;
	}

	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
		
	public static Vector<Pattern> getCategoricalData(BufferedReader bufReader, Vector<CategoricalPropertyDefinition> properties) throws Exception
	{
		String space = ",";		// default separator is ","
				
		return getCategoricalData(bufReader, properties, space);
	}
	
	/*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*/
	
	
	public static Vector<Pattern> getCategoricalData(BufferedReader bufReader, Vector<CategoricalPropertyDefinition> properties, String separator) throws Exception
	{
		Vector<Pattern> temp = new Vector<Pattern>();	
		StringTokenizer tokenizer;
		int numOfProp=properties.size();
		String str;	
		Pattern pattern;
		int i;
		       
        
        while ((str = bufReader.readLine()) != null)
		{
        	if (str==null || str.equals("")) continue;
        	
        	tokenizer = new StringTokenizer(str,separator);
        	
        	if (tokenizer.countTokens()<0) 
        		throw new ClusteringException("invalid number of properties");

        	if (numOfProp > tokenizer.countTokens())        		
        		throw new ClusteringException("all patterns must have same number of properties");        			
        	
       		
        	pattern = new Pattern();        	
        	pattern.properties = new Vector<Double>(numOfProp);
        	pattern.labels = new Vector<Character>(numOfProp);
        	i=0;
        	while (tokenizer.hasMoreTokens() && i<numOfProp)
        	{
        		String token = tokenizer.nextToken();
        		try{pattern.properties.add(i,new Integer(properties.get(i).categories.get(token)));}
        		catch(Exception e){e.printStackTrace();pattern.properties.add(i,null);}
        		pattern.labels.add(i, token.charAt(0));
        		i++;
        	}
        	if (tokenizer.hasMoreTokens())
        	{        		
        		try{pattern.correctCluster = tokenizer.nextToken();}
        		catch(Exception e){e.printStackTrace();}
        	}
        	temp.add(pattern);
        
		}
		
		
		
		return temp;
	}

}
