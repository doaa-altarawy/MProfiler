package org.biojava.bio.constants;


public enum DatasetType
{
	Generic('g'),
	Marcov('m'),
	Real('r'),
	All('-');
	
	public char code;
	
	private DatasetType(char code)
	{
		this.code = code;
	}
	
}
