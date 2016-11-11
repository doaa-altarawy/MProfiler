package org.biojava.bio.constants;

import java.util.Vector;

public enum Dataset
{
	unknown(OrganismCode.unknown),
	dm01g(OrganismCode.DM),
	dm02g(OrganismCode.DM),
	dm03g(OrganismCode.DM),
	dm04g(OrganismCode.DM),
	dm05g(OrganismCode.DM),
	dm06g(OrganismCode.DM),
	dm07g(OrganismCode.DM),
	dm08g(OrganismCode.DM),
	mus01g(OrganismCode.MM),
	mus02g(OrganismCode.MM),
	mus03g(OrganismCode.MM),
	mus04g(OrganismCode.MM),
	mus05g(OrganismCode.MM),
	mus06g(OrganismCode.MM),
	mus07g(OrganismCode.MM),
	mus08g(OrganismCode.MM),
	mus09g(OrganismCode.MM),
	mus10g(OrganismCode.MM),
	mus11g(OrganismCode.MM),
	mus12g(OrganismCode.MM),
	//hm50g(OrganismCode.HS),
	hm01g(OrganismCode.HS),
	hm02g(OrganismCode.HS),
	hm03g(OrganismCode.HS),
	hm04g(OrganismCode.HS),
	hm05g(OrganismCode.HS),
	hm06g(OrganismCode.HS),
	hm07g(OrganismCode.HS),
	hm08g(OrganismCode.HS),
	hm09g(OrganismCode.HS),
	hm10g(OrganismCode.HS),
	hm11g(OrganismCode.HS),
	hm12g(OrganismCode.HS),
	hm13g(OrganismCode.HS),
	hm14g(OrganismCode.HS),
	hm15g(OrganismCode.HS),
	hm16g(OrganismCode.HS),
	hm17g(OrganismCode.HS),
	hm18g(OrganismCode.HS),
	hm19g(OrganismCode.HS),
	hm20g(OrganismCode.HS),
	hm21g(OrganismCode.HS),
	hm22g(OrganismCode.HS),
	hm23g(OrganismCode.HS),
	hm24g(OrganismCode.HS),
	hm25g(OrganismCode.HS),
	hm26g(OrganismCode.HS),
	yst01g(OrganismCode.SC),
	yst02g(OrganismCode.SC),
	yst03g(OrganismCode.SC),
	yst04g(OrganismCode.SC),
	yst05g(OrganismCode.SC),
	yst06g(OrganismCode.SC),
	yst07g(OrganismCode.SC),
	yst08g(OrganismCode.SC),
	yst09g(OrganismCode.SC),
	yst10g(OrganismCode.SC),


	dm01m(OrganismCode.DM),
	dm02m(OrganismCode.DM),
	dm03m(OrganismCode.DM),
	dm04m(OrganismCode.DM),
	dm05m(OrganismCode.DM),
	dm06m(OrganismCode.DM),
	dm07m(OrganismCode.DM),
	dm08m(OrganismCode.DM),
	mus01m(OrganismCode.MM),
	mus02m(OrganismCode.MM),
	mus03m(OrganismCode.MM),
	mus04m(OrganismCode.MM),
	mus05m(OrganismCode.MM),
	mus06m(OrganismCode.MM),
	mus07m(OrganismCode.MM),
	mus08m(OrganismCode.MM),
	mus09m(OrganismCode.MM),
	mus10m(OrganismCode.MM),
	mus11m(OrganismCode.MM),
	mus12m(OrganismCode.MM),
	hm01m(OrganismCode.HS),
	hm02m(OrganismCode.HS),
	hm03m(OrganismCode.HS),
	hm04m(OrganismCode.HS),
	hm05m(OrganismCode.HS),
	hm06m(OrganismCode.HS),
	hm07m(OrganismCode.HS),
	hm08m(OrganismCode.HS),
	hm09m(OrganismCode.HS),
	hm10m(OrganismCode.HS),
	hm11m(OrganismCode.HS),
	hm12m(OrganismCode.HS),
	hm13m(OrganismCode.HS),
	hm14m(OrganismCode.HS),
	hm15m(OrganismCode.HS),
	hm16m(OrganismCode.HS),
	hm17m(OrganismCode.HS),
	hm18m(OrganismCode.HS),
	hm19m(OrganismCode.HS),
	hm20m(OrganismCode.HS),
	hm21m(OrganismCode.HS),
	hm22m(OrganismCode.HS),
	hm23m(OrganismCode.HS),
	hm24m(OrganismCode.HS),
	hm25m(OrganismCode.HS),
	hm26m(OrganismCode.HS),
	yst01m(OrganismCode.SC),
	yst02m(OrganismCode.SC),
	yst03m(OrganismCode.SC),
	yst04m(OrganismCode.SC),
	yst05m(OrganismCode.SC),
	yst06m(OrganismCode.SC),
	yst07m(OrganismCode.SC),
	yst08m(OrganismCode.SC),
	yst09m(OrganismCode.SC),
	yst10m(OrganismCode.SC),
	
	dm01r(OrganismCode.DM),
	dm02r(OrganismCode.DM),
	dm03r(OrganismCode.DM),
	dm04r(OrganismCode.DM),
	dm05r(OrganismCode.DM),
	dm06r(OrganismCode.DM),
	dm07r(OrganismCode.DM),
	dm08r(OrganismCode.DM),
	mus01r(OrganismCode.MM),
	mus02r(OrganismCode.MM),
	mus03r(OrganismCode.MM),
	mus04r(OrganismCode.MM),
	mus05r(OrganismCode.MM),
	mus06r(OrganismCode.MM),
	mus07r(OrganismCode.MM),
	mus08r(OrganismCode.MM),
	mus09r(OrganismCode.MM),
	mus10r(OrganismCode.MM),
	mus11r(OrganismCode.MM),
	mus12r(OrganismCode.MM),
	hm01r(OrganismCode.HS),
	hm02r(OrganismCode.HS),
	hm03r(OrganismCode.HS),
	hm04r(OrganismCode.HS),
	hm05r(OrganismCode.HS),
	hm06r(OrganismCode.HS),
	hm07r(OrganismCode.HS),
	hm08r(OrganismCode.HS),
	hm09r(OrganismCode.HS),
	hm10r(OrganismCode.HS),
	hm11r(OrganismCode.HS),
	hm12r(OrganismCode.HS),
	hm13r(OrganismCode.HS),
	hm14r(OrganismCode.HS),
	hm15r(OrganismCode.HS),
	hm16r(OrganismCode.HS),
	hm17r(OrganismCode.HS),
	hm18r(OrganismCode.HS),
	hm19r(OrganismCode.HS),
	hm20r(OrganismCode.HS),
	hm21r(OrganismCode.HS),
	hm22r(OrganismCode.HS),
	hm23r(OrganismCode.HS),
	hm24r(OrganismCode.HS),
	hm25r(OrganismCode.HS),
	hm26r(OrganismCode.HS),
	yst01r(OrganismCode.SC),
	yst02r(OrganismCode.SC),
	yst03r(OrganismCode.SC),
	yst04r(OrganismCode.SC),
	yst05r(OrganismCode.SC),
	yst06r(OrganismCode.SC),
	yst07r(OrganismCode.SC),
	yst08r(OrganismCode.SC),
	yst09r(OrganismCode.SC),
	yst10r(OrganismCode.SC);
	
	
	OrganismCode organism;
	
	private Dataset(OrganismCode organism)
	{
		this.organism = organism;
	}
	
	/**
	 * Get Dataset of the given Str, if no element it returns null
	 * instead of throwing exception like valueOf(str)
	 * @param str
	 * @return
	 */
	public static Dataset getValue(String str)
	{		
		try{return valueOf(str);}
		catch (Exception e) 
		{return null;}
	}
	
	public static Dataset getValue(String str, DatasetType type)
	{		
		try{return valueOf(str);}
		catch (Exception e) 
		{			
			try{return valueOf(str+type.code);}	
			catch (Exception e2) {return null;}
		}
	}
	
	public static Vector<Dataset> getByOrganism(OrganismCode org)
	{
		Vector<Dataset> datasets = new Vector<Dataset>();
		
		for (int i=0; i<values().length; i++)
		{
			if (values()[i].organism.equals(org))
				datasets.add(values()[i]);
		}
		
		return datasets;
	}
	
	public static Vector<Dataset> getByOrganismAndType(OrganismCode org, DatasetType type)
	{
		Vector<Dataset> datasets = new Vector<Dataset>();
		
		for (int i=0; i<values().length; i++)
		{
			if ((values()[i].organism.equals(org))
					&& (values()[i].name().charAt(values()[i].name().length()-1) == type.code))
				datasets.add(values()[i]);
		}
		
		return datasets;
	}
	public static Vector<Dataset> getByDatasetType(DatasetType type)
	{
		Vector<Dataset> datasets = new Vector<Dataset>();
		
		if (type.equals(DatasetType.All))
		{
			for (int i=0; i<values().length; i++)			
				datasets.add(values()[i]);
		}
		else
		for (int i=0; i<values().length; i++)
		{
			if (values()[i].name().charAt(values()[i].name().length()-2) == type.code)
				datasets.add(values()[i]);
		}
		
		return datasets;
	}
	
	public DatasetType getDatasetType()
	{
		for (int i =0; i< DatasetType.values().length; i++)
		{
				char c = name().charAt(name().length()-1);
				if (c == DatasetType.values()[i].code)
					return DatasetType.values()[i];
		}
		return null;
	}
}
