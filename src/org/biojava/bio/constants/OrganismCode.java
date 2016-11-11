package org.biojava.bio.constants;

public enum OrganismCode
{
	HS ("Homo sapiens", "HS"),
	HSI	("Homo sapiens (intergenic frequencies)", "HSI"),
	MM 	("Mus musculus", "MM"),
	MMI ("Mus musculus (intergenic frequencies)", "MMI"),
	CFI ("Canis familiaris (intergenic frequencies)", "CFI"),
	RN	("Rattus norvegicus", "RN"),
	DM ("Drosophila melanogaster", "DM"),
	SC	("S. cerevisiae", "SC"),
	FR ("(Taki)Fugu rubripes", "FR"),
	DR ("Danio rerio", "DR"),
	AT ("Arabidopsis thaliana", "AT"),
	GG ("Gallus gallus", "GG"),
	CE ("C. elegans", "CE"),
	AG ("Anopheles gambiae", "AG"),
	AN ("A.nidulans", "AN"),
	CI ("Ciona intestinalis", "CI"),
	XT ("Xenopus tropicalis", "XT"),
	BEC ("E.coli K12", "BEC"),
	BBS ("Bacillus subtilis", "BBS"),
	unknown("Unknown", "");
	
	private final String nameText;
	private final String weederCode;
	
	OrganismCode(String name, String weederCode)
	{
		this.nameText = name;
		this.weederCode = weederCode;
	}

	public String nameText() {return nameText;}
	public String weederCode() {return weederCode;}
	
	public static String[] getNames()
	{
		String names[] = new String[values().length];
		
		for (int i=0; i<values().length; i++)
		{
			names[i] = values()[i].nameText;
		}
		return names;
	}
	
	public static OrganismCode getOrganism(int i)
	{
		return values()[i];
	}
}
