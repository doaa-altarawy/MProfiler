package org.biojava.bio.constants;

import java.awt.Color;

public enum MotifFinder
{
	//BioProspector (Color.red),
	MEME (Color.yellow), 
	MDSCAN (Color.lightGray), 
	WEEDER (Color.white), 
	BP (Color.red), 
	MTSAMP (Color.cyan),
	SPACE (Color.gray),
	MITRA (Color.magenta),
	ANNSPEC (Color.darkGray), ////////////redundant color find another
	ALIGNACE (Color.blue), ////////////redundant color find another
	IMPROBIZER (Color.blue), ////////////redundant color find another
	MotifCut (Color.darkGray),		// not used
	MProfiler (Color.pink),
	Doaa (Color.pink),
	Target (Color.green),			// the correct or reference motif
	MotifVoter (Color.black);

	
	public Color color;
	
	private MotifFinder(Color color)
	{
		this.color = color;
	}
}
