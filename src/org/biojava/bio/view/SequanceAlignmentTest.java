/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 */


package org.biojava.bio.view;

import java.util.Iterator;

import org.biojava.bio.alignment.LocalAlignment;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.gui.WMPanel;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.SymbolList;

/* @author Doaa Altarawy */
public class SequanceAlignmentTest 
{

  
    public static void main(String[] args) throws Exception
    {
    	
    	    	
    	Sequence seq1 = DNATools.createDNASequence("agctgagcttttggcc", "Seq  ");
    	Sequence seq2 = DNATools.createDNASequence("gactt", "Motif");
    	
    	//Alphabet alpha = AlphabetManager.alphabetForName("ACGT");
    	FiniteAlphabet alpha = (FiniteAlphabet)seq1.getAlphabet();
    	
    	SubstitutionMatrix substitutionMatrix = new SubstitutionMatrix(alpha, 1, -1);
    	//substitutionMatrix.printMatrix();
    	LocalAlignment localAlign = new LocalAlignment(-2,1,1,1,substitutionMatrix);
    	Alignment align = localAlign.getAlignment(seq2, seq1);
    	
    	SimpleWeightMatrix wm = new SimpleWeightMatrix(alpha,8, DistributionFactory.DefaultDistributionFactory.DEFAULT);
    	Distribution dist = wm.getColumn(1);
    	
    	WMPanel.wmViewer(wm, "This is a test subMatrix");	    	
    	for (Iterator itr = align.getLabels().iterator(); itr.hasNext();)
    	{
    		String label = (String) itr.next();
    		SymbolList  s = align.symbolListForLabel(label);
    		System.out.println(label+": "+ s.seqString());
    	}
    	
    	System.out.println(localAlign.getAlignmentString());
      }

    
}
