package org.biojava.bio.alignment;

import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;

/*
 * Created on 13.06.2008 
 */

/**
 * extends Smith and Waterman efficient dynamic programing algorithm to
 * perform local sequence alignments, to give more usable output
 * 
 * @author Doaa Altarawy 
 */
public class LocalAlignment extends SmithWaterman 
{ 

  /**
   * Constructs the new SmithWaterman alignment object. Alignments are only
   * performed, if the alphabet of the given <code>SubstitutionMatrix</code>
   * equals the alpabet of both the query and the target <code>Sequence</code>.
   * The alignment parameters here are expenses and not scores as they are in
   * the <code>NeedlemanWunsch</code> object. scores are just given by
   * multipliing the expenses with <code>(-1)</code>. For example you could
   * use parameters like "-2, 5, 3, 3, 0". If the expenses for gap extension are
   * equal to the cost of starting a gap (delete or insert), no affine gap
   * penalties are used, which saves memory.
   * 
   * @param match
   *          expenses for a match
   * @param replace
   *          expenses for a replace operation
   * @param insert
   *          expenses for a gap opening in the query sequence
   * @param delete
   *          expenses for a gap opening in the target sequence
   * @param gapExtend
   *          expenses for the extension of a gap which was started earlier.
   * @param matrix
   *          the <code>SubstitutionMatrix</code> object to use.
   */
  public LocalAlignment(double match, double replace, double insert, double delete, double gapExtend, SubstitutionMatrix matrix) 
  {
    super(match, replace, insert, delete, gapExtend, matrix);    
  }
  
  /**
   * 
   * @param match
   * @param replace
   * @param insertDelete the same cost for inset and delete parameters
   * @param gapExtend
   * @param matrix
   */
  public LocalAlignment(double match, double replace, double insertDelete, double gapExtend, SubstitutionMatrix matrix) 
  {
    super(match, replace, insertDelete, insertDelete, gapExtend, matrix);    
  }
  
  /**
   * 
   * @param match
   * @param replace
   * @param insertDelete the same cost for inset and delete parameters
   * @param gapExtend
   * uses unity SubstitutionMatrix, with -1 for mismatch
   */
  public LocalAlignment(double match, double replace, double insertDelete, double gapExtend) 
  {
	  super(match, replace, insertDelete, insertDelete, gapExtend, getSubstitutionMatrix());    
  }
  
  private static SubstitutionMatrix getSubstitutionMatrix()
  {
	  Alphabet alphabet = AlphabetManager.alphabetForName("ACGT");
	  return new SubstitutionMatrix((FiniteAlphabet)alphabet, 1, -1);	    
  }
  
  public int getTargetStart()
  {
	  return super.targetStart;
  }
 
  
  public int getQueryStart()
  {
	  return super.queryStart;
  }
  
}
