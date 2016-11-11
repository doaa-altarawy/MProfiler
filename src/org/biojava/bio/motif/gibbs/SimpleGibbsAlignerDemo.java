package org.biojava.bio.motif.gibbs; 
 
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
 
 
 
public class SimpleGibbsAlignerDemo 
{
	
  public void runDemo(File file, boolean protein, int window, int trials) 
  {
        
    SequenceIterator it;
 
    try
    {
	    for(int i = 0; i < trials; i++)
	    {
		      BufferedReader br = new BufferedReader(new FileReader(file));
		      if(protein)
		        it =(SequenceIterator)SeqIOTools.fileToBiojava("fasta", "protein", br);
		      else
		        it =(SequenceIterator)SeqIOTools.fileToBiojava("fasta", "DNA", br);
		      
		      
		      //make an aligner with Heuristic stopping criteria
		      SimpleGibbsAligner gibbs = new SimpleGibbsAligner(window,
		          it, GibbsStoppingCriteria.HEURISTIC);

		      
		      //System.out.println("Starting Information Content (bits): "+gibbs.getInfoContent());

		      //start the aligner running
		      gibbs.iterate();
		 
		      //how many iterations till convergence?
		      System.out.println("Converged after "+gibbs.getIterations()+" iterations");
		      //What is the information content of the motif?
		      System.out.println("Information Content (bits): "+gibbs.getInfoContent());
		      
		      //get the sequences, offsets and window size to print out the motif
		      Sequence[] seqs = gibbs.getSequences();
		      int[] offSets = gibbs.getOffSets();
		      int wind = gibbs.getWindowSize();
		 
		      //print out the motif
		      for (int j = 0; j < offSets.length; j++) {
		        System.out.println(seqs[j].subStr(offSets[j],offSets[j]+wind -1));
		      }
		      System.out.println();
	    	}
    	}
	    catch (Exception e) {
			e.printStackTrace();
		}
    }
  
  
}