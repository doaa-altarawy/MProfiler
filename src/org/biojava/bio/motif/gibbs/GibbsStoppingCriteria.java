package org.biojava.bio.motif.gibbs;
 
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionTools;
 
 
/**
 * Defines the criteria under which Gibbs Sampling should stop
 */
public interface GibbsStoppingCriteria {
 
  /**
   * Uses a heuristic procedure to determine when to stop. If the information
   * content of the motif has failed to increase above its previous maximum for
   * 100 iterations then the method will return true. NOTE: it is expected that
   * the same SimpleGibbsSampler will be passed to the stop() method at each
   * call.
   */
  public static GibbsStoppingCriteria HEURISTIC = new Heuristic();
 
  /**
   * Returns true when the emission spectra of the last iteration equals that
   * of this iteration. Note that this may never return if convergence is not
   * reached. Thus the method has a built in stopping point of 10,000
   * iterations. NOTE: it is expected that the same SimpleGibbsSampler will be
   * passed to the stop() method at each call.
   */
  public static GibbsStoppingCriteria CONVERGE = new Converge();
 
 /**
  * This method should return true when stopping criteria have been reached.
  * @param sga the GibbsAligner that is being tested for stopping conditions
  * @return true if it should stop, false otherwise.
  */
  public boolean stop(SimpleGibbsAligner sga);
 
  /**
   * Implementation of GibbsStoppingCriteria
   */
  class Heuristic implements GibbsStoppingCriteria{
    double bestInfo = 0.0; //the level of conservation
    int bestIteration = 0; //the most conserved pattern
 
    public boolean stop(SimpleGibbsAligner sga){
      double info = sga.getInfoContent();
      if(info > bestInfo){
        bestInfo = info;
        bestIteration = sga.getIterations();
        return false; //don"t stop
      }else if(sga.getIterations() >= bestIteration+99){
        return true;
      }
      return false; //don"t stop
    }
  }// end of Heuristic
 
  /**
   * Implementation of GibbsStoppingCriteria
   */
  class Converge implements GibbsStoppingCriteria{
    Distribution[] previous = null; //the last pattern
 
    public boolean stop(SimpleGibbsAligner sga){
      if(previous == null) return false; //there is no previous yet.
      if(sga.getIterations() == 10000) return true; //max iterations.
      try{
        if (DistributionTools.areEmissionSpectraEqual(previous,sga.getPattern())){
          return true; // patterns have converged.
        }
        else {
          previous = sga.getPattern();
          return false; //don"t stop
        }
      }catch(BioException e){
        //this can"t really happen but...
        e.printStackTrace();
        return false;
      }
    }
  }// end of converge
 
 
}// end of GibbsStoppingCriteria 