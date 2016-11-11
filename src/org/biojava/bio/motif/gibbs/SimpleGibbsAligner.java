package org.biojava.bio.motif.gibbs;
 
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Vector;
import java.util.logging.Logger;


import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.SimpleAlignment;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
 
 
 
/**
 * A class that uses Gibbs Sampling to generate a local alignment of an over
 * represented motif.
 */
public class SimpleGibbsAligner 
{
  private static Logger log = Logger.getLogger(SimpleGibbsAligner.class.getName());
  private Sequence[] s; // sequence array.
  private int w; //window size.
  private int[] a; //starting indices.
  private int iterations = 0;
  private Distribution[] pattern; //the probabilistic pattern description.
  private Distribution background; //the probabilistic background description.
  private Random rand; //random number generator
  private Alphabet alphabet; //the alphabet in which the sampler operates.
  private GibbsStoppingCriteria criteria; //determines when to stop sampling.
 
  /**
   * Constructs the gibbs aligner to find a common motif in a collection
   * of sequences. It is assumed that all the sequences are constructed
   * from the same <code>Alphabet</code>. If this is not the case then calls
   * to iterate will throw exceptions. This class is designed to be single use
   * and is not thread safe. To use in a threaded environment each thread
   * should be given its own SimpleGibbsAligner.
   *
   * @param windowSize the expected size of the motif
   * @param it a collection of sequences in which to search for a motif.
   * @param criteria an object which specifies when sampling should stop.
   */
  public SimpleGibbsAligner(int windowSize,
                            SequenceIterator it,
                            GibbsStoppingCriteria criteria){
    w = windowSize;
    this.criteria = criteria;
    rand = new Random();
 
    //get the sequences
    Vector v = new Vector();
    while(it.hasNext()){
      try{
        v.add(it.nextSequence());
      }catch(BioException e){
        //cannot retrieve the sequence from the iterator, not likely to happen.
        e.printStackTrace();
      }
    }
    v.trimToSize();
    s = new Sequence[v.size()];
    v.copyInto(s);
 
    //Initialize the offsets
    a = new int[s.length];
    a = initIndices();
 
    //set the alphabet
    alphabet = s[0].getAlphabet();
  }
 
 
 
  /**
   * Initialize an array of random offsets.
   * @return the array of offsets
   */
  private int[] initIndices(){
    int[] indices = new int[s.length];
    for (int i = 0; i < indices.length; i++) {
      int index = rand.nextInt(s[i].length() - w-1);
      // as we are making offset indices to symbol lists
      // they must be from 1 not 0
      index++;
      indices[i] = index;
    }
    return indices;
  }
 
  /**
   * Iterates through a procedure of predictive updates and sampling until
   * the stopping criteria defined in the <code>stop()</code> method are met.
   * Once the method returns the <code>getXXX</code> methods can be used to
   * determine the results.
   */
  public void iterate(){
    try {
      //choose a sequence at random
      int index = rand.nextInt(s.length);
      do{
        //calculate pattern in all but the chosen sequence
        pattern = updatePattern(index, a);
        //Occasional try a phase shift
        if(rand.nextDouble() < 0.1){
          tryPhaseShift(index);
        }
        //calculate the background
        background = updateBackground(index);
        //sample the randomly chosen sequence to find the best start index a.
        a[index] = sampleSequence(index);
        //reportMatch(a[index], s[index]);
        iterations++;
        index = (++index)%s.length;
      }while(stop() == false);
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
  }
 
  /**
   * Determines when to stop iterating.
   * @return true if the StoppingCriteria says to stop and false otherwise.
   */
  protected boolean stop(){
    return criteria.stop(this);
  }
 
  /**
   * Produces a pattern to describe the motif covered by the window
   * @param excludeIndex the index of the sequence to be excluded from sampling.
   * @param offsets the matrix of offset positions
   * @return the updated motif pattern
   */
  private Distribution[] updatePattern(int excludeIndex, int[] offsets){
    Distribution[] d = null;
    StringBuffer buff = new StringBuffer("Sequances before Alignment:\n");
    
    Map label2Res = new HashMap(s.length);
    for (int i = 0; i < s.length; i++) 
    {//for each sequence
      if(i == excludeIndex) continue; //except this sequence
      SymbolList subSeq = s[i].subList(offsets[i],
                                       offsets[i] +w -1);//take the subsequence
      label2Res.put(new Integer(i),subSeq); //put it in the hashmap
      buff.append("Sequence "+i+": " +subSeq.seqString()+"\n");
    }
    Alignment al = new SimpleAlignment(label2Res);//make an alignment of subseqs (seq's don't change!!!!!)
    
    /* logging the aligned matrix*/    
    buff.append("Sequances After Alignment:\n");
    for (Iterator itr= al.getLabels().iterator(); itr.hasNext();)
    {
    	Object label = itr.next();
    	String seq = (String)al.symbolListForLabel(label).seqString();
    	buff.append("Sequence "+label+": " +seq + "\n");
    }    
    
    try {
      d = DistributionTools.distOverAlignment(al, false,1.0);//make the pattern
    }
    catch (IllegalAlphabetException ex) {
      ex.printStackTrace();
    }
 
    log.info(buff.toString());
    return d;
  }
 
  /**
   * produces a distribution to describe the background distribution
   * @param excludeIndex the index of the sequence to exclude
   * @return the updated background distribution.
   */
  private Distribution updateBackground(int excludeIndex){
    Distribution d = null;
 
    try {
      DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
      d = DistributionFactory.DEFAULT.createDistribution(alphabet);
      dtc.setNullModelWeight(1.0);
      dtc.registerDistribution(d);
 
      for (int i = 0; i < s.length; i++) {//for each sequence
        if(i == excludeIndex) continue; //except this sequence
        for(int j = 1; j <= s[i].length(); j++){//count each base
          if(j >= a[i] && j < a[i] + w-1) continue; //except these ones
          dtc.addCount(d, s[i].symbolAt(j), 1.0);
        }
      }
      dtc.train();
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
    return d;
  }
 
  /**
   * Attempts to prevent the pattern getting locked in a local optimum by
   * shifting the pattern one step to the left or right and seeing if it is
   * better than the current pattern. If the phase shift improves the model
   * the pattern and offsets will be updated.
   * @param excludeIndex the index of the sequence to be excluded.
   */
  private void tryPhaseShift(int excludeIndex){
    int[] newOffSets = new int[a.length];
    System.arraycopy(a,0,newOffSets,0,a.length); // copy offsets
    Distribution[] newPattern;
 
    if (rand.nextBoolean()) {//shift left
      for (int i = 0; i < newOffSets.length; i++) {
        if(i == excludeIndex) continue; //skip this sequence
        if(newOffSets[i] > 1) newOffSets[i]--;
      }
    }
    else {// shift right
      for (int i = 0; i < newOffSets.length; i++) {
        if(i == excludeIndex) continue; //skip this sequence
        if(newOffSets[i] < s[i].length() - w-2) newOffSets[i]++;
      }
    }
 
    newPattern = updatePattern(excludeIndex, newOffSets);
    if(getInfoContent(newPattern) > getInfoContent(pattern)){
      a = newOffSets;
      pattern = newPattern;
    }
  }
 
  /**
   * Determines a weighted distribution of offsets in the sequence to be
   * sampled and randomly selects an offset from that distribution to be used
   * in the next pattern update.
   * @param sequenceIndex the sequence to be sampled.
   * @return the selected offset
   */
  private int sampleSequence(int sequenceIndex){
    Distribution d = null;
    try {
      SymbolList seq = s[sequenceIndex];
      //make an alphabet of the possible offsets
      IntegerAlphabet.SubIntegerAlphabet alpha =
             IntegerAlphabet.getSubAlphabet(1, seq.length()-w-1);
      //make a distribution to hold the weighted probabilities of each offset.
      d = DistributionFactory.DEFAULT.createDistribution(alpha);
      DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();
      dtc.setNullModelWeight(1.0);
      dtc.registerDistribution(d);
 
      //score each subsequence
      for(int i = 1; i <= seq.length()-w-1; i++){
        double score = scoreSequence(seq.subList(i, i+w-1));
        //add the weight to the distribution of offsets
        dtc.addCount(d,alpha.getSymbol(i),score);
      }
      dtc.train();
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
 
    //sample the distribution of offsets
    int offset = ((IntegerAlphabet.IntegerSymbol)d.sampleSymbol()).intValue();
    return offset;
  }
 
  /**
   * Scores a potential motif against the pattern description and background
   * distribution.
   * @param sl the potential motif to score
   * @return the score
   */
  private double scoreSequence(SymbolList sl){
    double pMotif = 1.0;
    double pBackGround = 1.0;
 
    for(int i = 0; i < sl.length(); i++){
      Symbol s = sl.symbolAt(i+1); //+1 as we are indexing from zero this time
      try {
        pMotif *= pattern[i].getWeight(s); //probability of s at position i
        pBackGround *= background.getWeight(s); //probability of s in background
      }
      catch (IllegalSymbolException ex) {
        ex.printStackTrace();
      }
    }
    return pMotif/pBackGround;
  }
 
  /**
   * Determines the information content (in bits) of the motif including pseudo
   * counts.
   * @return the Information content.
   */
  public double getInfoContent(){
    return getInfoContent(pattern);
  }
 
  /**
   * determines the information content (in bits) of the specified pattern
   * including pseudo counts.
   * @param d the pattern of the motif
   * @return the information content
   */
  private double getInfoContent(Distribution[] d){
    double info = 0.0;
    for (int i = 0; i < d.length; i++) {
      info += DistributionTools.bitsOfInformation(d[i]);
    }
    return info;
  }
 
  
  /**
   * Returns the current <code>Alphabet</code> being used.
   * @return an <CODE>Alphabet</CODE>
   */
  public Alphabet getAlphabet(){
    return alphabet;
  }
 
  /**
   * Get the background distribution.
   * @return a <CODE>Distribution</CODE> of background frequencies.
   */
  public Distribution getBackground() {
    return background;
  }
 
  /**
   * The current iteration of the sampler
   * @return an int >= 0
   */
  public int getIterations() {
    return iterations;
  }
 
  /**
   * The current pattern at this iteration of the sampler
   * @return the pattern as a <CODE>Distribution[]</CODE>. 
   * Effectively a weight matrix.
   */
  public Distribution[] getPattern() {
    return pattern;
  }
 
  /**
   * Tje set of sequence offsets being used for this iteration of 
   * sampling
   * @return an array of ints &ge; 1
   */
  public int[] getOffSets(){
    return a;
  }
 
  /**
   * The set of <code>Sequence</code>s being sampled
   * @return  a <CODE>Sequence[]</CODE>
   */
  public Sequence[] getSequences(){
    return s;
  }
 
  /**
   * The size of the pattern being sampled for.
   * @return  an <code>int</code> &gt; 0
   */
  public int getWindowSize(){
    return w;
  }
}