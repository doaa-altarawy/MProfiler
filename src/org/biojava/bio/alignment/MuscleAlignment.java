package org.biojava.bio.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.Vector;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.AbstractSymbolList;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleAlignment;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.bio.seq.RichSequence.IOTools;

/**
 * Performs Muscle multiple alignment using an external exe file
 * @author Doaa1
 *
 */
public class MuscleAlignment extends AbstractSymbolList implements Alignment
{
	private LinkedHashMap labelToSymbolList;
	  private List labels;
	  private Alphabet alphabet;
	  private int length;

	
	  /**
	   * Generate an alignment from a list of SymbolLists.
	   * <p>
	   * *
	   * @param labelToResList  the label-to-symbol list mapping
	 * @throws Exception 
	   */
	  public MuscleAlignment(Vector<SymbolList> list, String n)  throws Exception 
	  {
	    if(list.isEmpty()) 
	       	return;
	    

	    Map<Integer, SymbolList> labelToResList = getMuscleExternalRun(list, n);
	    this.labels = Collections.unmodifiableList(new ArrayList(labelToResList.keySet()));
	    this.labelToSymbolList =new LinkedHashMap( labelToResList);

	    
	    
	    int length = -1;
	    List alphaList = new ArrayList();
	    for(Iterator li = labels.iterator(); li.hasNext(); ) 
	    {
	      Object label = li.next();
	      try 
	      {
	        SymbolList rl = symbolListForLabel(label);
	        alphaList.add(rl.getAlphabet());
	        if(length == -1) 
	          length = rl.length();
	        else 
	        {
	          if(rl.length() != length) 
	          {
	            StringBuffer sb = new StringBuffer();
	            for(Iterator labI = labels.iterator(); labI.hasNext(); ) 
	            {
	              Object lab = labI.next();
	              sb.append("\n\t" + lab + " (" + symbolListForLabel(lab).length() + ")");
	            }
	            throw new IllegalArgumentException(
	              "All SymbolLists must be the same length: " + sb.substring(0)
	            );
	          }
	        }
	      } 
	      catch (NoSuchElementException nsee) 
	      {
	        if(labelToSymbolList.containsKey(label))
	          throw new IllegalArgumentException("The symbol list associated with " + label + " is null");
	        else 
	          throw new BioError("Something is screwey - map is lying about key/values", nsee);	        
	      }
	    }

	    this.alphabet = AlphabetManager.getCrossProductAlphabet(alphaList);
	    this.length = length;
	  }
	  
	  private Map<Integer, SymbolList> getMuscleExternalRun(Vector<SymbolList> seq, String n) throws Exception
	  {
		 		
		String line;
		File dir = new File (".\\resources\\muscle\\");		
		File in = new File(dir+ "\\" + n + "_input.muscle");
		File out = new File(dir+ "\\" + n + "_output.muscle");
		String command = ".\\resources\\muscle\\muscle -stable -in "
			+ in.getName()+	" -out "+ out.getName();
			
		FileWriter f = new FileWriter(in);
		
		for (Iterator<SymbolList> itr=seq.iterator(); itr.hasNext();)
			f.append(">s\n"+itr.next().seqString()+"\n");
		f.flush();
		f.close();
				
		Process p = Runtime.getRuntime().exec(command, null, dir);
		
		BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	    while ((line = error.readLine()) != null){}
	       //System.out.println(line);
	    
		
	    SequenceIterator itr = IOTools.readFastaDNA(new BufferedReader(new FileReader(out)), null);		
		Map<Integer, SymbolList> list = new HashMap<Integer, SymbolList>();
		for (int i=0; itr.hasNext(); i++) 
		    list.put(i, itr.nextSequence()); 
		
		in.delete();
		out.deleteOnExit();	 	
	 	
	 	return list;
	  }
	  /*========================================================================*/
	  
	  public int length() {
	    return length;
	  }

	  public Alphabet getAlphabet() {
	    return alphabet;
	  }

	  public Symbol symbolAt(int index) {
	    try {
	      if (labels.size() == 1) {
	          return symbolAt(labels.get(0), index);
	      } else {
	          return alphabet.getSymbol(new ColAsList(index));
	      }
	    } catch (IllegalSymbolException ire) {
	      throw new BioError(

	        "Somehow my crossproduct alphabet is incompatible with column " + index, ire
	      );
	    }
	  }

	  public List getLabels() {
	    return labels;
	  }

	  public Symbol symbolAt(Object label, int column) {
	    return symbolListForLabel(label).symbolAt(column);
	  }

	  public Alignment subAlignment(Set labels, Location loc)
	  throws NoSuchElementException {
	    Map labelsToResList = new LinkedHashMap();
	    Iterator i;
	    if(labels != null) {
	      i = labels.iterator();
	    } else {
	      i = getLabels().iterator();
	    }
	    while(i.hasNext()) {
	      Object label = i.next();
	      SymbolList sym = symbolListForLabel(label);
	      if(loc != null) {
	        sym = loc.symbols(sym);
	      }
	      labelsToResList.put(label, sym);
	    }
	    return new SimpleAlignment(labelsToResList);
	  }

	  public SymbolList symbolListForLabel(Object label)
	  throws NoSuchElementException {
	    SymbolList rl = (SymbolList) labelToSymbolList.get(label);
	    if(rl == null) {
	      throw new NoSuchElementException("No symbol list associated with label " + label);
	    }
	    return rl;
	  }

	  public Iterator symbolListIterator() {
	    return new Alignment.SymbolListIterator(this);
	  }

	  /**
	   * Makes a column of the alignment behave like a list.
	   *
	   * @author Matthew Pocock
	   */
	  private final class ColAsList extends AbstractList implements Serializable {
	    private final int col;

	    public ColAsList(int col) {
	      this.col = col;
	    }

	    protected ColAsList() {
	      this.col = 0;
	    }

	    public Object get(int indx) {
	      return symbolAt(labels.get(indx), col);
	    }

	    public int size() {
	      return labels.size();
	    }
	  }
}
