package clustering.techniques.denistyBased;

import java.util.List;

import org.biojava.bio.motif.Motif;

public interface Clusterer
{
	public void buildClusterer(List<Motif> motifs) throws Exception;
	public int numberOfClusters();
	public List<Motif> getCluster(int i);
	public String toString();
	    
}
