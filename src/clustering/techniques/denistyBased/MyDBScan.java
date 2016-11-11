package clustering.techniques.denistyBased;


import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif;
import weka.clusterers.forOPTICSAndDBScan.DataObjects.DataObject;
import weka.core.Utils;


public class MyDBScan implements Clusterer
{

	/**
     * Specifies the radius for a range-query
     */
    private double epsilon = 0.8;

    /**
     * Specifies the density (the range-query must contain at least minPoints DataObjects)
     */
    private int minPoints = 2;

    private int numberOfGeneratedClusters;
       
    private MyDatabase<Motif> database;
    
    private int clusterID;

    private double elapsedTime;
  
    public MyDBScan(){}
    
    public MyDBScan(int minPoints, double epsilon)
    {
    	this.epsilon = epsilon;
    	this.minPoints = minPoints;    	
    }
    
   /**
     * Generate Clustering via DBScan
     * @param instances The instances that need to be clustered
     * @throws java.lang.Exception If clustering was not successful
     */
    public void buildClusterer(List<Motif> motifs) throws Exception 
    {
       
        long time_1 = System.currentTimeMillis();
        
        numberOfGeneratedClusters = 0;
        clusterID = 0;
        
        database = new MyDatabase<Motif>(motifs);

        Iterator iterator = database.dataObjectIterator();
        while (iterator.hasNext()) 
        {
            DataObject dataObject = (DataObject) iterator.next();
            if (dataObject.getClusterLabel() == DataObject.UNCLASSIFIED) 
            {
                if (expandCluster(dataObject)) 
                {
                    clusterID++;
                    numberOfGeneratedClusters++;
                }
            }
        }

        long time_2 = System.currentTimeMillis();
        elapsedTime = (double) (time_2 - time_1) / 1000.0;
    }

    /**
     * Assigns this dataObject to a cluster or remains it as NOISE
     * @param dataObject The DataObject that needs to be assigned
     * @return true, if the DataObject could be assigned, else false
     */
    private boolean expandCluster(DataObject dataObject) 
    {
        List seedList = database.epsilonRangeQuery(getEpsilon(), dataObject);
        /** dataObject is NO coreObject */
        if (getNumOfFinders(seedList) < getMinPoints())//if (seedList.size() < getMinPoints()) 
        {
            dataObject.setClusterLabel(DataObject.NOISE);
            return false;
        }

        /** dataObject is coreObject */
        for (int i = 0; i < seedList.size(); i++) 
        {
            DataObject seedListDataObject = (DataObject) seedList.get(i);
            /** label this seedListDataObject with the current clusterID, because it is in epsilon-range */
            seedListDataObject.setClusterLabel(clusterID);
            if (seedListDataObject.equals(dataObject)) 
            {
                seedList.remove(i);
                i--;
            }
        }

        /** Iterate the seedList of the startDataObject */
        for (int j = 0; j < seedList.size(); j++) 
        {
            DataObject seedListDataObject = (DataObject) seedList.get(j);
            List seedListDataObject_Neighbourhood = database.epsilonRangeQuery(getEpsilon(), seedListDataObject);

            /** seedListDataObject is coreObject */
            if (getNumOfFinders(seedListDataObject_Neighbourhood) >= getMinPoints()) //if (seedListDataObject_Neighbourhood.size() >= getMinPoints()) 
            {
                for (int i = 0; i < seedListDataObject_Neighbourhood.size(); i++) 
                {
                    DataObject p = (DataObject) seedListDataObject_Neighbourhood.get(i);
                    if (p.getClusterLabel() == DataObject.UNCLASSIFIED || p.getClusterLabel() == DataObject.NOISE) {
                        if (p.getClusterLabel() == DataObject.UNCLASSIFIED) 
                        {
                            seedList.add(p);
                        }
                        p.setClusterLabel(clusterID);
                    }
                }
            }
            seedList.remove(j);
            j--;
        }

        return true;
    }

    private int getNumOfFinders(List<MotifDataObject> motifs)
	{		
		Set<MotifFinder> finders = new HashSet<MotifFinder>();
		
		for (Iterator<MotifDataObject> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next().getMotif();
			finders.add(m.getFinder());			// This is set, elements doesn't repeat				
		}
		
		return finders.size();
	}
    
    public int numberOfClusters(){
        return numberOfGeneratedClusters;
    }
 
    public void setMinPoints(int minPoints) {
        this.minPoints = minPoints;
    }
    
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public int getMinPoints() {
        return minPoints;
    }

    public MyDatabase getDatabase()
    {
    	return database;
    }
  
     /**
     * Returns a description of the clusterer
     * 
     * @return a string representation of the clusterer
     */
    public String toString() 
    {
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append("DBScan clustering results\n" +
                "========================================================================================\n\n");
        stringBuffer.append("Clustered DataObjects: " + database.size() + "\n");
        stringBuffer.append("Epsilon: " + getEpsilon() + "; minPoints: " + getMinPoints() + "\n");
        stringBuffer.append("Number of generated clusters: " + numberOfGeneratedClusters + "\n");
        DecimalFormat decimalFormat = new DecimalFormat(".##");
        stringBuffer.append("Elapsed time: " + decimalFormat.format(elapsedTime) + "\n\n");

        for (Iterator<MotifDataObject> itr = database.dataObjectIterator(); itr.hasNext();) 
        {
            DataObject dataObject = itr.next();
            stringBuffer.append(Utils.padRight(((MotifDataObject)dataObject).getMotif().getFinder()+", "+ ((MotifDataObject)dataObject).getMotif().getCorrect(), 69) + "  -->  " +
                    ((dataObject.getClusterLabel() == DataObject.NOISE) ?
                    "NOISE\n" : dataObject.getClusterLabel() + "\n"));
        }
        return stringBuffer.toString() + "\n";
    }
   
    public List<Motif> getCluster(int i)
    {
    	List<Motif> motifs = new ArrayList<Motif>();
    	
    	for (Iterator<MotifDataObject> itr = database.dataObjectIterator(); itr.hasNext();) 
        {
            DataObject dataObject = itr.next();
            if (dataObject.getClusterLabel() == i)
            	motifs.add(((MotifDataObject)dataObject).getMotif());
        }
    	
    	return motifs;
    }
}
