/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    Copyright (C) 2004
 *    & Matthias Schubert (schubert@dbs.ifi.lmu.de)
 *    & Zhanna Melnikova-Albrecht (melnikov@cip.ifi.lmu.de)
 *    & Rainer Holzmann (holzmann@cip.ifi.lmu.de)
 */

package clustering.techniques.denistyBased;


import java.io.Serializable;

import weka.clusterers.forOPTICSAndDBScan.DataObjects.DataObject;
import weka.core.Instance;
import weka.core.RevisionHandler;

import org.biojava.bio.motif.Motif;


public class MotifDataObject implements DataObject, Serializable, RevisionHandler {

    /** for serialization */
    private static final long serialVersionUID = -4408119914898291076L;

	static final int UNCLASSIFIED = -1;
    static final int NOISE = Integer.MIN_VALUE;
    static final double UNDEFINED = Integer.MAX_VALUE;
   
    private Motif motif;

    /**
     * Holds the ID of the cluster, to which this DataObject is assigned
     */
    private int clusterID;

    /**
     * Holds the status for this DataObject (true, if it has been processed, else false)
     */
    private boolean processed;

    /**
     * Holds the coreDistance for this DataObject
     */
    private double c_dist;

    /**
     * Holds the reachabilityDistance for this DataObject
     */
    private double r_dist;

    /**
     * Holds the (unique) key that is associated with this DataObject
     */
    private String key;

    // *****************************************************************************************************************
    // constructors
    // *****************************************************************************************************************

    /**
     * Constructs a new DataObject. The original instance is kept as instance-variable
     * @param originalInstance the original instance
     */
    public MotifDataObject(Motif originalMotif, int key) 
    {
        motif = originalMotif;
        this.key = Integer.toString(key);
        clusterID = UNCLASSIFIED;
        processed = false;
        c_dist = UNDEFINED;
        r_dist = UNDEFINED;
    }

    // *****************************************************************************************************************
    // methods
    // *****************************************************************************************************************

    public boolean equals(DataObject dataObject) 
    {
        if (this == dataObject) return true;
      
        
        return false;
    }

    /**
     * Calculates the euclidian-distance between dataObject and this.dataObject
     * @param dataObject The DataObject, that is used for distance-calculation with this.dataObject
     * @return double-value The euclidian-distance between dataObject and this.dataObject
     *                      NaN, if the computation could not be performed
     */
    public double distance(DataObject dataObject) 
    {
       
        return 1 - motif.sim(((MotifDataObject)dataObject).getMotif());
    }

   
   
    /**
     * Sets the clusterID (cluster), to which this DataObject belongs to
     * @param clusterID Number of the Cluster
     */
    public void setClusterLabel(int clusterID) {
        this.clusterID = clusterID;
    }

    /**
     * Returns the clusterID, to which this DataObject belongs to
     * @return clusterID
     */
    public int getClusterLabel() {
        return clusterID;
    }

    /**
     * Marks this dataObject as processed
     * @param processed True, if the DataObject has been already processed, false else
     */
    public void setProcessed(boolean processed) {
        this.processed = processed;
    }

    /**
     * Gives information about the status of a dataObject
     * @return True, if this dataObject has been processed, else false
     */
    public boolean isProcessed() {
        return processed;
    }

    /**
     * Sets a new coreDistance for this dataObject
     * @param c_dist coreDistance
     */
    public void setCoreDistance(double c_dist) {
        this.c_dist = c_dist;
    }

    /**
     * Returns the coreDistance for this dataObject
     * @return coreDistance
     */
    public double getCoreDistance() {
        return c_dist;
    }

    /**
     * Sets a new reachability-distance for this dataObject
     */
    public void setReachabilityDistance(double r_dist) {
        this.r_dist = r_dist;
    }

    /**
     * Returns the reachabilityDistance for this dataObject
     */
    public double getReachabilityDistance() {
        return r_dist;
    }

    
	public Motif getMotif()
	{
		return motif;
	}

	public Instance getInstance()
	{
		// TODO Auto-generated method stub
		return null;
	}

	public String getKey()
	{
		return key;
	}

	public void setKey(String key)
	{
		this.key = key;
		
	}

	public String getRevision()
	{
		// TODO Auto-generated method stub
		return null;
	}

 }