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

import weka.clusterers.forOPTICSAndDBScan.DataObjects.DataObject;
import weka.clusterers.forOPTICSAndDBScan.OPTICS_GUI.OPTICS_Visualizer;
import weka.clusterers.forOPTICSAndDBScan.OPTICS_GUI.SERObject;
import weka.clusterers.forOPTICSAndDBScan.Utils.EpsilonRange_ListElement;
import weka.clusterers.forOPTICSAndDBScan.Utils.UpdateQueue;
import weka.clusterers.forOPTICSAndDBScan.Utils.UpdateQueueElement;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Utils;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.Iterator;
import java.util.List;


import org.biojava.bio.motif.Motif;


public class OPTICS implements Clusterer
{

    /**
     * Specifies the radius for a range-query
     */
    private double epsilon = 0.7;

    /**
     * Specifies the density (the range-query must contain at least minPoints DataObjects)
     */
    private int minPoints = 2;

    /**
     * Holds the number of clusters generated
     */
    public int numberOfGeneratedClusters;

   
    /**
     * The database that is used for OPTICS
     */
    private MyDatabase<Motif> database;

    /**
     * Holds the time-value (seconds) for the duration of the clustering-process
     */
    private double elapsedTime;

    /**
     * Flag that indicates if the results are written to a file or not
     */
    private boolean writeOPTICSresults = true;

    /**
     * Holds the ClusterOrder (dataObjects with their r_dist and c_dist) for the GUI
     */
    private FastVector resultVector;

    // *****************************************************************************************************************
    // constructors
    // *****************************************************************************************************************

    public OPTICS(int minPoints, double epsilon)
    {
    	this.epsilon = epsilon;
    	this.minPoints = minPoints;    	
    }
    // *****************************************************************************************************************
    // methods
    // *****************************************************************************************************************
   
    /**
     * Generate Clustering via OPTICS
     * @param instances The instances that need to be clustered
     * @throws java.lang.Exception If clustering was not successful
     */
    public void buildClusterer(List<Motif> motifs) throws Exception {
      
        resultVector = new FastVector();
        long time_1 = System.currentTimeMillis();

        numberOfGeneratedClusters = 0;

        database = new MyDatabase<Motif>(motifs);
        UpdateQueue seeds = new UpdateQueue();

        /** OPTICS-Begin */
        Iterator iterator = database.dataObjectIterator();
        while (iterator.hasNext()) {
            DataObject dataObject = (DataObject) iterator.next();
            if (!dataObject.isProcessed()) {
                expandClusterOrder(dataObject, seeds);
            }
        }

        ///////////
        extractDBSCANClustering();
        
        long time_2 = System.currentTimeMillis();
        elapsedTime = (double) (time_2 - time_1) / 1000.0;

        if (writeOPTICSresults) {
            String fileName = "";
            GregorianCalendar gregorianCalendar = new GregorianCalendar();
            String timeStamp = gregorianCalendar.get(Calendar.DAY_OF_MONTH) + "-" +
                    (gregorianCalendar.get(Calendar.MONTH) + 1) +
                    "-" + gregorianCalendar.get(Calendar.YEAR) +
                    "--" + gregorianCalendar.get(Calendar.HOUR_OF_DAY) +
                    "-" + gregorianCalendar.get(Calendar.MINUTE) +
                    "-" + gregorianCalendar.get(Calendar.SECOND);
            fileName = "OPTICS_" + timeStamp + ".TXT";

            FileWriter fileWriter = new FileWriter(fileName);
            BufferedWriter bufferedOPTICSWriter = new BufferedWriter(fileWriter);
            for (int i = 0; i < resultVector.size(); i++) {
                bufferedOPTICSWriter.write(format_dataObject((DataObject) resultVector.elementAt(i)));
            }
            bufferedOPTICSWriter.flush();
            bufferedOPTICSWriter.close();
        }

       new OPTICS_Visualizer(getSERObject(), "OPTICS Visualizer - Main Window");
    }

    /**
     * Expands the ClusterOrder for this dataObject
     * @param dataObject Start-DataObject
     * @param seeds SeedList that stores dataObjects with reachability-distances
     */
    private void expandClusterOrder(DataObject dataObject, UpdateQueue seeds) {
        List list = database.coreDistance(getMinPoints(), getEpsilon(), dataObject);
        List epsilonRange_List = (List) list.get(1);
        dataObject.setReachabilityDistance(DataObject.UNDEFINED);
        dataObject.setCoreDistance(((Double) list.get(2)).doubleValue());
        dataObject.setProcessed(true);

        resultVector.addElement(dataObject);

        if (dataObject.getCoreDistance() != DataObject.UNDEFINED) {
            update(seeds, epsilonRange_List, dataObject);
            while (seeds.hasNext()) {
                UpdateQueueElement updateQueueElement = seeds.next();
                DataObject currentDataObject = (DataObject) updateQueueElement.getObject();
                currentDataObject.setReachabilityDistance(updateQueueElement.getPriority());
                List list_1 = database.coreDistance(getMinPoints(), getEpsilon(), currentDataObject);
                List epsilonRange_List_1 = (List) list_1.get(1);
                currentDataObject.setCoreDistance(((Double) list_1.get(2)).doubleValue());
                currentDataObject.setProcessed(true);

                resultVector.addElement(currentDataObject);

                if (currentDataObject.getCoreDistance() != DataObject.UNDEFINED) {
                    update(seeds, epsilonRange_List_1, currentDataObject);
                }
            }
        }
    }

    /**
     * Returns the internal database
     * 
     * @return the internal database
     */
    private String database_Type = "weka.clusterers.forOPTICSAndDBScan.Databases.SequentialDatabase";
    private String database_distanceType = "clustering.techniques.denistyBased.MotifDataObject";

    public SERObject getSERObject() {
        SERObject serObject = new SERObject(resultVector,
                database.size(),
                1,
                getEpsilon(),
                getMinPoints(),
                writeOPTICSresults,
                database_Type,
                database_distanceType,
                numberOfGeneratedClusters,
                Utils.doubleToString(elapsedTime, 3, 3));
        return serObject;
    }
    /**
     * Wraps the dataObject into a String, that contains the dataObject's key, the dataObject itself,
     * the coreDistance and its reachabilityDistance in a formatted manner.
     * @param dataObject The dataObject that is wrapped into a formatted string.
     * @return String Formatted string
     */
    private String format_dataObject(DataObject dataObject) {
        StringBuffer stringBuffer = new StringBuffer();

        stringBuffer.append("(" + Utils.doubleToString(Double.parseDouble(dataObject.getKey()),
                (Integer.toString(database.size()).length()), 0) + ".) "
                + Utils.padRight(dataObject.toString(), 40) + "  -->  c_dist: " +

                ((dataObject.getCoreDistance() == DataObject.UNDEFINED) ?
                Utils.padRight("UNDEFINED", 12) :
                Utils.padRight(Utils.doubleToString(dataObject.getCoreDistance(), 2, 3), 12)) +

                " r_dist: " +
                ((dataObject.getReachabilityDistance() == DataObject.UNDEFINED) ?
                Utils.padRight("UNDEFINED", 12) :
                Utils.doubleToString(dataObject.getReachabilityDistance(), 2, 3)) + "\n");

        return stringBuffer.toString();
    }

    /**
     * Updates reachability-distances in the Seeds-List
     * @param seeds UpdateQueue that holds DataObjects with their corresponding reachability-distances
     * @param epsilonRange_list List of DataObjects that were found in epsilon-range of centralObject
     * @param centralObject
     */
    private void update(UpdateQueue seeds, List epsilonRange_list, DataObject centralObject) {
        double coreDistance = centralObject.getCoreDistance();
        double new_r_dist = DataObject.UNDEFINED;

        for (int i = 0; i < epsilonRange_list.size(); i++) {
            EpsilonRange_ListElement listElement = (EpsilonRange_ListElement) epsilonRange_list.get(i);
            DataObject neighbourhood_object = listElement.getDataObject();
            if (!neighbourhood_object.isProcessed()) {
                new_r_dist = Math.max(coreDistance, listElement.getDistance());
                seeds.add(new_r_dist, neighbourhood_object, neighbourhood_object.getKey());
            }
        }
    }

    public void extractDBSCANClustering()
    { 	
    	int clusterID = DataObject.NOISE;
    	
    	for (int i=0; i<resultVector.size(); i++)
    	{
    		DataObject dataObject = (DataObject)resultVector.elementAt(i);
    		if (dataObject.getReachabilityDistance() > getEpsilon())
    		{	
			 // UNDEFINED > e
    			if (dataObject.getCoreDistance() <= getEpsilon())
    			{
    				if (clusterID == DataObject.NOISE)
    					clusterID = 0;
    				else clusterID++;
    				dataObject.setClusterLabel(clusterID);
    				numberOfGeneratedClusters++;
    			}
    			else
    				dataObject.setClusterLabel(DataObject.NOISE);
    		}
    		else
    			dataObject.setClusterLabel(clusterID);
    	}
    }
    
    /**
     * Classifies a given instance.
     *
     * @param instance The instance to be assigned to a cluster
     * @return int The number of the assigned cluster as an integer
     * @throws java.lang.Exception If instance could not be clustered
     * successfully
     */
    public int clusterInstance(Instance instance) throws Exception {
        throw new Exception();
    }


    /**
     * Sets a new value for minPoints
     * @param minPoints MinPoints
     */
    public void setMinPoints(int minPoints) {
        this.minPoints = minPoints;
    }

    /**
     * Sets a new value for epsilon
     * @param epsilon Epsilon
     */
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    /**
     * Returns the value of epsilon
     * @return double Epsilon
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * Returns the value of minPoints
     * @return int MinPoints
     */
    public int getMinPoints() {
        return minPoints;
    }

    public int numberOfClusters(){
        return numberOfGeneratedClusters;
    }

    /**
     * Returns the flag for writing actions
     * @return writeOPTICSresults (flag)
     */
    public boolean getWriteOPTICSresults() {
        return writeOPTICSresults;
    }

    /**
     * Sets the flag for writing actions
     * @param writeOPTICSresults Results are written to a file if the flag is set
     */
    public void setWriteOPTICSresults(boolean writeOPTICSresults) {
        this.writeOPTICSresults = writeOPTICSresults;
    }

    /**
     * Returns the resultVector
     * @return resultVector
     */
    public FastVector getResultVector() {
        return resultVector;
    }
          
    /**
     * Returns a description of the clusterer
     * 
     * @return the clusterer as string
     */
    public String toString() {
        StringBuffer stringBuffer = new StringBuffer();
        stringBuffer.append("OPTICS clustering results\n" +
                "============================================================================================\n\n");
        stringBuffer.append("Clustered DataObjects: " + database.size() + "\n");
        stringBuffer.append("Write results to file: " + (writeOPTICSresults ? "yes" : "no") + "\n");
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

