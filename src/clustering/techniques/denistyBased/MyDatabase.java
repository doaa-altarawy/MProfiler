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
import weka.clusterers.forOPTICSAndDBScan.Databases.Database;
import weka.clusterers.forOPTICSAndDBScan.Utils.EpsilonRange_ListElement;
import weka.clusterers.forOPTICSAndDBScan.Utils.PriorityQueue;
import weka.clusterers.forOPTICSAndDBScan.Utils.PriorityQueueElement;
import weka.core.Instances;
import weka.core.RevisionHandler;
import weka.core.RevisionUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;

import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.MotifFinders;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.tools.Sortable;


public class MyDatabase<T>
{
  
    /**
     * Internal, sorted Treemap for storing all the DataObjects
     */
    private Vector<MotifDataObject> data;
    
   
    // *****************************************************************************************************************
    // constructors
    // *****************************************************************************************************************

    /**
     * Constructs a new sequential database and holds the original instances
     * @param instances
     */
    public MyDatabase(Collection<T> instances)
    {
        data = new Vector<MotifDataObject>(instances.size());
        
        int i = 0;
        for (Iterator<T> itr = instances.iterator(); itr.hasNext();)
        {
        	data.add(new MotifDataObject((Motif)itr.next(), i++));
        }
        
    }

    // *****************************************************************************************************************
    // methods
    // *****************************************************************************************************************

    
   
    /**
     * Performs an epsilon range query for this dataObject
     * @param epsilon Specifies the range for the query
     * @param queryDataObject The dataObject that is used as query-object for epsilon range query
     * @return List with all the DataObjects that are within the specified range
     */
    public List epsilonRangeQuery(double epsilon, DataObject queryDataObject) 
    {
        ArrayList epsilonRange_List = new ArrayList();
        Iterator iterator = dataObjectIterator();
        while (iterator.hasNext()) {
            DataObject dataObject = (DataObject) iterator.next();
            double distance = queryDataObject.distance(dataObject);
            if (distance < epsilon) {
                epsilonRange_List.add(dataObject);
            }
        }

        return epsilonRange_List;
    }

    /**
     * Emits the k next-neighbours and performs an epsilon-range-query at the parallel.
     * The returned list contains two elements:
     * At index=0 --> list with all k next-neighbours;
     * At index=1 --> list with all dataObjects within epsilon;
     * @param k number of next neighbours
     * @param epsilon Specifies the range for the query
     * @param dataObject the start object
     * @return list with the k-next neighbours (PriorityQueueElements) and a list
     *         with candidates from the epsilon-range-query (EpsilonRange_ListElements)
     */
    public List k_nextNeighbourQuery(int k, double epsilon, DataObject dataObject) {
        Iterator iterator = dataObjectIterator();

        List return_List = new ArrayList();
        List nextNeighbours_List = new ArrayList();
        List epsilonRange_List = new ArrayList();

        PriorityQueue priorityQueue = new PriorityQueue();
        List temp = new ArrayList<DataObject>();
        
        while (iterator.hasNext()) {
            DataObject next_dataObject = (MotifDataObject) iterator.next();
            double dist = dataObject.distance(next_dataObject);

            if (dist <= epsilon) epsilonRange_List.add(new EpsilonRange_ListElement(dist, next_dataObject));

            if (getNumOfFinders(temp) < k) {
                priorityQueue.add(dist, next_dataObject);
                temp.add(next_dataObject);
            } else {
                if (dist < priorityQueue.getPriority(0)) {
                    priorityQueue.next(); //removes the highest distance
                    priorityQueue.add(dist, next_dataObject);
                    temp.add(next_dataObject);
                }
            }
        }

        while (priorityQueue.hasNext()) {
            nextNeighbours_List.add(0, priorityQueue.next());
        }

        return_List.add(nextNeighbours_List);
        return_List.add(epsilonRange_List);
        return return_List;
    }

    /**
     * Returns an array, element[k] is the k-th distant element from the 
     * diven dataobject, k is number of finders not number of points  
     * @param k
     * @param dataObject
     * @return
     */
    private double[] k_nextNeighbourDist(DataObject dataObject, int k) 
    {
        Iterator iterator = dataObjectIterator();
      
        Sortable<Motif>[] nextNeighbours = new Sortable[data.size()];
        double[] dist = new double[k+1];       
        
        int i=0;
        while (iterator.hasNext()) 
        {
            MotifDataObject next_dataObject = (MotifDataObject) iterator.next();
            double d = dataObject.distance(next_dataObject);
            nextNeighbours[i++] = new Sortable<Motif>(next_dataObject.getMotif(), d, false);
        }
        Arrays.sort(nextNeighbours);
        
        int n=0;
        for (i=0; i<k+1; i++)
        {
        	while (n<nextNeighbours.length && getNumOfFinders(nextNeighbours, n)<i)
        		n++;
        	if (n<nextNeighbours.length)
        		dist[i] = nextNeighbours[n].value;
        	else
        		dist[i] = 2.0;
        }
                
        return dist;
    }
    
    public double[][] getK_Distance(int k)
    {
    	double[][] list = new double[k][data.size()];
    	
    	int i=0;
    	for (Iterator itr=dataObjectIterator(); itr.hasNext();i++)
		{
    		MotifDataObject m = (MotifDataObject)itr.next();
    		
    		double[] temp = k_nextNeighbourDist(m, k);
    		
    		for (int j=0; j<k; j++)
    			list[j][i] = temp[j+1];    		  		
		}
    	
    	for (i=0; i<k; i++)
    		Arrays.sort(list[i]);
    	
    	return list;
    }
    /**
     * Calculates the coreDistance for the specified DataObject.
     * The returned list contains three elements:
     * At index=0 --> list with all k next-neighbours;
     * At index=1 --> list with all dataObjects within epsilon;
     * At index=2 --> coreDistance as Double-value
     * @param minPoints minPoints-many neighbours within epsilon must be found to have a non-undefined coreDistance
     * @param epsilon Specifies the range for the query
     * @param dataObject Calculate coreDistance for this dataObject
     * @return list with the k-next neighbours (PriorityQueueElements) and a list
     *         with candidates from the epsilon-range-query (EpsilonRange_ListElements) and
     *         the double-value for the calculated coreDistance
     */
    public List coreDistance(int minPoints, double epsilon, DataObject dataObject) {
        List list = k_nextNeighbourQuery(minPoints, epsilon, dataObject);

        if (((List) list.get(1)).size() < minPoints)  //(getNumOfFinders(((List) list.get(1))) < minPoints) 
       	{
            list.add(new Double(DataObject.UNDEFINED));
            return list;
        } else {
            List nextNeighbours_List = (List) list.get(0);
            PriorityQueueElement priorityQueueElement =
                    (PriorityQueueElement) nextNeighbours_List.get(nextNeighbours_List.size() - 1);
            if (priorityQueueElement.getPriority() <= epsilon) {
                list.add(new Double(priorityQueueElement.getPriority()));
                return list;
            } else {
                list.add(new Double(DataObject.UNDEFINED));
                return list;
            }
        }
    }

    private int getNumOfFinders(List<MotifDataObject> motifs)
	{		
		Set<MotifFinder> finders = new HashSet<MotifFinder>();
		
		for (Iterator<MotifDataObject> itr=motifs.iterator(); itr.hasNext();)
		{
			Object o = itr.next();
			Motif m = null;
			if (o instanceof MotifDataObject)
				m = ((MotifDataObject)o).getMotif();
			else if (o instanceof EpsilonRange_ListElement)
				m = ((MotifDataObject)((EpsilonRange_ListElement)o).getDataObject()).getMotif();
			finders.add(m.getFinder());			// This is set, elements doesn't repeat				
		}
		
		return finders.size();
	}
    
    private int getNumOfFinders(Sortable<Motif>[] motifs, int n)
	{		
		Set<MotifFinder> finders = new HashSet<MotifFinder>();
		
		for (int i=0; i<n; i++)
			finders.add(motifs[i].obj.getFinder());			// This is set, elements doesn't repeat				
				
		return finders.size();
	}
    /**
     * Returns the size of the database (the number of dataObjects in the database)
     * @return size
     */
    public int size() {
        return data.size();
    }

   
    /**
     * Returns an iterator over all the dataObjects in the database
     * @return iterator
     */
    public Iterator dataObjectIterator() {
        return data.iterator();
    }

         
}
