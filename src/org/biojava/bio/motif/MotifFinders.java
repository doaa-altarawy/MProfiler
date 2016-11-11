package org.biojava.bio.motif;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.tools.LogSets;
import org.biojava.bio.tools.MotifTools;
import org.biojava.bio.tools.Sortable;
import clustering.dataTypes.Pattern;
import clustering.distances.MotifDistance;
import clustering.techniques.denistyBased.Clusterer;
import clustering.techniques.denistyBased.MyDBScan;
import clustering.techniques.denistyBased.OPTICS;
import clustering.techniques.fuzzySoftClustering.Fuzzy;
import clustering.techniques.hardClustering.ABF;
import clustering.techniques.hardClustering.HardClustering;
import clustering.techniques.hierarchialClustering.HierarchialClusteringMotif;
import clustering.util.Util;

/**
 * Implementation of MotifVoter: an Ensemble method for motif finding
 * @author Doaa Altarawy
 *
 */
public class MotifFinders
{
	public static enum  Method {DBScan, OPTIC, KMeans, Hierarcial, MotifVoter, Clique, MotiffProfiler};
	
	// Input
	Dataset dataset;		
	Vector<Motif> inputMotifs;				// vector of all given motifs
	
	Map<Motif, Map<Motif, Double>> sim;
	Map<Motif, Vector<Motif>> orderedSim;
	
	// Results
	Vector<Motif> accepted;
	List<Motif> notAccepted;	
	Motif finalMotif;
	
	Map<Integer, Vector<Motif>> clusters ;
	int MAX_NUM_OF_FINDERS;
	int minpts;
	double epsilon;
	int numOfClusters;
		
	Method type = Method.MotifVoter; 
	
	/* Constructors */	
	public MotifFinders(Dataset dataset, Vector<Motif> motifs)
	{
		this.dataset = dataset;
		this.inputMotifs = motifs;
		MAX_NUM_OF_FINDERS = MotifTools.getNumOfFinders(inputMotifs);
	}
	
	public MotifFinders(Dataset dataset, Vector<Motif> motifs, int m, double e, Method type, int numOfClusters)
	{
		this.dataset = dataset;
		this.inputMotifs = motifs;		
		minpts = m;
		epsilon = e;
		this.type = type;
		this.numOfClusters = numOfClusters;
		MAX_NUM_OF_FINDERS = MotifTools.getNumOfFinders(inputMotifs);
	}	
	/*----------------------------------------------------------------------------*/
	
	public Motif findMotif() throws Exception
	{
		 long time_1 = System.currentTimeMillis();
		 computeSortedSimilarity(inputMotifs);
		 MAX_NUM_OF_FINDERS = MotifTools.getNumOfFinders(inputMotifs);		 
		 int size = 30;//inputMotifs.size()/5;
		 double sim = 0.6;
		 
		// step 2
		if (type==Method.MotifVoter)
		{						
			getAcceptedSetMV(inputMotifs);		
		}
		else if (type==Method.KMeans)
		{
			Vector<Motif> centers = getInitCenters(inputMotifs, size, sim);
			computeSortedSimilarity(centers);
			getClusteringKMeans(centers, numOfClusters);			 
		}
		else if (type==Method.Hierarcial)
		{
			Vector<Motif> centers = getInitCenters(inputMotifs, size, sim);
			computeSortedSimilarity(centers);
			getHierarchialClustering(centers, numOfClusters);  
		}
		else if (type==Method.Clique)
		{
			Vector<Motif> centers = getInitCenters(inputMotifs, size, sim);
			computeSortedSimilarity(centers);
			getAcceptedSetClique(centers);
		}
		else if (type==Method.DBScan || type==Method.OPTIC)
		{
			Vector<Motif> centers = getInitCenters(inputMotifs, size, sim);
			computeSortedSimilarity(centers);
			getClusteresDBScan(inputMotifs, type); 
		}
		else if (type==Method.MotiffProfiler)
		{			
			Vector<Motif> centers = getInitCenters(inputMotifs, size, sim);
			computeSortedSimilarity(centers);
			getAcceptedSetMProfiler(centers);		
		}
		
		long time_2 = System.currentTimeMillis();
        System.out.println("ElapsedTime ="+ ((double) (time_2 - time_1) / 1000.0));
		      
		return finalMotif;
	}	
	
	/*-----------------------------------------------------------------------------*/
	
	private List<Motif> getClusteringFuzzy(int numOfClusters) throws Exception
	{
		if (numOfClusters==-1)
			numOfClusters = MotifTools.getNumOfFinders(inputMotifs);
		
		Vector<Pattern> data = new Vector<Pattern>(inputMotifs.size());
		
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			Vector<Motif> v = new Vector<Motif>(1);
			v.add(m);
			Pattern pattern = new Pattern(m.getFinder().name(), v);
			data.add(pattern);
		}
		
		Fuzzy fuzzy = new Fuzzy();
		fuzzy.partition(data, new MotifDistance(), 2, numOfClusters, 1);
		
		
		clusters = new HashMap<Integer, Vector<Motif>>(numOfClusters+1);
		for (int n=0; n<numOfClusters; n++)
		{
			Vector<Motif> list = new Vector<Motif>();
			int i=0;
			for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext(); i++)
			{
				Motif m = itr.next();
				if (fuzzy.result.membership[n][i]>0.2)
				{
					list.add(m);
					System.out.println("cluster #"+n+" included: "+fuzzy.result.membership[n][i]);					

				}
				else 
					System.out.println("cluster #"+n+" wont include: "+fuzzy.result.membership[n][i]);					
			}
			Motif m = MotifTools.extractSitesMultiFinder(list);
			list.add(0, m);
			clusters.put(n, list);
		}		
		
		///
		computeSortedSimilarity(inputMotifs);
		clusters.put(numOfClusters, getBestCluster(clusters));
		accepted = new Vector<Motif>(clusters.get(clusters.size()-1));
		accepted.remove(0);
		finalMotif = MotifTools.extractSitesMultiFinder(accepted);
		return null;
	}
	
	/*-----------------------------------------------------------------------------*/
	
	
	private void getClusteresDBScan(Vector<Motif> inputMotifs, Method type) throws Exception
	{
				
		Clusterer dbScan = null;
		if (type == Method.DBScan)
			dbScan = new MyDBScan(minpts, epsilon);
		else if (type== Method.OPTIC)
			dbScan = new OPTICS(minpts, epsilon);
		
		dbScan.buildClusterer(inputMotifs);
		System.out.println(dbScan.toString());
		
		clusters = new HashMap<Integer, Vector<Motif>>(dbScan.numberOfClusters()+1);
		for (int i=0; i<dbScan.numberOfClusters(); i++)
		{
			List<Motif> cluster = dbScan.getCluster(i);
			Motif finalM = MotifTools.extractSitesMultiFinder(cluster);
			Vector<Motif> temp = new Vector<Motif>();
			temp.add(finalM);
			temp.addAll(cluster);
			clusters.put(i, temp);
			System.out.println("Clutser "+i+" has "+cluster.size()+" Motifs.");
		}		
		
		///
		computeSortedSimilarity(inputMotifs);
		Vector<Motif> bestCluster = getBestCluster(clusters);
		
		if (bestCluster!=null)
		{
			clusters.put(dbScan.numberOfClusters(), bestCluster);
			accepted = new Vector<Motif>(clusters.get(clusters.size()-1));
			accepted.remove(0);
			finalMotif = MotifTools.extractSitesMultiFinder(accepted);
		}
		
	}
	
	/*-----------------------------------------------------------------------------*/
	private List<Motif> getClusteringKMeans(Vector<Motif> inputMotifs, int numOfClusters) throws Exception
	{
		if (numOfClusters==-1)
			numOfClusters = MotifTools.getNumOfFinders(inputMotifs);
		
		Vector<Pattern> data = new Vector<Pattern>(inputMotifs.size());
		
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			Vector<Motif> v = new Vector<Motif>(1);
			v.add(m);
			Pattern pattern = new Pattern(m.getFinder().name(), v);
			data.add(pattern);
		}
		
		HardClustering kmeans = new ABF();//KMeans();
		kmeans.partition(data, new MotifDistance(), numOfClusters);
		
		
		clusters = new HashMap<Integer, Vector<Motif>>(numOfClusters+1);
		for (int n=0; n<numOfClusters; n++)
		{
			Vector<Motif> list = kmeans.getCluster(n);
			//list = getMaxConfidenceScore(list);
			Motif m = MotifTools.extractSitesOneFinder(list);
			list.add(0, m);
			clusters.put(n, list);
		}		
		
		///
		computeSortedSimilarity(inputMotifs);
		Vector<Motif> bestCluster = getBestCluster(clusters);
		//bestCluster = getBestClusterDistSqr(clusters);
		
		if (bestCluster!=null)
		{
			clusters.put(numOfClusters, bestCluster);
			accepted = new Vector<Motif>(clusters.get(clusters.size()-1));
			finalMotif = accepted.get(0);
			accepted.remove(0);
		}
		return null;
	}
	/*-----------------------------------------------------------------------------*/
	
	private List<Motif> getHierarchialClustering(Vector<Motif> inputMotifs, int numOfClusters) throws Exception
	{
		LogSets logSet = new LogSets(dataset, 0.0, 0.0, "Hirar"); 
		
		if (numOfClusters==-1)
			numOfClusters = 1;
		
		Vector<Pattern> data = new Vector<Pattern>(inputMotifs.size());
		
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			Vector<Motif> v = new Vector<Motif>(1);
			v.add(m);
			Pattern pattern = new Pattern(m.getFinder().name(), v);
			data.add(pattern);
		}
		
		HierarchialClusteringMotif clustering = new HierarchialClusteringMotif();	////////
		clustering.partition(data, new MotifDistance(), data.size(), numOfClusters);
		
		
		clusters = new HashMap<Integer, Vector<Motif>>(numOfClusters+1);
		int n = 0;
		
		//for (int i=1; i<clustering.result.numOfIterations; i++)
		{
			//System.out.println("Iteration num:"+i);
			for (Iterator<Vector<Motif>> itr=clustering.result.getClusters().iterator(); itr.hasNext();)
			{
				Vector<Motif> list = itr.next();
				//if (list.size()<3) continue;
				Vector<Motif> temp = new Vector<Motif>(inputMotifs);
				temp.removeAll(list);
				//clusterCentersDistanceIndex(list, temp, logSet, false, i); 
				Motif m = MotifTools.extractSitesOneFinder(list);
				list.add(0, m);
				clusters.put(n++, list);
			}		
		}
		///
		computeSortedSimilarity(inputMotifs);
		Vector<Motif> bestCluster = getBestClusterSilhouette(clusters);
		
		if (bestCluster!=null)
		{
			clusters.put(numOfClusters, bestCluster);
			accepted = new Vector<Motif>(bestCluster);//new Vector<Motif>(clusters.get(clusters.size()-1));
			finalMotif = accepted.get(0);
			accepted.remove(0);
		}
	
		logSet.writeToFile();
		//numOfClusters = 0;
		return null;
	}
	/*-----------------------------------------------------------------------------*/
	/*-----------------------------------------------------------------------------*/
	
	StringBuffer log = new StringBuffer();
	private Vector<Motif> getInitCenters(Vector<Motif> motifs, int size, double maxSim) throws Exception
	{
		List<Motif> x;
		int numOfFinders = 0;		
		Vector<Motif> centers = new Vector<Motif>(motifs.size());
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{					
			x = new ArrayList<Motif>(motifs.size());			
			Motif z = itr.next(); 
			x.add(z);
			
			System.out.println("Finder z:"+z.getFinder().name()+"\n----------");
				
			for (Iterator<Motif> itr2 = orderedSim.get(z).iterator(); itr2.hasNext();)	
			{
				
				Motif p = itr2.next();
				if (z.equals(p)) //skip same motif
					continue;
				
				x.add(p);
			
				System.out.print("Sim="+sim.get(z).get(p));
				numOfFinders = MotifTools.getNumOfFinders(x);
				if (numOfFinders < 3)
				{
					System.out.print("Finders still <3..\n");
					continue;
				}
				
				// Heuristic
			//	if (x.size()< ((MAX_NUM_OF_FINDERS)))
			//		continue;
				if (x.size()> size)
				{
					System.out.println("Max Size exceeded..\n");
					break;
				}
												
				
				Motif c = MotifTools.extractSitesMultiFinder(x);
				c.setFinder(z.finder);
				if (c.getAllSites().size()<2) // make it <2 if using PWM
				{
					System.out.println("sites size still<2");
					continue;
				}
				double val = centers.isEmpty()?-1:c.sim(centers.lastElement());
				
				if (centers.isEmpty() || val<maxSim)
				{
					System.out.println("Added.., sim="+val);
					centers.add(c);							
				}
				else
					System.out.println("Not Added.., sim="+val);
			}
			
		}
		
		//FileConverterTools.writeMotifVoter(centers, new File(".\\"+dataset.name()+"_centers.voter"));
		System.out.println("******\nCenters size="+centers.size());
		
		return centers;
	}
	/*-----------------------------------------------------------------------------*/
	private List<Motif> getAcceptedSetMProfiler(Vector<Motif> motifs) throws Exception
	{
		Vector<Motif> maxSet = null;		
		double maxWeight = -1.0;		
		List<Motif> x, p_x;		
		double weight = -1.0;		
		LogSets logSet = new LogSets(dataset, 0.25, 0.2, "MProfiler"); 
		
		int n = 0;
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{					
			x = new ArrayList<Motif>(motifs.size());
			p_x = new ArrayList<Motif>(motifs);
			Motif z = itr.next(); 
			x.add(z);
			p_x.remove(z);
			System.out.println("n="+ ++n);
				
			for (Iterator<Motif> itr2 = orderedSim.get(z).iterator(); itr2.hasNext();)	
			{
				
				Motif p = itr2.next();
				if (z.equals(p)) //skip same motif
					continue;
				
				x.add(p);
				p_x.remove(p);				
				
								
				// Heuristic
				if (x.size()< MAX_NUM_OF_FINDERS)
					continue;
				if (x.size()> 20)//((MAX_NUM_OF_FINDERS)*8))
					break;
								
				weight = clusterCentersDistanceIndex(x, p_x, null, false, n, true); 
			
				if (weight>maxWeight)
				{
					System.out.println("Max changed to:"+weight);						
					maxWeight = weight;					
					maxSet = new Vector<Motif>(x);
				}				
			}
			
		}
				
		accepted = maxSet; 	
		finalMotif = MotifTools.extractSitesOneFinder(accepted);
		
		// write log file
		logSet.writeToFile();
				
		return accepted;
	}
	/*-----------------------------------------------------------------------------*/
	private List<Motif> getAcceptedSetMV(Vector<Motif> motifs) throws Exception
	{
		Vector<Motif> maxSet = new Vector<Motif>();		
		double maxWeight = -1.0;		
		List<Motif> x, p_x;		
		double weight = -1.0;		
		LogSets logSet = new LogSets(dataset, 0.25, 0.2, "MV_confi_all"); 
		
		int n = 0;
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{					
			x = new ArrayList<Motif>(motifs.size());
			p_x = new ArrayList<Motif>(motifs);
			Motif z = itr.next(); 
			x.add(z);
			p_x.remove(z);
			System.out.println("n="+ ++n);
				
			//MAX_NUM_OF_FINDERS = MotifTools.getNumOfFinders(orderedSim.get(z));
			for (Iterator<Motif> itr2 = orderedSim.get(z).iterator(); itr2.hasNext();)	
			{
				
				Motif p = itr2.next();
				if (z.equals(p)) //skip same motif
					continue;
				
				x.add(p);
				p_x.remove(p);				
				
				
				int numOfFinders = MotifTools.getNumOfFinders(x);
				if ( numOfFinders  < MAX_NUM_OF_FINDERS)
					continue;
				
				// Heuristic
				//if (x.size()< MAX_NUM_OF_FINDERS)
					//continue;
			//	if (x.size()> (MAX_NUM_OF_FINDERS*8))
				//	break;
								
				weight = clusterCentersDistanceIndex(x, p_x, null, false, n, false); 
			
				if (weight>maxWeight)
				{
					System.out.println("Max changed to:"+weight);						
					maxWeight = weight;					
					maxSet = new Vector<Motif>(x);
				}				
			}
			
		}
				
		accepted = maxSet; 	
		finalMotif = accepted.size()>0?MotifTools.extractSitesMV(accepted):new Motif(MotifFinder.MotifVoter);
		
		// write log file
		logSet.writeToFile();
				
		return accepted;
	}
	/*-----------------------------------------------------------------------------*/
	private List<Motif> getAcceptedSetClique(Vector<Motif> inputMotifs) throws Exception
	{
		Vector<Motif> maxSet = null;
		int maxNumOfFinders = 0;
		double maxWeight = -1.0;		
		List<Motif> x, p_x;
		int numOfFinders = 0;
		double weight = -1.0;
		LogSets logSet = new LogSets(dataset, 0.25, 0.25, "clique_Profile_old");
		
		System.out.println("MAX_NUM_OF_FINDERS="+MAX_NUM_OF_FINDERS);
		
		try{
		/*
		Sortable<List<Motif>>[] initialClique = new Sortable[inputMotifs.size()];
		Map<List<Motif>, Double> bestInitialClique = new HashMap<List<Motif>, Double>(inputMotifs.size());
		
		int i=0;
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			List<Motif> clique = new ArrayList<Motif>();
			clique = getInitialClique(itr.next());
			List<Motif> rest = new ArrayList<Motif>(inputMotifs);
			rest.removeAll(clique);
			double w = clusterCentersDistanceIndex(clique, rest, logSet, false, 0);
			numOfFinders = MotifTools.getNumOfFinders(clique);
			//if (maxNumOfFinders<numOfFinders)
			//	maxNumOfFinders = numOfFinders;
			initialClique[i++] = new Sortable<List<Motif>>(clique, w);
		}
		//numOfFinders = numOfFinders / inputMotifs.size();
		Arrays.sort(initialClique);
		//int count = 0;
		for (i=0; i<initialClique.length; i++)
		{			
			//if (getNumOfFinders(initialClique[i].obj)>= 4)
				bestInitialClique.put(initialClique[i].obj, initialClique[i].value);
			//if (++count > initialClique.length/2)	// take only half top weight sets
			//	break;
		}
		*/
		int i = 0;
		//for (Iterator<List<Motif>> itr=bestInitialClique.keySet().iterator(); itr.hasNext();)
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
		{
			//List<Motif> z = itr.next();
			Motif z = itr.next();
			
			System.out.println("i="+i++);
			
			//List returnList = getMaxClique(new ArrayList<Motif>(z), logSet, bestInitialClique.get(z));
			List returnList = getMaxClique(inputMotifs, z, logSet, 0);
			x = (List<Motif>)returnList.get(0);
			
			//p_x = (List<Motif>)returnList.get(1);
		
			
									
			weight = (Double)returnList.get(2);//clusterCentersDistanceIndex(x, p_x, logSet, false, bestInitialClique.get(z)); 
			
			if (weight>maxWeight)// && numOfFinders>=maxNumOfFinders)
			{				
				maxWeight = weight;
				//maxNumOfFinders = numOfFinders;
				maxSet = new Vector<Motif>(x);			
			}
		}
		
		accepted = maxSet; //getMaxConfidenceScore(maxSet);		
		finalMotif = MotifTools.extractSitesOneFinder(accepted);
		
		}
		catch (Exception e) {e.printStackTrace();}
		finally{
			Motif center = MotifTools.extractSitesOneFinder(inputMotifs);
			
			double xmean = meanDist(inputMotifs, center);		
			double xvar = varianceDist(inputMotifs, center, xmean);			
			logSet.logSet(center, inputMotifs, xmean, xvar, 0, 0, 1, 0, false, -1, -1); // log the whole set			
			logSet.writeToFile();
		}
		
		//FileOutputStream logFile = new FileOutputStream("log2.txt");
		//logFile.write(log.toString().getBytes());
		//logFile.flush();
		//logFile.close();
		
		return accepted;
	}	
	
	/*-----------------------------------------------------------------------------*/
	//private List getMaxClique(List<Motif> clique, LogSets logSet, double source) throws Exception
	private List getMaxClique(Vector<Motif> inputMotifs, Motif start, LogSets logSet, double source) throws Exception
	{			
		Map<Motif, Integer> motifPointer = new HashMap<Motif, Integer>(inputMotifs.size());
		List<Motif> clique = new ArrayList<Motif>();
		clique.add(start);
		List<Motif> maxSet = new ArrayList<Motif>(clique);
		List<Motif> p_x = new ArrayList<Motif>(inputMotifs);
		p_x.removeAll(clique);	
	//	Motif center = clique.get(0);//getPerwiseIntersection(clique);
		double maxWeight = -1;
		double weight;// = silhouetteWidth(start, clique, p_x);
						
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
			motifPointer.put(itr.next(), 0);
		
		
	//	double separatness[] = new double[2];
	//	double meanDist[] = new double[2];
		//int current = 0; // next will be: current = (current+1)%2;
		
		
		while (true)
		{
			if (clique.size()> 25 ) // heuristic
				break;
			
			Motif pointingMotif = null;
			Motif maxMotif = null;
			double maxSim = 0;
			
			for (Iterator<Motif> itr=clique.iterator(); itr.hasNext();)
			{
				Motif m = itr.next();
				Vector<Motif> orderedMotifs = orderedSim.get(m);
				Integer ptr = motifPointer.get(m);
				if (ptr >= orderedMotifs.size()) // if similar motifs list to m finished 
					continue;
				Motif p = orderedMotifs.get(ptr);	// next similar motif to m
				while (clique.contains(p))	// check if p already taken, find another
				{
					motifPointer.put(m, ++ptr);	// get next pointer
					if (ptr < orderedMotifs.size())
						p = orderedMotifs.get(ptr);	// next similar motif to m
					else
					{p=null; break;}					
				}
				if (p==null) continue;
				
				double d = sim.get(m).get(p); //center.sim(p);
				if ( d > maxSim)
				{
					maxSim = d;
					maxMotif = p;
					pointingMotif = m;
				}
			}
			
			if (Util.equalZero(maxSim)) // all sets finished and maxSize wasn't assigned to any point
				break;
						
			
			clique.add(maxMotif);
			p_x.remove(maxMotif);
			motifPointer.put(pointingMotif, motifPointer.get(pointingMotif)+1);	// get next pointer
			//center = getPerwiseIntersection(clique); //getCenter(clique);
					
									
			if (clique.size()< MAX_NUM_OF_FINDERS)
				continue;
			
			//Motif x_center = extractSites(clique);
			//Motif p_center = extractSites(p_x);
			//separatness[current] = separatness(x_center, p_center);
			//meanDist[current] = meanDist(clique, x_center);
			
			//int numOfFinders = MotifTools.getNumOfFinders(clique);
			//if ( numOfFinders  < 4)
			//	continue;
			
			weight = clusterCentersDistanceIndex(clique, p_x, logSet, false, source, true);
			
			/*
			int prev = (current+1)%2;
			
			if ((separatness[prev]-separatness[current] > 0.01)
					&& meanDist[prev]-meanDist[current] > 0.02)
			{
				System.out.println("Set Choosen..");
				System.out.println("Finders:"+numOfFinders+", size="+clique.size());	
				clique.remove(maxMotif);
				return clique;
			}
			current = prev;
			*/
			
			if (weight>maxWeight)
			{
				//System.out.println("LocalMax changed to:"+weight);
				//System.out.println("Size="+clique.size());	
				maxWeight = weight;
				maxSet = new Vector<Motif>(clique);			
			}		
			
		}
		
		List returnLists = new ArrayList(3);
		returnLists.add(maxSet);
		returnLists.add(p_x);
		returnLists.add(maxWeight);
		return returnLists;
	}
	
	/*-----------------------------------------------------------------------------*/
	private List<Motif> getInitialClique(Motif start) throws Exception
	{		
		List<Motif> clique = new ArrayList<Motif>();
		Map<Motif, Integer> motifPointer = new HashMap<Motif, Integer>(inputMotifs.size());	
		Map<MotifFinder, Integer> finders = MotifTools.getEmptyFindersMap();
		int numOfFinders = MAX_NUM_OF_FINDERS; /////////////
		
		for (Iterator<Motif> itr=inputMotifs.iterator(); itr.hasNext();)
			motifPointer.put(itr.next(), 0);
		
		clique.add(start);		
		finders.put(start.getFinder(), finders.get(start.getFinder())+1);
		Motif center = start;
		
		while (numOfFinders < MAX_NUM_OF_FINDERS)
		{			
			Motif pointingMotif = null;
			Motif maxMotif = null;
			double maxSim = 0;
			
			for (Iterator<Motif> itr=clique.iterator(); itr.hasNext();)
			{
				Motif m = itr.next();
				Vector<Motif> orderedMotifs = orderedSim.get(m);
				Integer ptr = motifPointer.get(m);
				if (ptr >= orderedMotifs.size()) // if similar motifs list to m finished 
					continue;
				Motif p = orderedMotifs.get(ptr);	// next similar motif to m
				while (clique.contains(p) || finders.get(p.getFinder())>=1)	// check if p already taken, find another
				{
					motifPointer.put(m, ++ptr);	// get next pointer
					if (ptr < orderedMotifs.size())
						p = orderedMotifs.get(ptr);	// next similar motif to m
					else
						{p=null; break;}					
				}
				if (p==null) continue;
				
				double d = center.sim(p); //sim.get(m).get(p);
				if ( d > maxSim)
				{
					maxSim = d;
					maxMotif = p;
					pointingMotif = m;
				}
			}
			
			if (Util.equalZero(maxSim)) // all sets finished and maxSize wasn't assigned to any point
				break;
			
			clique.add(maxMotif);
			finders.put(maxMotif.getFinder(), finders.get(maxMotif.getFinder())+1);
			motifPointer.put(pointingMotif, motifPointer.get(pointingMotif)+1);	// get next pointer
			center = MotifTools.getCenter(clique);//getPerwiseIntersection(clique);
										
			numOfFinders = MotifTools.getNumOfFinders(clique);					
		}
		
		//System.out.println("Final # of Finders:"+numOfFinders);
		return clique;
	}
	
	 
	/*-----------------------------------------------------------------------------*/
	/*--------------------------- Validation Functions ----------------------------*/
	public double weightRatio(List<Motif> x, List<Motif> p)
	{
		double a = weight(x);
		double b = weight(p);
		double weight;
		
		if (b==0)
			weight = Double.NaN;
		else				
			weight = a / b;
		
		return weight;
	}
	
	public double weight(List<Motif> motifs)
	{
		return weight(motifs, true);
	}
	
	public double weight(List<Motif> motifs, boolean addsim)
	{
		return weight(motifs, sim(motifs), addsim);
	}
	
	public double weight(List<Motif> motifs, double simm, boolean addSim)
	{
		double var = var(motifs, simm);		
		double total = 0.0;
								
		if (var != 0)
		{
			if (addSim)
				total = ( simm / Math.sqrt(var));
			else
				total = ( 1.0 / Math.sqrt(var));
		}
		
		return total;
	}
	
	public double var(List<Motif> motifs, double simm)
	{
		double sum = 0.0;		
				
		Object[] list = motifs.toArray();
		//int count = 0;
		
		for (int i=0; i<list.length; i++)
		{			
			//Map<Motif, Double> mSim = sim.get(list[i]);
			for (int j=i+1; j<list.length; j++)
			{				
				sum += Math.pow(sim.get(list[i]).get(list[j]) - simm, 2.0);				
			}
		}
		
		sum = 2.0 * sum + (double)motifs.size() * Math.pow((1.0-simm), 2.0); // symmetric matrix ^ 2, plus diagonal (1-simm)^2
		
		//sum = sum / (double)count;
				
		return sum;
	}
	
	/*
	 * public double var(List<Motif> motifs, double simm)
	{
		double sum = 0.0;		
				
		List<Motif> list = new ArrayList<Motif>(motifs);
		//int count = 0;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			list.remove(m);
			Map<Motif, Double> mSim = sim.get(m);
			for (Iterator<Motif> itr2=list.iterator(); itr2.hasNext();)
			{
				Motif temp = itr2.next();
				sum += Math.pow(mSim.get(temp) - simm, 2.0);
				//count ++;
			}
		}
		
		sum = 2.0 * sum + (double)motifs.size() * Math.pow((1.0-simm), 2.0); // symmetric matrix ^ 2, plus diagonal (1-simm)^2
		
		//sum = sum / (double)count;
				
		return sum;
	}
	
	 */
	
	/*-----------------------------------------------------------------------------*/
	public double clusterCentersDistanceIndex(List<Motif> x, List<Motif> p, LogSets logSet, boolean saveSet, double i, boolean singleFinder) throws Exception
	{
		double weight = 0;
		
		Motif x_Center = singleFinder?MotifTools.extractSitesOneFinder(x):MotifTools.extractSitesMultiFinder(x);
		
		//Motif p_Center = singleFinder?MotifTools.getPerwiseIntersection(p):MotifTools.getCenter(p);
		
		//double separatness = separatness(x_Center, p_Center);
		/*
		//weight = weight(x, false) * separatness / (x_compact1>x_compact2?x_compact1:x_compact2);
		//double d1 = maxDist(x, x_Center);
		double xmean = meanDist(x, x_Center);
		double pmean = meanDist(p, p_Center);
		//double d3 = avgSimDifferenceWithCenter(x);
		double xvar = varianceDist(x, x_Center, xmean);
		double pvar = varianceDist(p, p_Center, pmean);
		if (Util.equalZero(xvar) || Util.equalZero(pmean) || xvar<0.0009)	// those cases are noise, ignore them
			weight = 0;	
		else weight =  separatness * xmean / pmean * pvar/ xvar;// / (Util.equalZero(var)?1:var);
		
		if (logSet!=null)
			logSet.logSet(x, p, xmean, xvar, pmean, pvar, separatness, weight, saveSet, i);
		*/
		PositionWeightMatrix pwm = new PositionWeightMatrix(x_Center, i+dataset.name()+"_"+type.name());
		
		if (singleFinder && pwm.columns()>50)
			return 0;
		//double map = pwm.getMAPScore();
		
		double xmean = sim(x);
		double pmean = sim(p);		
		double xvar = Math.sqrt(var(x, xmean));
		double pvar = Math.sqrt(var(p, pmean));
		if (Util.equalZero(xvar) || Util.equalZero(pmean))	// those cases are noise, ignore them
			weight = 0;	
		//else weight =  xmean/ xvar / (-map);//pmean * pvar;// / (Util.equalZero(var)?1:var);
		else weight =  xmean/ xvar / pmean * pvar;// / (Util.equalZero(var)?1:var);
		
		
		if (logSet!=null)
			//logSet.logSet(x_Center, x, xmean, xvar, pwm.numOfWords, pwm.columns(), 0, weight, saveSet, pwm.totalInfoContent, -map);
			logSet.logSet(x_Center, x, xmean, xvar, pmean, pvar, 0, weight, saveSet, 0, 0);
				
		return weight;
	}
	 
	/*-----------------------------------------------------------------------------*/
	
	public double silhouetteWidth(Motif motif, List<Motif> x, List<Motif> p)
	{
		double width = 0;
		double x_sim = 0;
		double p_sim = 0;
		Map<Motif, Double> simm = sim.get(motif);
		
		for (Iterator<Motif> itr=x.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			Double d = simm.get(m);
			if (d!=null)
				x_sim =+ 1.0 - d;
		}
	
		for (Iterator<Motif> itr=p.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			Double d = simm.get(m);
			if (d!=null)
				p_sim =+ 1.0 - d;
		}
		
		x_sim = x_sim / (double)x.size();
		p_sim = p_sim / (double)p.size();
		
	
		width = (p_sim - x_sim) / (x_sim>p_sim?x_sim:p_sim);
		
		return width;
	}
	
	/*-----------------------------------------------------------------------------*/
	public double separatness(Motif c1, Motif c2)
	{
		return 1.0 - c1.sim(c2);
	}
	/*-----------------------------------------------------------------------------*/

	public double compactness_MaxDist(List<Motif> x) throws Exception
	{
		return maxDist(x, MotifTools.getCenter(x));
	}
	
	public double meanDist(List<Motif> x) throws Exception
	{
		return meanDist(x, MotifTools.getCenter(x));
	}
	
	public double maxDist(List<Motif> x, Motif center) throws Exception
	{
		double compact = 0;
		
		for (Iterator<Motif> itr=x.iterator(); itr.hasNext();)
		{
			double d = 1.0 - center.sim(itr.next());
			if (d>compact)
				compact = d;
		}
		
		return compact;
	}
	
	public double minDist(List<Motif> x, Motif center) throws Exception
	{
		double compact = Double.MAX_VALUE;
		
		for (Iterator<Motif> itr=x.iterator(); itr.hasNext();)
		{
			double d = 1.0 - center.sim(itr.next());
			if (d<compact)
				compact = d;
		}
		
		return compact;
	}
	
	public double meanDist(List<Motif> x, Motif x_Center) throws Exception
	{
		double mean = 0.0;
		
		for (Iterator<Motif> itr=x.iterator(); itr.hasNext();)
		{
			mean += x_Center.sim(itr.next());
			
		}
		mean = mean / (double)x.size();
		
		return mean;
	}
	
	/*-----------------------------------------------------------------------------*/
	
	public double varianceDist(List<Motif> x, Motif x_Center, double meanDist)
	{		
		double var = 0;
		
		for (Iterator<Motif> itr=x.iterator(); itr.hasNext();)
		{
			double sim = x_Center.sim(itr.next());
			var += Math.pow(sim - meanDist, 2.0);
			
		}
		var = var / (double) x.size();
		//var = Math.sqrt(var);
		
		return var;
	}
	
	public double avgSimDifferenceWithCenter(List<Motif> motifs)
	{
		double sum = 0;		
		double total = 0;
		double simm = sim(motifs);
		//System.out.println("sim(X)="+simm+ ", X.size="+motifs.size());
		List<Motif> list = new ArrayList<Motif>(motifs);
		int count = 0;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			list.remove(m);
			Map<Motif, Double> mSim = sim.get(m);
			for (Iterator<Motif> itr2=list.iterator(); itr2.hasNext();)
			{
				Motif temp = itr2.next();
				sum += Math.abs(mSim.get(temp) - simm);		
				count++;
			}
		}
		
		total = sum / (double)count;
		
		
		return total;
	}
	/*-----------------------------------------------------------------------------*/
	/*-----------------------------------------------------------------------------*/
	public void computeSortedSimilarity(Vector<Motif> inputMotifs) throws IOException
	{
		Sortable<Motif>[][] similarities = new Sortable[inputMotifs.size()][inputMotifs.size()];
		sim = new HashMap<Motif, Map<Motif,Double>>(inputMotifs.size());
		orderedSim = new HashMap<Motif, Vector<Motif>>(inputMotifs.size());
		
		for (int i=0; i<inputMotifs.size(); i++)
		{			
			for (int j=0; j<inputMotifs.size(); j++)
			{				
				double d = inputMotifs.get(i).sim(inputMotifs.get(j));
				similarities[i][j] = new Sortable<Motif>(inputMotifs.get(j), d);
			}
			Arrays.sort(similarities[i]);	
			Map<Motif, Double> map = new HashMap<Motif, Double>(inputMotifs.size()-1);
			Vector<Motif> v = new Vector<Motif>(inputMotifs.size()-1);
			for (int k=0; k<similarities[i].length; k++)
			{
				map.put(similarities[i][k].obj, similarities[i][k].value);
				if (similarities[i][k].value>0)	// ok, by E. Wijaya
					v.add(similarities[i][k].obj);
			}
			sim.put(inputMotifs.get(i), map);	
			orderedSim.put(inputMotifs.get(i), v);
		}			
	
	}
		
	
	public Vector<Motif> getBestCluster(Map<Integer, Vector<Motif>> clusters)
	{
				
		int bestClusterIndex = -1;
		double bestWeight = Integer.MAX_VALUE;
		
		for (Iterator<Integer> clusterItr= clusters.keySet().iterator(); clusterItr.hasNext();)
		{
			int currentClusterIndex = clusterItr.next();
			List<Motif> currentCluster = new ArrayList<Motif>(clusters.get(currentClusterIndex));
			currentCluster.remove(0);	// its center
			List<Motif> restOfClusters = new ArrayList<Motif>();
			for (Iterator<Integer> allClusterItr= clusters.keySet().iterator(); allClusterItr.hasNext();)
			{
				int i = allClusterItr.next();
				if (i!= currentClusterIndex)
				{
					List<Motif> temp = new ArrayList<Motif>(clusters.get(i));
					temp.remove(0);	// its center
					restOfClusters.addAll(temp);
				}
			}
			
			double a = weight(currentCluster);
			double b = weight(restOfClusters);
			double weight = a / b;
			
			System.out.println("Weight of cluster #"+currentClusterIndex+" To others="+weight);
			System.out.println("W(X)="+a+", W(P-X)="+b);
			
			
			if (weight<bestWeight && !Util.equal(weight, Util.ZERO))
			{
				bestWeight = weight;
				bestClusterIndex = currentClusterIndex;					
			}
		}// clusterItr		
		
		System.out.println("Clutser num "+bestClusterIndex+" was choosen as best cluster.");
		return clusters.get(bestClusterIndex);
	}
	
	public Vector<Motif> getBestClusterSilhouette(Map<Integer, Vector<Motif>> clusters)
	{
				
		int bestClusterIndex = -1;
		double bestWeight = Integer.MAX_VALUE;
		
		for (Iterator<Integer> clusterItr= clusters.keySet().iterator(); clusterItr.hasNext();)
		{
			int currentClusterIndex = clusterItr.next();
			List<Motif> currentCluster = new ArrayList<Motif>(clusters.get(currentClusterIndex));
			currentCluster.remove(0);	// its center
			List<Motif> restOfClusters = new ArrayList<Motif>();
			for (Iterator<Integer> allClusterItr= clusters.keySet().iterator(); allClusterItr.hasNext();)
			{
				int i = allClusterItr.next();
				if (i!= currentClusterIndex)
				{
					List<Motif> temp = new ArrayList<Motif>(clusters.get(i));
					temp.remove(0);	// its center
					restOfClusters.addAll(temp);
				}
			}
			
			double weight = 0;
			
			for (Iterator<Motif> m=currentCluster.iterator(); m.hasNext();)
				weight+= silhouetteWidth(m.next(), currentCluster, restOfClusters);
			weight = weight / currentCluster.size();
			
			System.out.println("silhouetteWidth #"+currentClusterIndex+"="+weight);
					
			if (weight<bestWeight && !Util.equal(weight, Util.ZERO))
			{
				bestWeight = weight;
				bestClusterIndex = currentClusterIndex;					
			}
		}// clusterItr		
		
		System.out.println("Clutser num "+bestClusterIndex+" was choosen as best cluster.");
		return clusters.get(bestClusterIndex);
	}
	
	/*-----------------------------------------------------------------------------*/

	public Vector<Motif> getBestClusterDist(Map<Integer, Vector<Motif>> clusters)
	{
		int bestClusterIndex = -1;
		double bestDistsqr = Integer.MAX_VALUE;
		int maxNumOfFinders = 0;
		
		for (Iterator<Integer> clusterItr= clusters.keySet().iterator(); clusterItr.hasNext();)
		{
			int currentClusterIndex = clusterItr.next();
			double distSqr = 0;
			Vector<Motif> currentClutser = new Vector<Motif>(clusters.get(currentClusterIndex));
			Motif center = currentClutser.get(0);
			currentClutser.remove(0);
			int numOfFinders = MotifTools.getNumOfFinders(currentClutser);
			for (Iterator<Motif> itr2=currentClutser.iterator(); itr2.hasNext();)
				distSqr += 1.0 - center.sim(itr2.next());
					
			distSqr = distSqr / (double)currentClutser.size();
			
			System.out.println("--------------------");
			System.out.println("Dist sqr of cluster #"+currentClusterIndex+" s="+distSqr+"("+numOfFinders+" finders)");
			
			
			if (distSqr<bestDistsqr && !Util.equalZero(distSqr))
			{
				bestDistsqr = distSqr;
				maxNumOfFinders = numOfFinders;
				bestClusterIndex = currentClusterIndex;		
				System.out.println("Best Cluster changed to: "+bestClusterIndex);
				System.out.println("Max Finders changed to: "+maxNumOfFinders);
			}
		}// clusterItr		
		
		System.out.println("Clutser num "+bestClusterIndex+" was choosen as best cluster.");
		return clusters.get(bestClusterIndex);
	}
	
	/*-----------------------------------------------------------------------------*/

	
	/*-----------------------------------------------------------------------------*/
		
	public double sim(List<Motif> motifs)
	{
		double sum = 0;
		Object[] list = motifs.toArray();
		//int count = 0;
		
		for (int i=0; i<list.length; i++)
		{			
			//Map<Motif, Double> mSim = sim.get(list[i]);
			for (int j=i+1; j<list.length; j++)
			{
				sum += sim.get(list[i]).get(list[j]);				
			}
		}
		
		sum = 2 * sum + motifs.size();	// symmetric matrix ^ 2, plus diagonal sim = 1 * motifs.size()
		
		double total = sum / Math.pow(motifs.size(), 2); // (double)count;
		return total;
	}

	/*
	 * public double sim(List<Motif> motifs)
	{
		double sum = 0;
		List<Motif> list = new ArrayList(motifs);
		//int count = 0;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			list.remove(m);
			Map<Motif, Double> mSim = sim.get(m);
			for (Iterator<Motif> itr2=list.iterator(); itr2.hasNext();)
			{
				Motif temp = itr2.next();
				//if (!m.getFinder().equals(temp.getFinder()))
				{	
					sum += mSim.get(temp);	
					//count++;
				}
			}
		}
		
		sum = 2 * sum + motifs.size();	// symmetric matrix ^ 2, plus diagonal sim = 1 * motifs.size()
		
		double total = sum / Math.pow(motifs.size(), 2); // (double)count;
		return total;
	}
	 */
	/*-----------------------------------------------------------------------------*/
	
	/* Setters and getters */
	public Dataset getDataset()
	{
		return dataset;
	}

	public void setDataset(Dataset dataset)
	{
		this.dataset = dataset;
	}

	public List<Motif> getAccepted()
	{
		return accepted;
	}

	public void setAccepted(Vector<Motif> accepted)
	{
		this.accepted = accepted;
	}

	public Motif getFinalMotif()
	{
		return finalMotif;
	}

	public Map<Integer, Vector<Motif>> getClusters()
	{
		return clusters;
	}
	
	
		
	/*-----------------------------------------------------------------------------*/
	
	
}
