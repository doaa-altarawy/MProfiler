/**
 * 
 */
package org.biojava.bio.view;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;
import org.biojava.bio.BioException;
import org.biojava.bio.charts.BarChart;
import org.biojava.bio.charts.LinesSeriesChart;
import org.biojava.bio.motif.MotifFinders;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.motif.MotifVoter;
import org.biojava.bio.motif.MotifFinders.Method;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.DatasetType;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.constants.OrganismCode;
import org.biojava.bio.tools.EvalEntry;
import org.biojava.bio.tools.LogSets;
import org.biojava.bio.tools.MotifEvaluation;
import org.biojava.bio.tools.FileConverterTools;
import org.biojava.bio.tools.MotifTools;
import org.biojava.bio.tools.RunExternalExe;
import org.biojava.bio.tools.Sortable;
import org.biojava.bio.tools.TompaDataset;
import org.biojava.bio.tools.MotifEvaluation.MEAURES;
import org.biojava.bio.tools.Filter;
import org.biojava.bio.constants.FileNames;
import org.biojava.utils.ChangeVetoException;

import clustering.techniques.denistyBased.MotifDataObject;
import clustering.techniques.denistyBased.MyDatabase;

/**
 * @author Doaa Altarawy
 *
 */
public class EvaluationController
{	
			
	public String log;
	public int numOfEvaluations = 1;
	public static int MAX_NUM_OF_MOTIFS_PER_FIDER = 30;
			
	public void runAndEvaluate(MotifFinder finder, String fileName, int numOfMotifs, int motifWidth, int numOfTrials, OrganismCode organism, Dataset dataset, String outFileDir) throws Exception
	{
		String out = RunExternalExe.runMotifFinder(finder, fileName, numOfMotifs, motifWidth, numOfTrials, organism, dataset, outFileDir);			
		evaluateTompaMotif(finder, out);
		
	}
	
	public void convertFilesAndMerge(File dir, Dataset dataset, boolean readReverse) throws ChangeVetoException, IOException, BioException
	{
		if (!dir.isDirectory())			
			return;
		File[] files = dir.listFiles(new Filter("voter"));
		
		for (int i=0; i<files.length; i++)
			files[i].delete();
		
		files = dir.listFiles(new Filter("al_out"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractAligACE(files[0], dataset, readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("annspec_out"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractANNSpec(files[0], readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("bpout"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractBP(files[0], readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("fasta.wee"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractWeeder(files[0], dataset, readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("fn"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractSPACE(files[0]);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("mds_out"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractMDScan(files[0], dataset, readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("meme_out"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractMEME(files[0], dataset, readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("msp_out"));
		if (files!=null && files.length>0)
			try{FileConverterTools.extractMotifSampler(files[0], dataset, readReverse);}
			catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("txt"));
		for (int i=0; i<files.length; i++)
		{
			if (files[i].getName().contains("final-output"))
				try{FileConverterTools.extractMITRA(files[i], readReverse);}
				catch (Exception e) {e.printStackTrace();}			
		}
		try{
			String outFileName = FileConverterTools.mergeFiles(dir.listFiles(new Filter("voter")), dataset);
			removeRedundantSites(outFileName, dataset);
		}
		catch (Exception e) {e.printStackTrace();}
		
		files = dir.listFiles(new Filter("txt"));
		for (int i=0; i<files.length; i++)
		{
			if (!files[i].getName().contains("final-output"))
				try{FileConverterTools.readMotifVoterWeb(files[i], readReverse, false);}
				catch (Exception e) {e.printStackTrace();}
		}
	}
	/**
	 * get the scores of each given motif alone
	 * @param name
	 * @param inputFileName
	 * @return
	 * @throws Exception 
	 */
	public void evaluateTompaMotif(MotifFinder finder, String inputFileName) throws Exception
	{
		String output = finder.name()+"\n-------------------\n";
		
		Vector<Motif> predicted = FileConverterTools.readTompa(new File(inputFileName), null);
			//String dataset;			
			//dataset= predicted.firstElement().getDataset().name();
			//dataset = dataset.substring(0, dataset.length()-1);	// remove g,r or m
			
			predicted.firstElement().setFinder(finder);
			output += evaluateMotifList(predicted, predicted.firstElement().getDataset());
			
				
		log = output;		
	}
	
	public static String evaluateMotifList(Collection<Motif> motifs, Dataset dataset) throws Exception
	{
		Motif target = TompaDataset.getAnswer(dataset);
		return evaluateMotifList(motifs, target); 
	}
	
	public static String evaluateMotifList(Collection<Motif> motifs, Motif target) throws Exception
	{
		return evaluateMotifList(motifs, target, 0);
	}
	public static String evaluateMotifList(Collection<Motif> motifs, Motif target, int evalNum) throws Exception
	{
		String output = "";
		
		MotifEvaluation[] e = new MotifEvaluation[motifs.size()];
						
		int i = 0;
		for (Iterator<Motif> itr = motifs.iterator(); itr.hasNext(); i++)
		{
			Motif predict = itr.next();
			String title = (predict.getFinder()==null)?"":predict.getFinder().name();
			e[i] = new MotifEvaluation((i+1)+"-"+title, target, predict);
			
			output += e[i].toString();
		}			
		//System.out.println(output);
		
		
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(evalNum));
		
		createPredefinedCharts(e, evalNum);
		
		return output;
	}
	
	public void evaluateTompaMotifEach(String inputFileName) throws Exception
	{
		String output = "";
		
		Vector<Motif> predicted = FileConverterTools.readTompa(new File(inputFileName), null);
		MotifEvaluation[] e = new MotifEvaluation[predicted.size()];
		
		//String dataset;			
		//dataset= predicted.firstElement().getDataset().name();
		//dataset = dataset.substring(0, dataset.length()-1);	// remove g,r or m
					
		Motif target = TompaDataset.getAnswer(predicted.firstElement().getDataset());
		
		int i = 0;			
		for (Iterator<Motif> itr = predicted.iterator(); itr.hasNext(); i++)
		{				
				
			e[i] = new MotifEvaluation("Num "+(i+1), target, itr.next());
							
			output += e[i].toString();
		}			
		System.out.println(output);
		
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		
		createPredefinedCharts(e);
			
				
		log = output;
	}
	/** Evaluate accumulative scores for:
	 * i = 0
	 * i = 0+1
	 * i = 0+1+2
	 * .
	 * .
	 * assuming all motifs are for the same Dataset
	 * @param name
	 * @param inputFileName
	 * @return true if succeeded
	 * @throws Exception 
	 */	
	public void evaluateTompaMotifAcc(String inputFileName) throws Exception
	{
		String output = "";
		
		Vector<Motif> predicted = FileConverterTools.readTompa(new File(inputFileName), null);
		MotifEvaluation[] e = new MotifEvaluation[predicted.size()];
		
		//String dataset;			
		//dataset= predicted.firstElement().getDataset().name();
		//dataset = dataset.substring(0, dataset.length()-1);	// remove g,r or m
					
		Motif target = TompaDataset.getAnswer(predicted.firstElement().getDataset());
		
		int i = 0;
		Motif totalMotif = new Motif(predicted.firstElement().getDataset());
		for (Iterator<Motif> itr = predicted.iterator(); itr.hasNext(); i++)
		{
			totalMotif.merge(itr.next());
				
			e[i] = new MotifEvaluation("Top "+(i+1), target, totalMotif);
							
			output += e[i].toString();
		}			
		System.out.println(output);
		
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		
		createPredefinedCharts(e);			
			
				
		log = output;
	}
	
	/**
	 * get scores of the sum of all given motifs
	 * @param name
	 * @param out
	 * @return
	 * @throws Exception 
	 */
	public void evaluateTompaMotifTotal(String inputFileName) throws Exception
	{
		String output = "";
		
		Map<Dataset, Motif> predicted = FileConverterTools.readTompaTotals(new File(inputFileName));
		MotifEvaluation[] e = new MotifEvaluation[1];
		MotifEvaluation temp;
		
		
		for (Iterator<Dataset> itr = predicted.keySet().iterator(); itr.hasNext();)
		{
			Dataset dataset = itr.next();			
			Motif predict = predicted.get(dataset);
			Motif target = TompaDataset.getAnswer(dataset);
			if (e[0] == null)
			{
				e[0] = new MotifEvaluation("Total", target, predict);
				output += e[0].toString();
			}
			else
			{
				temp = new MotifEvaluation("Total" , target, predict);
				e[0].mergeDatasetEvaluation(temp);
				output += temp.toString();
			}				
		}	
		
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		
		output += e[0].toString();
		System.out.println(output);
		
		createPredefinedCharts(e);
			
				
		log = output;
	}

	public List<Motif> runMotifVoter(String file, Dataset dataset, boolean multipleDataSets, Method type, int minpts, double epsilon, int numOfClusters) throws Exception
	{
				
		if (!multipleDataSets)
			return runMotifVoter(file, file.substring(0, file.lastIndexOf('.'))+ ".doaa.voter", dataset, type, minpts, epsilon, numOfClusters);
		
		File dir = new File(file);
		File[] files = dir.listFiles(new Filter("voter"));
		
		File parent = new File(dir.getCanonicalPath()+ "\\" +type.name());
		parent.mkdir();
		for (int i=0; i<files.length; i++)
		{			
			String datasetName = files[i].getName();
			if (datasetName.contains("_"))
				datasetName = datasetName.substring(0, datasetName.indexOf("_"));
			else 
				datasetName = datasetName.substring(0, datasetName.indexOf("."));
			Dataset dataset2 = Dataset.getValue(datasetName);
			String fileName = parent+ "\\" + files[i].getName().substring(0, files[i].getName().lastIndexOf('.'))+ ".doaa.voter";
			runMotifVoter(files[i].getCanonicalPath(), fileName, dataset2, type, minpts, epsilon, numOfClusters);
		}
		
		return null;
	}
	
	public List<Motif> runMotifVoter(String infile, String out, Dataset dataset, MotifFinders.Method type, int minpts, double epsilon, int numOfClusters) throws Exception
	{
		System.out.println("Running MV for input file:" + infile);
		System.out.println("------------------------------------");

		List<Motif> list = new ArrayList<Motif>();
				
		Vector<Motif> motifs = FileConverterTools.readMotifVoter(new File(infile), false);
		Motif target = TompaDataset.getAnswer(dataset);
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			m.setDataset(dataset);
			m.setCorrect(m.countIntersect(target));			
		}
		
		
		MotifFinders motifVoter = new MotifFinders(dataset, motifs, minpts, epsilon, type, numOfClusters);
				
		motifVoter.findMotif();
		
		if (motifVoter.getFinalMotif()!=null)
			list.add(motifVoter.getFinalMotif());
		if (motifVoter.getAccepted()!=null) 
			list.addAll(motifVoter.getAccepted());
					
		FileConverterTools.writeMotifVoter(list, new File(out));
		if	(type==Method.DBScan || type==Method.Hierarcial || type==Method.KMeans)// Clusters
		{
			Map<Integer, Vector<Motif>> clusters = motifVoter.getClusters();
			for (Iterator<Integer> itr = clusters.keySet().iterator(); itr.hasNext();)
			{
				int n = itr.next();
				List<Motif> result = clusters.get(n);
				if (result!=null) evaluateMotifList(result, target, n);
			}
			numOfEvaluations = motifVoter.getClusters().size();
		}
		else
		{									
			evaluateMotifList(list, target, 0);
		}
		
		list.add(target);
		return list;
	}
	
	public List<Motif> removeRedundantSites(String fileName, Dataset dataset) throws Exception
	{
		Vector<Motif> motifs = FileConverterTools.readMotifVoter(new File(fileName), false);
		Vector<Motif> clean = new Vector<Motif>(motifs.size());
		Motif target = TompaDataset.getAnswer(dataset);
		Map<MotifFinder, Integer> count = new HashMap<MotifFinder, Integer>(MotifFinder.values().length);
		Map<MotifFinder, Motif> finderMotif = new HashMap<MotifFinder, Motif>(MotifFinder.values().length);
		
		for (int i=0; i<MotifFinder.values().length; i++)
		{
			count.put(MotifFinder.values()[i], 0);
			finderMotif.put(MotifFinder.values()[i], new Motif(MotifFinder.values()[i]));
		}
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();
			if (count.get(m.getFinder())> MAX_NUM_OF_MOTIFS_PER_FIDER)
				continue;
			//if (finderMotif.get(m.getFinder()).sim(m)>0.9)	// skip very similar motifs of the same finder
				//continue;			
			count.put(m.getFinder(), count.get(m.getFinder())+1);
			m.setDataset(dataset);
			m.removeRedundant();
			finderMotif.get(m.getFinder()).merge(m);
			clean.add(m);			
		}
				
		
		String outFileName = fileName.substring(0, fileName.indexOf('.'))+ ".clean.voter";						
		FileConverterTools.writeMotifVoter(clean, new File(outFileName));
		evaluateMotifList(clean, target, 0);
		
		return clean;
	}
	
	public List<Motif> extractSites(String fileName, Dataset dataset) throws Exception
	{
		Vector<Motif> motifs = FileConverterTools.readMotifVoter(new File(fileName), false);
		motifs.get(0).setDataset(dataset);
		Vector<Motif> list = new Vector<Motif>(motifs.size());
		Motif target = TompaDataset.getAnswer(dataset);
		
		list.add(MotifTools.extractSitesMultiFinder(motifs));
		list.addAll(motifs);
		
		String outFileName = fileName.substring(0, fileName.indexOf('.'))+ ".extract.voter";						
		FileConverterTools.writeMotifVoter(list, new File(outFileName));
		evaluateMotifList(list, target, 0);
		
		list.add(target);
		return list;
	}
	
	public List<Double> computeSimilarity(String fileName, Dataset dataset) throws Exception
	{
		List<Double> list = new ArrayList<Double>(2);
		Vector<Motif> motifs = FileConverterTools.readMotifVoter(new File(fileName), false);
				
		MotifFinders motifVoter = new MotifFinders(dataset, motifs);
		//motifVoter.computeSortedSimilarity(motifs);
		double sim = motifVoter.compactness_MaxDist(motifs);//motifVoter.sim(motifs);
		motifVoter.computeSortedSimilarity(motifs);
		double weight = motifVoter.weight(motifs);
		
		list.add(sim);
		list.add(weight);
		return list;
	}
			
	public String drawKDistanceGraph(File file, Dataset dataset, int k) throws IOException, BioException
	{	
		Vector<Motif> motifs = FileConverterTools.readMotifVoter(file, false);
		
		MyDatabase<Motif> db = new MyDatabase<Motif>(motifs);
		
		double[][] dist = new double[k][motifs.size()];
		
		
		dist = db.getK_Distance(k);
		
		LinesSeriesChart chart = new LinesSeriesChart();
		String[] series = new String[k];
		Double[] x_values = new Double[motifs.size()];
		Double[][] y_values = new Double[k][motifs.size()];
		
		for (int i=0; i<k; i++)
		{
			series[i] = (i+1)+"-Dist";
				
			for (int j=0; j<dist[i].length; j++)
			{				
				x_values[j] = (double)j;
				y_values[i][j] = dist[i][j];
			}
		}
		chart.createXYChartImage("K-Dist", "Motif num", "Distance", FileNames.LINE_CHART, series, x_values, y_values);
	
		return FileNames.LINE_CHART;
	}
	
	public Vector<Motif> evaluateMotifVoterEach(File inputFile, Dataset dataset) throws Exception
	{
		String output = "";
		
		Vector<Motif> predicted = FileConverterTools.readMotifVoter(inputFile, false);
		MotifEvaluation[] e = new MotifEvaluation[predicted.size()];
		
		//String dataset;			
		//dataset= predicted.firstElement().getDataset().name();
		//dataset = dataset.substring(0, dataset.length()-1);	// remove g,r or m
					
		Motif target = TompaDataset.getAnswer(dataset);
		
		int i = 0;			
		for (Iterator<Motif> itr = predicted.iterator(); itr.hasNext(); i++)
		{				
			Motif m = itr.next();
			e[i] = new MotifEvaluation((i+1)+"-"+m.getFinder().name(), target, m);
							
			output += e[i].toString();
		}			
		System.out.println(output);
		
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		
		createPredefinedCharts(e);
		
		log = output;
		
		predicted.add(target);
		return predicted;			
		
	}
	
	public Vector<Motif> evaluateMotifVoterTotal(File inputFile, Dataset dataset) throws Exception
	{
		String output = "";
		
		Vector<Motif> predicted = FileConverterTools.readMotifVoter(inputFile, false);
		//for (Iterator<Motif> i=predicted.iterator(); i.hasNext();)
		//	predicted.firstElement().merge(i.next());
		predicted.firstElement().setDataset(dataset);
		predicted.set(0, MotifTools.extractSitesMultiFinder(predicted)); /////////////
		
		MotifEvaluation[] e = new MotifEvaluation[1];
		
		//String dataset;			
		//dataset= predicted.firstElement().getDataset().name();
		//dataset = dataset.substring(0, dataset.length()-1);	// remove g,r or m
					
		Motif target = TompaDataset.getAnswer(dataset);
		
		
		
		e[0] = new MotifEvaluation(predicted.firstElement().getFinder().name(), target, predicted.firstElement());
					
		System.out.println(output);
		
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		
		createPredefinedCharts(e);
		
		log = output;
		
		Vector<Motif> motifs = new Vector<Motif>();
		motifs.add(predicted.firstElement());
		motifs.add(target);
		return motifs;			
		
	}
	
	/**
	 * Evaluate all .voter files in the given directory
	 * dataset is taken from each filename.
	 * each file produces one evaluation of the first element in the file only 
	 * @param inputDir
	 * @return
	 * @throws Exception 
	 */
	public MotifEvaluation[] evaluateMotifVoterDirFirstEle(File[] files, OrganismCode org, DatasetType type, boolean readReverse) throws Exception
	{
		String output = "";
		
		Vector<Dataset> datasets = Dataset.getByOrganismAndType(org, type);	
		Map<Dataset, MotifEvaluation> eval = new HashMap<Dataset, MotifEvaluation>(datasets.size()+1);
		//for (int i=0; i<datasets.length; i++)
		//	flag.put(datasets[i], false);
		//MotifEvaluation[] e = new MotifEvaluation[datasets.length+1];
		
		String datasetCode = datasets.firstElement().toString();
		datasetCode = datasetCode.substring(0, datasetCode.length()-3);
		eval.put(Dataset.unknown, new MotifEvaluation("Total "+datasetCode));					
		
		int n;
		for (n=0; n<files.length; n++)
		{
			if (files[n].isDirectory()) continue;
			String datasetName = files[n].getName();	
			if (datasetName.contains("_"))
				datasetName = datasetName.substring(0, datasetName.indexOf('_'));
			else datasetName = datasetName.substring(0, datasetName.indexOf('.'));
			Dataset dataset = Dataset.getValue(datasetName);
			
			Motif predicted;
			if (datasets.contains(dataset))
			{
				Vector<Motif> temp = FileConverterTools.readMotifVoter(files[n], readReverse); 
				if (temp!=null && temp.size()>0) predicted = temp.firstElement();
				else continue;
			}
			else
				continue;
			Motif target = TompaDataset.getAnswer(dataset);
			
			eval.put(dataset, new MotifEvaluation(dataset.name(), target, predicted));
			eval.get(Dataset.unknown).mergeDatasetEvaluation(eval.get(dataset));
			output += eval.get(dataset).toString();
									
			
		}
		output += eval.get(Dataset.unknown).toString();
		
		MotifEvaluation[] e = new MotifEvaluation[datasets.size()+1];
		e[datasets.size()] = eval.get(Dataset.unknown);
		int i=0;
		for (Iterator<Dataset> itr=datasets.iterator(); itr.hasNext(); i++)
		{
			Dataset d = itr.next();
			if (eval.get(d)== null)
			{
				Motif target = TompaDataset.getAnswer(d);				
				e[i] = new MotifEvaluation(d.name(), target, new Motif(d));
				e[datasets.size()].mergeDatasetEvaluation(e[i]);
			}
			else
				e[i] = eval.get(d);
		}
		
		System.out.println(output);	
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		log = output;
		createPredefinedCharts(e);
		
		return e;
	}
	
	/**
	 * evaluates all 4 organisms in Tompa et. al benchmark.
	 * It takes the final output of a method, and claculate its evaluation (nCC,...)
	 * All outputs are written in an excel file and also charts are displayed
	 * @param files: all .voter files of the 4 species
	 * @param type: generic, marcov or real
	 * @param readReverse
	 * @throws Exception
	 */
	public void evaluateMotifVoterDirAllDataset(File[] files, DatasetType type, boolean readReverse) throws Exception
	{
		OrganismCode[] organisms = {OrganismCode.DM, OrganismCode.HS, OrganismCode.MM, OrganismCode.SC};
		MotifEvaluation[][] evals = new MotifEvaluation[organisms.length][];
		
		int count = 0;
		
		for (int i=0; i<organisms.length; i++)
		{
			evals[i] = evaluateMotifVoterDirFirstEle(files, organisms[i], type, readReverse);
			count+= evals[i].length;			
		}
		
		MotifEvaluation[] totalEvals = new MotifEvaluation[count+1];
		int n=0, j;
		totalEvals[totalEvals.length-1] = new MotifEvaluation("Total");
		for (int i=0; i<evals.length; i++)
		{
			for (j=0; j<evals[i].length-1; j++)
			{
				totalEvals[n++] = evals[i][j];
			}
			totalEvals[totalEvals.length-evals.length+i-1] = evals[i][j];
			totalEvals[totalEvals.length-1].mergeDatasetEvaluation(evals[i][j]);
		}
		
		FileConverterTools.toExcelFile(totalEvals, FileNames.getStatisticsExcelFileName(0));
		
		createPredefinedCharts(Arrays.copyOfRange(totalEvals, totalEvals.length-evals.length-1, totalEvals.length));
	}
	/**
	 * Evaluate all .voter files in the given directory
	 * dataset is taken from each filename.
	 * each file produces one evaluation that is the MEREGE of all its motifs 
	 * @param inputDir
	 * @return
	 * @throws Exception 
	 */
	
	public void evaluateMotifVoterDirMerged(File[] files, OrganismCode org, DatasetType type, boolean readReverse) throws Exception
	{
		String output = "";
		
		Vector<Dataset> datasets = Dataset.getByOrganismAndType(org, type);	
		Map<Dataset, MotifEvaluation> eval = new HashMap<Dataset, MotifEvaluation>(datasets.size()+1);
		//for (int i=0; i<datasets.length; i++)
		//	flag.put(datasets[i], false);
		//MotifEvaluation[] e = new MotifEvaluation[datasets.length+1];
		
		int n;
		for (n=0; n<files.length; n++)
		{
			String datasetName = files[n].getName();	
			if (datasetName.contains("_"))
				datasetName = datasetName.substring(0, datasetName.indexOf('_'));
			else datasetName = datasetName.substring(0, datasetName.indexOf('.'));
			Dataset dataset = Dataset.getValue(datasetName);
			
			Vector<Motif>  predicted;
			if (datasets.contains(dataset))
			{
				predicted = FileConverterTools.readMotifVoter(files[n], readReverse);
				for (Iterator<Motif> itr=predicted.iterator(); itr.hasNext();)
					predicted.firstElement().merge(itr.next());				
			}
			else
				continue;
			Motif target = TompaDataset.getAnswer(dataset);
			
			eval.put(dataset, new MotifEvaluation((n+1)+"-"+dataset, target, predicted.firstElement()));
			if (eval.get(Dataset.unknown) == null)
			{
				eval.put(Dataset.unknown, new MotifEvaluation("Total", target, predicted.firstElement()));					
			}
			else
			{					
				eval.get(Dataset.unknown).mergeDatasetEvaluation(eval.get(dataset));
				output += eval.get(dataset).toString();
			}						
			
		}
		output += eval.get(Dataset.unknown).toString();
		
		MotifEvaluation[] e = new MotifEvaluation[datasets.size()+1];
		e[datasets.size()] = eval.get(Dataset.unknown);
		int i=0;
		for (Iterator<Dataset> itr=datasets.iterator(); itr.hasNext(); i++)
		{
			Dataset d = itr.next();
			if (eval.get(d)== null)
			{
				Motif target = TompaDataset.getAnswer(d);				
				e[i] = new MotifEvaluation(n++ +"-"+d.name(), target, new Motif(d));
				e[datasets.size()].mergeDatasetEvaluation(e[i]);
			}
			else
				e[i] = eval.get(d);
		}
		
		System.out.println(output);	
		FileConverterTools.toExcelFile(e, FileNames.getStatisticsExcelFileName(0));
		log = output;
		createPredefinedCharts(e);
		
		
	}
	
	
	public void mergeComputeSets(File parent) throws Exception
	{
		File[] dir = parent.listFiles();
		double[][] MV = new double[6][EvalEntry.headerTotal.length];
		double[][] Clique = new double[6][EvalEntry.headerTotal.length];
		double[][] InitCenter = new double[6][EvalEntry.headerTotal.length];
		double[][] MProfiler = new double[6][EvalEntry.headerTotal.length];
		
	
		int count = 0;
		
		for (int i=0; i<dir.length; i++)
		{
			if (!dir[i].isDirectory()) continue;
			
			File[] files = dir[i].listFiles(new Filter("xls"));
			File evalFile = null;
			for (int j=0; j<files.length; j++)
				if (files[j].getName().equals("SetsEval.xls"))
					{ evalFile = files[j]; break;}
			if (evalFile==null) continue;
			System.out.println("Current File: "+evalFile);
			count++;
			int numOfSheets = new HSSFWorkbook(new POIFSFileSystem(new FileInputStream(evalFile))).getNumberOfSheets();
			for (int j=0; j<numOfSheets; j++)
			{
				Object[][] entries = EvalEntry.read(j, evalFile);
				String method = (String)entries[0][0];
				
				if (method.contains("MV"))
				{
					for (int n=0; n<MV.length; n++)
					{
						MV[n][0] += (Double)entries[n][5];
						MV[n][1] += (Double)entries[n][7];
						MV[n][2] += (Double)entries[n][8];
						MV[n][3] += (Double)entries[n][6];
						MV[n][4] += (Double)entries[n][3];
						MV[n][5] += (Double)entries[n][4];
											
					}
				}
				else if (method.contains("non_non"))
				{
					for (int n=0; n<MV.length; n++)
					{
						Clique[n][0] += (Double)entries[n][5];
						Clique[n][1] += (Double)entries[n][7];
						Clique[n][2] += (Double)entries[n][8];
						Clique[n][3] += (Double)entries[n][6];
						Clique[n][4] += (Double)entries[n][3];
						Clique[n][5] += (Double)entries[n][4];
						
					}
				}
				else if (method.contains("MAP"))
				{
					for (int n=0; n<MV.length; n++)
					{
						InitCenter[n][0] += (Double)entries[n][5];
						InitCenter[n][1] += (Double)entries[n][7];
						InitCenter[n][2] += (Double)entries[n][8];
						InitCenter[n][3] += (Double)entries[n][6];
						InitCenter[n][4] += (Double)entries[n][3];
						InitCenter[n][5] += (Double)entries[n][4];
					}
				}
				else if (method.contains("MProfiler"))
				{
					for (int n=0; n<MV.length; n++)
					{
						MProfiler[n][0] += (Double)entries[n][5];
						MProfiler[n][1] += (Double)entries[n][7];
						MProfiler[n][2] += (Double)entries[n][8];
						MProfiler[n][3] += (Double)entries[n][6];
						MProfiler[n][4] += (Double)entries[n][3];
						MProfiler[n][5] += (Double)entries[n][4];						
						
					}
				}
				
			}
		}
		
		System.out.println("Count="+count);
		
		for (int n=0; n<MV.length; n++)
		{
			for (int j=0; j<MV[n].length; j++)
			{
				Clique[n][j] = Clique[n][j]/(double)count;
				MV[n][j] = MV[n][j]/(double)count;
				InitCenter[n][j] = InitCenter[n][j]/(double)count;
				MProfiler[n][j] = MProfiler[n][j]/(double)count;
			}
		}
		
		LogSets log = new LogSets(EvalEntry.headerTotal, "");
		log.initExcelFile(0+"");
		for (int i=0; i<MV.length; i++)
		{
			log.addToExcel(0, new Object[]{"MV", (double)(i)/10.0, MV[i][0], MV[i][1], MV[i][2], MV[i][3], MV[i][4], MV[i][5]});
		}
		
		log.initExcelFile(1+"");
		for (int i=0; i<MV.length; i++)
		{
			log.addToExcel(1, new Object[]{"Clique", (double)(i)/10.0, Clique[i][0], Clique[i][1], Clique[i][2], Clique[i][3], Clique[i][4], Clique[i][5]});
		}
		
		log.initExcelFile(2+"");
		for (int i=0; i<MV.length; i++)
		{
			log.addToExcel(2, new Object[]{"MAP", (double)(i)/10.0, InitCenter[i][0], InitCenter[i][1], InitCenter[i][2], InitCenter[i][3], InitCenter[i][4], InitCenter[i][5]});
		}
		
		log.initExcelFile(3+"");
		for (int i=0; i<MV.length; i++)
		{
			log.addToExcel(3, new Object[]{"MProfiler", (double)(i)/10.0, MProfiler[i][0], MProfiler[i][1], MProfiler[i][2], MProfiler[i][3], MProfiler[i][4], MProfiler[i][5]});
		}
		log.writeToFile(parent.getCanonicalPath()+"\\SetsTotalEval.xls");
	}
	
	public void compareSetsAllDir(File dir) throws Exception
	{
		File[] list = dir.listFiles();
		
		for (int i=0; i<list.length; i++)
			if (list[i].isDirectory())
				compareSets(list[i]);		
	}
	
	public void compareSets(File dir) throws Exception
	{
		File[] files = dir.listFiles(new Filter("xls"));
		LogSets log = new LogSets(EvalEntry.headerSubTotal, "");
		
		int n = 0;
		for (int i=0; i<files.length; i++)
		{
			if (files[i].getName().contains("SetsEval"))
				continue;
			System.out.println("Current File: "+files[i]);
			Object[][] entries = EvalEntry.read(files[i]);			
			Object[][] eval = calcSetsStatistics(files[i].getName(), entries);
			log.logSheet(n++, eval);
			
		}
		
		log.writeToFile(dir.getCanonicalPath()+"\\SetsEval.xls");
	}
	
	private Object[][] calcSetsStatistics(String name, Object[][] data)
	{
		Object[][] percentEval = new Object[7][EvalEntry.headerSubTotal.length-2];
		String title = name.substring(0, name.indexOf('.'));
		int allDataset = -1;
		
		// for each evaluation percent
		for (int j=0; j<percentEval.length-1; j++)
		{
			double nCC = ((double)(j)/10.0);
			
			for (int n=1; n<percentEval[0].length; n++)
				percentEval[j][n] = 0.0;
			percentEval[j][0] = title;
			percentEval[j][1] =  nCC;
			
			int count = 0;		
			int nSpCount = 0;
			int nPPVCount = 0;
			int nPCCount = 0;
			for (int i=0; i<data.length; i++) // for each row in data
			{
				if (data[i][0]==null) break; // end of data
				
				//if ((j==0 && (Double)data[i][5]<= nCC)||( j>0 && (Double)data[i][5]> nCC))	
				if ((Double)data[i][5]> nCC)
				{									
					if (data[0].length>16 && data[i][16] instanceof String)
					{
						allDataset = i;
						continue;
					}
					count++;
					//percentEval[j][1] =  ((Double)percentEval[j][1]) + ((Double)data[i][1]);
					percentEval[j][2] =  ((Double)percentEval[j][2]) + ((Double)data[i][2]);
					////percentEval[j][3] =  ((Double)percentEval[j][3]) + ((Double)data[i][3]);
					/////percentEval[j][4] =  ((Double)percentEval[j][4]) + ((Double)data[i][4]);
					percentEval[j][5] =  ((Double)percentEval[j][5]) + 1;
					////percentEval[j][6] =  ((Double)percentEval[j][6]) + ((Double)data[i][6]);
					percentEval[j][7] =  ((Double)percentEval[j][7]) + (Double)(data[0].length>16?data[i][11]:data[i][7]);
					percentEval[j][8] =  ((Double)percentEval[j][8]) + (Double)(data[0].length>16?data[i][12]:data[i][8]);
					//percentEval[j][9] =  ((Double)percentEval[j][9]) + ((Double)data[i][13]);
					//percentEval[j][10] =  ((Double)percentEval[j][10]) + ((Double)data[i][14]);					
					
				}	
				if ((Double)data[i][3]> nCC) // nSp
				{								
					nSpCount++;
					percentEval[j][3] =  ((Double)percentEval[j][3]) + 1;				
					
				}	
				if ((Double)data[i][4]> nCC) // nPPV
				{								
					nPPVCount++;
					percentEval[j][4] =  ((Double)percentEval[j][4]) + 1;				
					
				}	
				if ((Double)data[i][6]> nCC) // nPC
				{								
					nPCCount++;
					percentEval[j][6] =  ((Double)percentEval[j][6]) + 1;				
					
				}	
							
			}
			if (count!=0)
			{
				//percentEval[j][1] =  ((Double)percentEval[j][1]) / (double) count;
				percentEval[j][2] =  ((Double)percentEval[j][2]) / (double) count;
				///////percentEval[j][3] =  ((Double)percentEval[j][3]) / (double) count;
				///////percentEval[j][4] =  ((Double)percentEval[j][4]) / (double) count;
				percentEval[j][5] =  ((Double)percentEval[j][5]) / (double) data.length;
				//////percentEval[j][6] =  ((Double)percentEval[j][6]) / (double) count;
				percentEval[j][7] =  ((Double)percentEval[j][7]) / (double) count;
				percentEval[j][8] =  ((Double)percentEval[j][8]) / (double) count;
				//percentEval[j][9] =  ((Double)percentEval[j][9]) / (double) count;
				//percentEval[j][10] =  ((Double)percentEval[j][10]) / (double) count;	
			}
			if (nSpCount!=0) percentEval[j][3] =  ((Double)percentEval[j][3]) / (double) data.length;
			if (nPPVCount!=0) percentEval[j][4] =  ((Double)percentEval[j][4]) / (double) data.length;
			if (nPCCount!=0) percentEval[j][6] =  ((Double)percentEval[j][6]) / (double) data.length;
		}
		
		if (allDataset!=-1)
		{
			percentEval[percentEval.length-1][0] = (String)data[allDataset][16];
			percentEval[percentEval.length-1][1] =  ((Double)data[allDataset][1]);
			percentEval[percentEval.length-1][2] =  ((Double)data[allDataset][2]);
			percentEval[percentEval.length-1][3] =  ((Double)data[allDataset][3]);
			percentEval[percentEval.length-1][4] =  ((Double)data[allDataset][4]);
			percentEval[percentEval.length-1][5] =  ((Double)data[allDataset][5]);
			percentEval[percentEval.length-1][6] =  ((Double)data[allDataset][6]);
			percentEval[percentEval.length-1][7] = ((Double)data[allDataset][11]);
			percentEval[percentEval.length-1][8] = ((Double)data[allDataset][12]);			
		}
		
		return percentEval;
	}
	/**
	 * Max number of series = 14
	 * 
	 * @param eval Motif evaluation class
	 * @param includeValue a 14 booleans indicate presence or absence of the value
	 * @throws IOException 
	 */
	public static void createBarChart(MotifEvaluation[] eval, MotifEvaluation.MEAURES[] measure, String chartName, boolean intValues) throws IOException
	{
		String[] series = new String[measure.length];
		String[] xAxis = new String[eval.length];
		Double[][] values = new Double[series.length][xAxis.length];
		
		
		for (int i=0; i<series.length; i++)
		{
			series[i] = measure[i].name();
			for (int j=0; j<xAxis.length; j++)
			{
				xAxis[j] = eval[j].getName();
				values[i][j] = eval[j].getMeasure(measure[i]);
			}
		}
		
		BarChart chart = new BarChart();
		chart.createChartImage("", chartName, series, xAxis, values, intValues);			
	}

	private static void createPredefinedCharts(MotifEvaluation[] e) throws Exception
	{
		createPredefinedCharts(e, 0);
	}
	
	private static void createPredefinedCharts(MotifEvaluation[] e, int evalNum) throws Exception
	{
		String[] charts = FileNames.getCharts(evalNum);
		
		MEAURES[] m = new MEAURES[]{
				MEAURES.nTP, MEAURES.nFP, MEAURES.nFN  
			};			
			createBarChart(e, m, charts[0], true);
			
			m = new MEAURES[]{
					MEAURES.nSn, MEAURES.nPPV, MEAURES.nCC, MEAURES.nPC  
				};			
			createBarChart(e, m, charts[1], false);
			
			m = new MEAURES[]{
					MEAURES.sTP, MEAURES.sFP, MEAURES.sFN  
				};			
			createBarChart(e, m, charts[2], true);
			
			m = new MEAURES[]{
					MEAURES.sSn, MEAURES.sPPV, MEAURES.sASP  
				};			
			createBarChart(e, m, charts[3], false);
			
			m = new MEAURES[]{
					MEAURES.nSn, MEAURES.nPPV, MEAURES.nCC, MEAURES.nPC, MEAURES.sSn, MEAURES.sPPV, MEAURES.sASP 
				};			
			createBarChart(e, m, charts[4], false);
	}
	
	public String getLog()
	{
		return log;
	}
}
