package org.biojava.bio.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.Vector;

import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.FileNames;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

@SuppressWarnings("deprecation")
public class TompaDataset
{
	static Map<Dataset, Map<String, Site>> datasetsSequences = new HashMap<Dataset, Map<String,Site>>();
	static String pathSep = File.separator;

	public static Motif getAnswer(Dataset dataset) throws IOException, BioException
	{
		Vector<Motif> targets = FileConverterTools.readTompa(new File(FileNames.ANSWER_FILE), dataset.getDatasetType());


		Motif target = null;

		for (Iterator<Motif> itr=targets.iterator(); itr.hasNext();)
		{
			Motif temp = itr.next();
			if (temp.getDataset()==null)
				System.out.print("can't happen!");
			if (temp.getDataset().equals(dataset))
			{
				target = temp;
				break;
			}
		}
		target.setTotalLength(getTotalDatasetLenght(dataset.name()));
		target.setFinder(MotifFinder.Target);

		if (target == null)
			System.out.println("Can't find the dataset:"+dataset.name());


		return target;
	}

	public static Map<String, Site> getSequences(Dataset dataset) throws Exception
	{
		if (datasetsSequences.get(dataset)==null)
		{
			SequenceIterator seqItr = TompaDataset.getDatasetSequences(dataset);
			Map<String, Site> sequences = new HashMap<String, Site>();
			while (seqItr.hasNext())
			{
				Sequence seq = seqItr.nextSequence();
				String source = seq.getName().substring(seq.getName().indexOf('_')+1);
				sequences.put(source, new Site(seq.seqString(), source, -seq.length(), -1, "+"));
			}
			return sequences;
		}
		else
			return datasetsSequences.get(dataset);
	}

	public static SequenceIterator getDatasetSequences(Dataset dataset) throws IOException
	{
		return getSequences(new File(getDatasetSequenceFileName(dataset)));
	}

	public static String getDatasetSequenceFileName(Dataset dataset) throws IOException
	{

		String fileName = FileNames.TOMPA_SEQUENCES + dataset.getDatasetType().name()+ pathSep
		                     + dataset.name().substring(0, dataset.name().length()-3)
		                     + pathSep + dataset.name() + ".fasta";
		return new File(fileName).getCanonicalPath();
	}

	public static SequenceIterator getSequences(File file) throws FileNotFoundException
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		return SeqIOTools.readFastaDNA(reader);
	}

	public static long[] getDatasetLengths(String dataset) throws BioException, IOException
	{
		long[] lens = null;
		Vector<Long> v = new Vector<Long>();
		String temp;
		BufferedReader br = new BufferedReader(new FileReader(FileNames.LENGHTH_FILE));

		dataset = dataset.substring(0, dataset.length()-1);
		boolean found = false;

		while ((temp = br.readLine())!=null)
		{
			if (temp.contains(">data set"))
			{
				found = (temp = br.readLine()).equalsIgnoreCase(dataset);
			}
			else if (temp.contains(">length"))
			{
				if (found)
				{
					temp = br.readLine();
					StringTokenizer st = new StringTokenizer(temp, ",");
					while (st.hasMoreTokens())
						v.add(Long.parseLong(st.nextToken()));
					break;
				}
			}
		}

		lens = new long[v.size()];
		int i = 0;
		for (Iterator<Long> itr=v.iterator(); itr.hasNext();)
			lens[i++] = itr.next().longValue();

		return lens;
	}

	public static long getTotalDatasetLenght(String dataset) throws BioException, IOException
	{
		return getTotalDatasetLenght(getDatasetLengths(dataset));
	}

	public static long getTotalDatasetLenght(long[] lens)
	{
		long total = 0;

		if (lens!=null)
		{
			for (int i=0; i<lens.length; i++)
				total+= lens[i];
		}

		return total;
	}
}
