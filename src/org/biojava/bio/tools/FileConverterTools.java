package org.biojava.bio.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.hssf.usermodel.HSSFRichTextString;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;
import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.DatasetType;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.tools.MotifEvaluation.MEAURES;
import org.biojava.utils.ChangeVetoException;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.biojavax.bio.seq.RichSequence.IOTools;
import org.biojavax.bio.seq.io.FastaHeader;


/**
 * 
 * This class reads different specific data files and convert them to FASTA format
 * or EMBL format.
 * 
 * @author Doaa Altarawy
 *
 */
public class FileConverterTools
{

	/**
	 * This main is for testing only..
	 * @param args
	 */
	public static void mainn(String[] args)
	{
		try{			
			//Vector<Motif> m = readTompa(new File("../tompa.out"));
			//writeTompa(m, new File("../myTest2.out"));
			createLengthFile(new File("E:\\Doaa Store\\Master\\Bioinformatics\\Coding\\Data Sets\\Assessment Data\\MChain"), new File(".\\resources\\length.tompa"));
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
		
	/**
	 * intermediate method to convert the output of Bioprospector to Tompa
	 * @param in
	 * @param out
	 * @throws IOException
	 * @throws BioException
	 */
	public static Vector<Motif> covertBioProspectorToTompaFile(File in, File out, Dataset dataset) throws IOException, BioException
	{
		Vector<Motif> motifs = readTompa(in, null);
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			Motif m = itr.next();			
			if (dataset==null || dataset.equals(Dataset.unknown))
				m.setDataset(Dataset.getValue(m.getComment().substring(m.getComment().lastIndexOf('\\')+1, m.getComment().lastIndexOf('.')-1)));
			else m.setDataset(dataset);
			m.setFinder(MotifFinder.BP);
			
			for (Iterator<String> seqItr=m.getSites().keySet().iterator(); seqItr.hasNext();)
			{
				String seq = seqItr.next();
				List<Site> list = m.getSites().get(seq);
				List<Site> sites = new ArrayList<Site>();
				for (Iterator<Site> site=list.iterator(); site.hasNext();)
				{
					Site s = site.next();
					s.source = s.source.substring(s.source.lastIndexOf('_')+1);					
					sites.add(s);
				}				
				m.getSites().put(seq, sites);
			}
		}
		
		writeTompa(motifs, out);
		return motifs;
	}
	
	/**
	 * Run once
	 * reads all files in the given directory and converts them 
	 * to fasta format
	 * 
	 * @param directory a directory of files
	 * @throws Exception
	 */
	public static void fromParenthesesToFASTA(File directory) throws Exception
	{
		File file[];
		if (directory.isDirectory())
		{
			file = directory.listFiles();
		}
		else 
		{
			file = new File[1];
			file[0] = directory;
		}
		
		for (int i=0; i<file.length; i++)
		{
			BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(file[i]));
			//RichSequenceIterator rsi = IOTools.readFastaDNA(br, null);
			String path = file[i].getCanonicalPath();
			path = path.substring(0, path.lastIndexOf('.')) + ".fasta";
			
			FileOutputStream outStream = new FileOutputStream(path);
			String buff;
			FastaHeader header = new FastaHeader();
					
			while ((buff = br.readLine()) !=null)
			{
				buff = buff.substring(1, buff.length()-2).replaceAll(" ", "");
				while (buff.indexOf('?')!=-1)
				{
					buff = buff.replace(buff.charAt(buff.indexOf('?')), ' ');
				}
								
				Sequence seq = DNATools.createDNASequence(buff, "");
				System.out.println(buff);
				IOTools.writeFasta(outStream, seq, null, header);
			}			
		}
	}
	
	
	/**
	 * Reads a file in Tompa's format
	 * and returns a vector of Motifs ordered as in the file
	 * @param in
	 * @return
	 * @throws IOException 
	 * @throws BioException 
	 */
	public static Vector<Motif> readTompa(File in, DatasetType type) throws IOException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(in));		
		StringTokenizer st;
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains(">data set"))
			{				
				temp = br.readLine();		
				// if non defined dataset, put the value in comment field
				m = new Motif();
				Dataset dataset = Dataset.getValue(temp, type);
				if (dataset==null)
					m.setComment(temp);
				else m.setDataset(dataset);
				if (!(temp = br.readLine()).contains(">instances"))
						throw new BioException("Invalid file format: "+ temp +".");				
				
			}
			else if (temp.contains(">instances"))
			{
				while ((temp = br.readLine())!=null && !temp.contains(">data set"))
				{
					st = new StringTokenizer(temp, ",");
					if (st.countTokens()<3)
					{
						temp = br.readLine();
						break;		// anything not in the correct format is ignored
					}
					String source = st.nextToken();
					int start = Integer.parseInt(st.nextToken());
					String dna = st.nextToken();
					int end = start + dna.length()-1 ;
					m.addSite(dna, source, start, end, "+");
				}
				motifs.add(m);
			}
			else
			{
				br.close();
				throw new BioException("Invalid file format: "+ temp);
			}
		}
		
		br.close();
		return motifs;
	}
	
	/**
	 * Read Tompa file into a map of Motifs: <dataset name, totalMotif>
	 * by merging motifs of the same dataset
	 * @param file
	 * @return
	 * @throws IOException
	 * @throws BioException
	 */
	public static Map<Dataset, Motif> readTompaTotals(File file) throws IOException, BioException
	{
		Vector<Motif> motifs = FileConverterTools.readTompa(file, null);		
		Map<Dataset, Motif> predictedPerDataset = new HashMap<Dataset, Motif>();
		
		for (Iterator<Motif> itr = motifs.iterator(); itr.hasNext();)
		{
			Motif predict = itr.next();
			//String dataset;			
			//dataset= predict.getDataset();
			//dataset = dataset.substring(0, dataset.length()-1);	// remove g,r or m
			
			Motif mt = predictedPerDataset.get(predict.getDataset());
			if (mt == null)
				mt = predict;
			else
				mt.merge(predict);
			predictedPerDataset.put(predict.getDataset(), mt);
		}
		return predictedPerDataset;
	}
	
	/**
	 * Save a vector of Motifs to Tompa format
	 * Not that MotifFinder field is lost in Tompa file format
	 * @param motifs
	 * @param outFile
	 * @throws IOException
	 */
	public static void writeTompa(Vector<Motif> motifs, File outFile) throws IOException
	{
		FileWriter out = new FileWriter(outFile);
		Motif m;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			m = itr.next();
			out.append(">data set\n");
			out.append(m.getDataset()+"\n");
			out.append(">instances\n");
			
			for (Iterator<Site> s=m.getAllSites().iterator(); s.hasNext();)
			{
				Site site = s.next();
				out.append(site.source+","+site.location.getMin()+","+site.dna+"\n");
			}
		}
		
		out.flush();
		out.close();
	}
	
	/**
	 *  Run once to generate file containing:
	 *	for each dataset: lengths of its sequences
	 * @throws IOException 
	 * @throws BioException 
	 * @throws NoSuchElementException 
	 *
	 */
	public static void createLengthFile(File datasetDir, File out) throws IOException, NoSuchElementException, BioException
	{
		File file[];
		FileOutputStream outStream = new FileOutputStream(out);
		
		if (datasetDir.isDirectory())
		{
			file = datasetDir.listFiles();
		}
		else 
		{
			file = new File[1];
			file[0] = datasetDir;
		}
		
		for (int i=0; i<file.length; i++)
		{
			// get dataset name without '.fasta' 
			String datasetName = file[i].getName();
			datasetName = datasetName.substring(0, datasetName.lastIndexOf('.'));
			
			BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(file[i]));
			StringBuffer buff = new StringBuffer(">data set\n");
			buff.append(datasetName+"\n");
			buff.append(">length\n");
			
			
			RichSequenceIterator itr = IOTools.readFastaDNA(br, null);
			while (itr.hasNext())
			{
				buff.append(itr.nextSequence().length());
				if (itr.hasNext())
					buff.append(",");
				else	
					buff.append("\n");	
			}			
			outStream.write(buff.toString().getBytes());			
			br.close();
			
		}
		outStream.close();
	}
	
	/**
	 * Save statistics in the MotifEvaluation class to the given file
	 * @param e MotifEvaluation array
	 * @param fileName
	 * @throws IOException
	 */
	public static void toExcelFile(MotifEvaluation e[], String fileName) throws IOException
	{
		// create a new file
		FileOutputStream out = new FileOutputStream(fileName);
		// create a new workbook
		HSSFWorkbook wb = new HSSFWorkbook();
		// create a new sheet
		HSSFSheet sheet = wb.createSheet();
		// declare a row object reference
		HSSFRow row = null;
		// declare a cell object reference
		HSSFCell cell = null;
		// create 3 cell styles
		HSSFCellStyle headerStyle = wb.createCellStyle();		
		HSSFFont f = wb.createFont();		
		f.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);	
		headerStyle.setFont(f);
		headerStyle.setAlignment(HSSFCellStyle.ALIGN_CENTER);
		
		// create header row
		for (short i=1; i<MotifEvaluation.MEAURES.values().length+1; i++)
		{
			row = sheet.createRow(0);
			cell = row.createCell(i);
			cell.setCellValue(new HSSFRichTextString(MotifEvaluation.MEAURES.values()[i-1].name()));
			cell.setCellStyle(headerStyle);
		}
		
		
		for (short rownum = 1; rownum < e.length+1; rownum++)
		{
		    // create a row
			row = sheet.createRow(rownum);
		    		
			cell = row.createCell((short)0);
			cell.setCellValue(new HSSFRichTextString(e[rownum-1].getName()));
			cell.setCellStyle(headerStyle);
			
		    for (short cellnum = 1; cellnum < MotifEvaluation.MEAURES.values().length+1; cellnum++)
		    {
		        cell = row.createCell(cellnum);		       
		    	cell.setCellValue(e[rownum-1].getMeasure(MEAURES.values()[cellnum-1]));
		    }
		}

		wb.write(out);
		out.close();

	}
	
	/**
	 * Reads an excel file to a 2D array of Objects
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	@SuppressWarnings({ "unchecked", "deprecation" })
	public static Object[][] readExcel(String fileName) throws IOException
	{
		InputStream input = new FileInputStream(fileName);
        POIFSFileSystem fs = new POIFSFileSystem( input );
        HSSFWorkbook wb = new HSSFWorkbook(fs);
        HSSFSheet sheet = wb.getSheetAt(0);
        Object tableData[][] = new Object[sheet.getLastRowNum()+1][sheet.getRow(0).getLastCellNum()+1];
        
        int i=0;
        for (Iterator rows = sheet.rowIterator(); rows.hasNext(); i++ ) 
        {           
            HSSFRow row = (HSSFRow) rows.next();
            tableData[i] = new Object[sheet.getRow(0).getLastCellNum()+1];
            int j=0;
            for(Iterator cells = row.cellIterator(); cells.hasNext(); j++) 
            {
                HSSFCell cell = (HSSFCell) cells.next();               
                switch ( cell.getCellType() ) {
                    case HSSFCell.CELL_TYPE_NUMERIC:
                        tableData[i][j] = cell.getNumericCellValue();
                        break;
                    case HSSFCell.CELL_TYPE_STRING: 
                    	tableData[i][j] = cell.getStringCellValue();
                        break;
                    default:
                        System.out.println( "unsuported sell type" );
                        break;
                }
            }
            
        }
        return tableData;
	}
	
	/**
	 * save motifs in MtifVoter file format
	 * the field "MOTIF consencusString" is left empty, it has no effect
	 * @param motifs
	 * @param outFile
	 * @throws IOException
	 */
	public static void writeMotifVoter(Collection<Motif> motifs, File outFile) throws IOException
	{
		FileWriter out = new FileWriter(outFile);
		Motif m;
		
		for (Iterator<Motif> itr=motifs.iterator(); itr.hasNext();)
		{
			m = itr.next();
			out.append("MOTIF FINDER: "+m.getFinder()+"\n");	
			out.append("MOTIF: c\n");		// no consensus to store
			//out.append((m.getFinder()!=null)?m.getFinder().name():""+"\n");
			out.append("INSTANCES:\n");
			
			for (Iterator<Site> s=m.getAllSites().iterator(); s.hasNext();)
			{
				Site seq = s.next();				
				out.append(seq.source+","+seq.location.getMin()+","+seq.dna+","+seq.getStrandString()+"\n");
			}
			out.append("\n");
		}
		
		out.flush();
		out.close();
	}
	
	/**
	 * Read MotifVoter file format to a vector of motifs
	 * @param inFile
	 * @return
	 * @throws IOException
	 * @throws BioException
	 */
	public static Vector<Motif> readMotifVoter(File inFile, boolean readReverse) throws IOException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("MOTIF FINDER:"))
			{				
				st = new StringTokenizer(temp.substring(temp.indexOf(':')), " ");
				st.nextToken();
				m = new Motif(MotifFinder.valueOf(st.nextToken()));

				//skip the MOTIF consensus line
				if ((temp = br.readLine()).contains("MOTIF:"))
					temp = br.readLine();		
				
				if (!temp.contains("INSTANCES:"))
						throw new BioException("Invalid file format: "+temp);				
				
			}
			else if (temp.contains("INSTANCES:"))
			{
				while ((temp = br.readLine())!=null && !temp.contains("MOTIF FINDER:"))
				{
					st = new StringTokenizer(temp, ",");
					if (st.countTokens()<4)
					{
						temp = br.readLine();
						break;		// anything not in the correct format is ignored
					}
					String source = st.nextToken();
					int start = Integer.parseInt(st.nextToken());
					String dna = st.nextToken();
					int end = start + dna.length()-1 ;
					String type = st.nextToken();
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
				}
				motifs.add(m);
			}
			else
			{
				br.close();
				throw new BioException("Invalid file format: "+ temp);
			}
		}
		
		br.close();
		return motifs;
	}
	
	/**
	 * Read MotifVoter Web output file to a vector of motifs
	 * @param inFile
	 * @return String: outFileName
	 * @throws IOException
	 * @throws BioException
	 */
	public static String readMotifVoterWeb(File inFile, boolean readReverse, boolean option2) throws IOException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("FINAL MotifVoter INSTANCES"))
			{					
				m = new Motif(MotifFinder.MotifVoter);
				while ((temp = br.readLine())!=null && !temp.contains("ALIGNED"))
				{
					st = new StringTokenizer(temp, ",");
					if (st.countTokens()<4)
					{
						temp = br.readLine();
						break;		// anything not in the correct format is ignored
					}
					String source = st.nextToken();
					int start = Integer.parseInt(st.nextToken());
					String dna = st.nextToken();
					int end = start + dna.length()-1 ;
					String type = st.nextToken();
					type = type.contains("-")?"-":"+";
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);				
				
			}
			else if (temp.contains("INSTANCES BY STAND ALONE PROGRAMS"))
			{
				while ((temp = br.readLine())!=null)
				{					
					m = new Motif(MotifFinder.valueOf(temp));
					while ((temp = br.readLine())!=null && !temp.contains("****"))
					{
						st = new StringTokenizer(temp, ",");
						if (st.countTokens()<3)
						{
							temp = br.readLine();
							break;		// anything not in the correct format is ignored
						}
						String source = st.nextToken();
						int start = Integer.parseInt(st.nextToken());
						String dna = st.nextToken();
						String type;
						if (option2)
						{
							type = st.nextToken(); 
						}
						else
						{
							type = dna.contains("-")?"-":"+";
							dna = dna.substring(0, dna.indexOf("\t"));
						}
						
						int end = start + dna.length()-1 ;
						
						if (readReverse || (!readReverse && type.equals("+")))
							m.addSite(dna, source, start, end, type);
					}
					if (m.getAllSites().size()>0)	// don't store empty motif
						motifs.add(m);
				}
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		
		Vector<Motif> motifVoterIn = new Vector<Motif>(motifs);		
		if (motifVoterIn.size()>1)
			motifVoterIn.remove(0);
		File outFile;
		if (option2)
			outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".MVW2.voter");
		else 
			outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".MVW.voter");
		//writeMotifVoter(motifVoterIn, new File(outFile+".voter"));
		//System.out.println("Converted file:" + outFile.getCanonicalPath()+".voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		
		return outFile.getCanonicalPath();
	}
		
	public static Vector<Motif> extractMITRA(File inFile, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("Orientation Start"))
			{
				m = new Motif(MotifFinder.MITRA);
				while ((temp = br.readLine())!=null && !temp.contains("---------") && !temp.contains("Source Information"))
				{
					st = new StringTokenizer(temp, " ");
					if (st.countTokens()<5)
						continue;
					
					String source = st.nextToken();
					source = source.substring(source.indexOf('_')+1);
					String type = st.nextToken();;
					type = type.equals("normal")?"+":"-";
					int start = Integer.parseInt(st.nextToken());
					int end = Integer.parseInt(st.nextToken());
					String dna = st.nextToken();				
					dna = dna.replaceAll("-", "");
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".MITRA.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	public static Vector<Motif> extractBP(File inFile, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains(">"))
			{
				m = new Motif(MotifFinder.BP);
				while (temp!=null && !temp.contains("***********"))
				{
					st = new StringTokenizer(temp, "\t");					
					if (st.countTokens()<5)
					{
						temp = br.readLine();
						continue;
					}
					
					String source = st.nextToken();
					source = source.substring(source.indexOf('_')+1);
					st.nextToken();	//ignore					
					String type = st.nextToken();
					type = type.contains("f")?"+":"-";
					int start = Integer.parseInt(st.nextToken());
					String dna = st.nextToken();
					int end = start + dna.length()-1 ;
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
					temp = br.readLine();
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".bp.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	public static Vector<Motif> extractWeeder(File inFile, Dataset dataset, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		long[] len = TompaDataset.getDatasetLengths(dataset.name());
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("Best occurrences:"))
			{
				m = new Motif(MotifFinder.WEEDER);
				temp = br.readLine(); //skip header
				while ((temp = br.readLine())!=null && !temp.contains("Frequency Matrix"))
				{
					st = new StringTokenizer(temp, " ");					
					if (st.countTokens()<4)
					{
						temp = br.readLine();
						continue;
					}
					
					String source = Integer.parseInt(st.nextToken())-1+"";
					String type = st.nextToken();
					String dna = st.nextToken();
					dna = dna.replace("[", "");
					dna = dna.replace("]", "");
					int start = Integer.parseInt(st.nextToken()) - (int)len[Integer.parseInt(source)] - 1;
					int end = start + dna.length()-1 ;
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);					
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".weeder.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	public static Vector<Motif> extractMEME(File inFile, Dataset dataset, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		long[] len = TompaDataset.getDatasetLengths(dataset.name());
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("sites sorted by position p-value"))
			{
				m = new Motif(MotifFinder.MEME);
				br.readLine(); //skip header
				br.readLine();
				br.readLine();
				while ((temp = br.readLine())!=null && !temp.contains("--------------"))
				{
					st = new StringTokenizer(temp, " ");					
					if (st.countTokens()<5)
					{
						temp = br.readLine();
						continue;
					}
					String source = st.nextToken();
					source = source.substring(source.indexOf('_')+1);
					
					String type = st.nextToken();
					String s;
					if (type.contains("+") || type.contains("-"))
						s = st.nextToken();
					else 
					{
						s = type;
						type = "+";
					}
					int start = Integer.parseInt(s) - (int)len[Integer.parseInt(source)] - 1;					
					st.nextToken(); //skip
					st.nextToken(); //skip
					String dna = st.nextToken();					
					int end = start + dna.length()-1 ;
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);					
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".MEME.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	public static Vector<Motif> extractMDScan(File inFile, Dataset dataset, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		Motif m = null;
		long[] len = TompaDataset.getDatasetLengths(dataset.name());
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains(">"))
			{
				m = new Motif(MotifFinder.MDSCAN);
				while (temp!=null && !temp.contains("Mtf") && !temp.contains("Total time"))
				{
					st = new StringTokenizer(temp, "\t");					
					if (st.countTokens()<3)
					{
						temp = br.readLine();
						continue;
					}
					
					String source = st.nextToken();
					source = source.substring(source.indexOf('_')+1);
										
					String type = st.nextToken();
					String s = type.substring(type.indexOf(" ")+2);
					int start = Integer.parseInt(s) - (int)len[Integer.parseInt(source)] - 1;
					type = type.contains("f")?"+":"-";
					String dna = st.nextToken();
					int end = start + dna.length()-1 ;
										
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
					temp = br.readLine();
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".MDScan.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
			
	public static Vector<Motif> extractAligACE(File inFile, Dataset dataset, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		long[] len = TompaDataset.getDatasetLengths(dataset.name());
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("Motif"))
			{
				m = new Motif(MotifFinder.ALIGNACE);
				while ((temp = br.readLine())!=null && !temp.contains("****"))
				{
					st = new StringTokenizer(temp, "\t");
					if (st.countTokens()<4)
						continue;
					
					String dna = st.nextToken();
					String source = st.nextToken();		
					String s = st.nextToken();
					int start = (Integer.parseInt(s)+1) - (int)len[Integer.parseInt(source)] - 1;
					int end = start + dna.length()-1 ;
					String type = st.nextToken();				
					type = type.equals("1")?"+":"-";
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
				}
				if (m.getAllSites().size()>0)	// don't store empty motif	
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".AlignACE.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	public static Vector<Motif> extractMotifSampler(File inFile, Dataset dataset, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;
		long[] len = TompaDataset.getDatasetLengths(dataset.name());
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("#id"))
			{
				m = new Motif(MotifFinder.MTSAMP);
				while ((temp = br.readLine())!=null && !temp.contains("#id"))
				{
					st = new StringTokenizer(temp, "\t");
					if (st.countTokens()<9)
						continue;
					
					String source = st.nextToken();
					source = source.substring(source.indexOf('_')+1);
					st.nextToken();	// skip
					st.nextToken(); // skip
					String s = st.nextToken();
					int start = Integer.parseInt(s)+ - (int)len[Integer.parseInt(source)] - 1;
					st.nextToken();	// skip
					st.nextToken(); // skip
					String type = st.nextToken();
					st.nextToken(); // skip
					
					String dna = st.nextToken();
					dna = dna.substring(dna.indexOf("site")+6, dna.length()-3);	
					
					int end = start + dna.length()-1 ;
					
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".MTSAMP.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	public static Vector<Motif> extractSPACE(File inFile) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;		
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("Significance Score"))
			{
				m = new Motif(MotifFinder.SPACE);
				while ((temp = br.readLine())!=null && !temp.contains("===="))
				{
					st = new StringTokenizer(temp, ",");
					if (st.countTokens()<3)
						continue;
					
					String source = st.nextToken();
					int start = Integer.parseInt(st.nextToken()) -1;
					String dna = st.nextToken();
					int end = start + dna.length()-1 ;

					m.addSite(dna, source, start, end, "+");
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".SPACE.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	public static Vector<Motif> extractANNSpec(File inFile, boolean readReverse) throws IOException, ChangeVetoException, BioException
	{
		Vector<Motif> motifs = new Vector<Motif>();
		String temp;		
		BufferedReader br = new BufferedReaderIgnoreEmpty(new FileReader(inFile));		
		StringTokenizer st;		
		Motif m = null;
		
		temp = br.readLine();
		while (temp!=null)
		{
			if (temp.contains("INFORMATION_CONTENT"))
			{
				m = new Motif(MotifFinder.ANNSPEC);
				br.readLine();
				while ((temp = br.readLine())!=null && !temp.contains("ALIGNMENT_CONS"))
				{
					st = new StringTokenizer(temp, " ");
					if (st.countTokens()<7)
						continue;
					st.nextToken(); //skip
					st.nextToken();	//skip
					st.nextToken();	//skip
					String type = st.nextToken();
					type = type.contains("'")?"-":"+";
					String dna = st.nextToken();					
					int start = Integer.parseInt(st.nextToken());					
					int end = start + dna.length()-1 ;
					String source = st.nextToken();
					source = source.substring(source.indexOf('_')+1);
					
					// Neglect file if ANNSpec combined 2 sequences
					try{Integer.parseInt(source);}
					catch (Exception e) {return null;}
					
					if (readReverse || (!readReverse && type.equals("+")))
						m.addSite(dna, source, start, end, type);
				}
				if (m.getAllSites().size()>0)	// don't store empty motif
					motifs.add(m);
				temp = br.readLine();				
			}
			else
				temp = br.readLine();
			
		}
		
		br.close();
		File outFile = new File(inFile.getParentFile().getCanonicalPath()+"\\"+inFile.getParentFile().getName()+".ANNSpec.voter");
		writeMotifVoter(motifs, outFile);
		System.out.println("Converted file:" + outFile.getCanonicalPath());
		return motifs;
	}
	
	private static class BufferedReaderIgnoreEmpty extends BufferedReader
	{

		public BufferedReaderIgnoreEmpty(Reader in)
		{
			super(in);			
		}
		
		public String readLine() throws IOException
		{
			String temp = super.readLine();
			while (temp!= null && temp.trim().equals(""))
				temp = super.readLine();
			return temp;			
		}
		
	}
	
	public static String mergeFiles(File[] inFiles, Dataset dataset) throws IOException
	{
		int BSIZE = 1024;		
		String outFileName = inFiles[0].getParentFile().getCanonicalPath()+"\\"+inFiles[0].getParentFile().getName()+".voter";
		File temp = new File(outFileName);
		temp.delete();
		FileChannel  outStream = new FileOutputStream(outFileName).getChannel();
		
		
		for (int i=0; i<inFiles.length; i++)
		{			
			if(inFiles[i].isDirectory()) continue;
			FileChannel in = new FileInputStream(inFiles[i]).getChannel();
		    ByteBuffer buffer = ByteBuffer.allocate(BSIZE);
		    while (in.read(buffer) != -1) {
		      buffer.flip(); // Prepare for writing
		      outStream.write(buffer);
		      buffer.clear(); // Prepare for reading
		    }
		    in.close();
		}
		outStream.close();
		System.out.println("Files Merged into: "+outFileName);
		return outFileName;
	}
	
	public static String mergeFiles(File inFilesDir, Dataset dataset) throws IOException
	{
		return mergeFiles(inFilesDir.listFiles(), dataset);
	}
	
	public static void copy(File sourceFile, File destDir) throws IOException 
	{
	     FileChannel in = null, out = null;
	     try 
	     {          
	          in = new FileInputStream(sourceFile).getChannel();
	          out = new FileOutputStream(destDir+"\\"+sourceFile.getName()).getChannel();
	 
	          long size = in.size();
	          MappedByteBuffer buf = in.map(FileChannel.MapMode.READ_ONLY, 0, size);
	 
	          out.write(buf);	 
	     } 
	     finally 
	     {
	          if (in != null)          
	        	  in.close();
	          if (out != null)     
	        	  out.close();
	     }
	}
}
