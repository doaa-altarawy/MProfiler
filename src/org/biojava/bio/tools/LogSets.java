package org.biojava.bio.tools;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.hssf.usermodel.HSSFRichTextString;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.motif.Motif;

/*-----------------------------------------------------------------------------*/
public class LogSets
{
	Motif target;
	String root = "log";
	File parent;
	File good;
	File bad;
	double Sn, PPV;
	FileOutputStream excelFile;
	HSSFWorkbook wb = new HSSFWorkbook();
	short rowCount = 0;
	String[] header;
	String outFileName;

	public LogSets(String[] header, String outFileName)
	{
		this.header = header;
		this.outFileName = outFileName;
	}

	public LogSets(Dataset dataset, double Sn, double PPV, String outFileName) throws IOException, BioException
	{
		header = EvalEntry.headerAll;
		this.outFileName = outFileName;
		this.Sn = Sn;
		this.PPV = PPV;
		target = TompaDataset.getAnswer(dataset);
		(new File(root)).mkdir();
		parent = new File(root + File.separator + dataset.name() + "_all_sets");
		parent.mkdir();
		good = new File(parent + File.separator + "good");
		good.mkdir();
		bad = new File(parent + File.separator +"bad");
		bad.mkdir();
		initExcelFile("sheet0");
		File[] fileList = good.listFiles();

		// clear folder
		for (int i=0; i<fileList.length; i++)
			fileList[i].delete();
		fileList = bad.listFiles();

		// clear folder
		for (int i=0; i<fileList.length; i++)
			fileList[i].delete();
	}

	/*********************************************************************************************/

	public void logSet(Motif center, List<Motif> x, double xmean, double xvar, double pmean, double pvar, double separatness, Object weight, boolean saveSet, double num1, double num2) throws Exception
	{
		if (x.size()==0) return;

		MotifEvaluation e = new MotifEvaluation("", target, center);

		if (saveSet)
			if (e.nSn > Sn && e.nPPV > PPV)
			{
				File file = new File(getNextVoterGoodFileName());
				List<Motif> list = new ArrayList<Motif>();
				list.add(center);
				list.addAll(x);
				FileConverterTools.writeMotifVoter(list, file);
			}
			else
			{
				File file = new File(getNextVoterBadFileName());
				List<Motif> list = new ArrayList<Motif>();
				list.add(center);
				list.addAll(x);
				FileConverterTools.writeMotifVoter(list, file);
			}

		EvalEntry values = new EvalEntry(header);

		values.add(rowCount);
		values.add(x.size());
		values.add(MotifTools.getNumOfFinders(x));
		values.add(e.nSn);
		values.add(e.nPPV);
		values.add(e.nCC);
		values.add(e.nPC);
	//	values.add(0);//MotifTools.sim(x));
		double a = 0;//weight(x);
		double b = 0;//weight(p_x);
	//	values.add(a);
	//	values.add(b);
		//values.add(0);//a/b;
		values.add(xmean);
		values.add(xvar);
		values.add(pmean);
		values.add(pvar);
		values.add(separatness);
		values.add(weight);
		values.add(num1);
		values.add(num2);

		addToExcel(values);
	}

	public void logSet(EvalEntry values) throws IOException
	{
		addToExcel(values);
	}
	public void logSet(Object[] values) throws IOException
	{
		addToExcel(values);
	}

	public void logSheet(int sheet, Object[][] values) throws IOException
	{
		initExcelFile(sheet+"");

		for (int i=0; i<values.length; i++)
			addToExcel(sheet, values[i]);
	}


	public void initExcelFile(String sheetName)
	{
		// create a new sheet
		HSSFSheet sheet = wb.createSheet(sheetName);
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

		rowCount = 0;
		// create header row
		for (short i=0; i<header.length; i++)
		{
			row = sheet.createRow(rowCount);
			cell = row.createCell(i);
			cell.setCellValue(new HSSFRichTextString(header[i]));
			cell.setCellStyle(headerStyle);
		}
		rowCount++;
	}

	private void addToExcel(Object[] values) throws IOException
	{
		addToExcel(0, values);
	}

	public void addToExcel(int sheetNum, Object[] values) throws IOException
	{
		// create a new sheet
		HSSFSheet sheet = wb.getSheetAt(sheetNum);
		HSSFRow row = sheet.createRow(rowCount);
		HSSFCell cell = null;

	    for (short cellnum = 0; cellnum < values.length; cellnum++)
	    {
	        cell = row.createCell(cellnum);
	        Object obj = values[cellnum];
	    	if (obj instanceof String)
	    		cell.setCellValue((String)obj);
	    	else if (obj instanceof Integer)
		    		cell.setCellValue(new Double(((Integer)obj).doubleValue()));
	    	else if (obj instanceof Double)
	    		cell.setCellValue((Double)obj);
	    	else if (obj instanceof Float)
	    		cell.setCellValue(new Double(((Float)obj).doubleValue()));
	    	else if (obj instanceof Short)
	    		cell.setCellValue(new Double(((Short)obj).doubleValue()));

	    }
	    rowCount++;
	}
	private void addToExcel(EvalEntry values) throws IOException
	{
		// create a new sheet
		HSSFSheet sheet = wb.getSheetAt(0);
		HSSFRow row = sheet.createRow(rowCount);
		HSSFCell cell = null;

	    for (short cellnum = 0; cellnum < values.size(); cellnum++)
	    {
	        cell = row.createCell(cellnum);
	        Object obj = values.get(cellnum);
	    	if (obj instanceof String)
	    		cell.setCellValue((String)obj);
	    	else if (obj instanceof Integer)
		    		cell.setCellValue(new Double(((Integer)obj).doubleValue()));
	    	else if (obj instanceof Double)
	    		cell.setCellValue((Double)obj);
	    	else if (obj instanceof Float)
	    		cell.setCellValue(new Double(((Float)obj).doubleValue()));
	    	else if (obj instanceof Short)
	    		cell.setCellValue(new Double(((Short)obj).doubleValue()));
	    }
	    rowCount++;
	}

	public void writeToFile() throws IOException
	{
		writeToFile(parent + "\\" + target.getDataset().name()+ "_"+outFileName +".xls");
	}

	public void writeToFile(String fileName) throws IOException
	{
		excelFile = new FileOutputStream(fileName);
		wb.write(excelFile);
		excelFile.close();
	}

	private String getNextVoterGoodFileName()
	{
		return good + "\\" + target.getDataset().name() + "_" + rowCount + ".set.voter";
	}

	private String getNextVoterBadFileName()
	{
		return bad + "\\" + target.getDataset().name() + "_" + rowCount + ".set.voter";
	}

}
