package org.biojava.bio.tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;

/**
 * Evaluation Entry to be saved in Excel file
 * @author Doaa1
 *
 */
public class EvalEntry
{
	public static final String[] headerAll_old = {"Set #", "Size", "# of Finders", "nSn", "nPPV", "nCC", "nPC", 
		"Sim(X)", "w(X)", "w(p_x)", "A(X)", "X_Mean", "X_Var", "P_Mean", "P_Var", 
		"Separateness", "F(x)", "init set"};
	
	public static final String[] headerAll = {"Set #", "Size", "# of Finders", "nSn", "nPPV", "nCC", "nPC", 
		"X_Mean", "X_Var", "#sites", "width", 
		"Separateness", "F(x)", "PWM IC", "MAPScore"};
	
		
	public static final String[] headerSubTotal = {"Method", "nCC Threshold", "Avg # of Finders", ">  nSn",
		">  nPPV", "# of sets > nCC", "> nPC", "Avg X_Mean", "Avg X_Var", "Avg P_Mean", "Avg P_Var"};
	
	public static final String[] headerTotal = {"Method", "nCC Threshold",  "# of sets > nCC", "XMean", "XVar", " > nPC", " > nSn", " > nPPV"};
	
	public Vector<Object> data;
	
	public EvalEntry(String[] header)
	{
		data = new Vector<Object>(header.length);
	}
	
	public void add(Object value)
	{
		data.add(value);
	}
	
	public Object get(int i)
	{
		return data.get(i);		
	}
	
	public int size()
	{
		return data.size(); 
	}
	
	public static Object[][] read(File file) throws Exception
	{
		return read(0, file);
	}
		
	public static Object[][] read(int sheetNum, File file) throws Exception
	{
		InputStream input = new FileInputStream(file);
        POIFSFileSystem fs = new POIFSFileSystem( input );
        HSSFWorkbook wb = new HSSFWorkbook(fs);
        HSSFSheet sheet = wb.getSheetAt(sheetNum);
        Object tableData[][] = new Object[sheet.getLastRowNum()+1][sheet.getRow(0).getLastCellNum()+1];
        
        int i=0;
        Iterator rows = sheet.rowIterator();
        rows.next(); // skip header
        for (; rows.hasNext(); i++) 
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
}
