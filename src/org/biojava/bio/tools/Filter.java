package org.biojava.bio.tools;

import java.io.File;
import java.io.FilenameFilter;

import javax.swing.filechooser.FileFilter;

public class Filter extends FileFilter implements FilenameFilter
{
    String type;
    
    public Filter(String type)
    {
 	   this.type = "."+type;
    }
    
    public boolean accept(File file) 
    {
 	   if (file.isDirectory())
 		      return true;
 		   
 	   return file.getName().endsWith(type);
    }
    public String getDescription() {
        return "*"+type;
    }

	public boolean accept(File dir, String name)
	{
		return name.endsWith(type);
	}
}
