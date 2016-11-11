package org.biojava.bio.tools;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class UnZip 
{
   final static int BUFFER = 2048;
   public static String unZip(File inFile, String oldExtenstion, String newExtention) 
   {
	  String fileName = null;
	  try
      {    	 
         BufferedOutputStream dest = null;
         FileInputStream fis = new FileInputStream(inFile);
         ZipInputStream zis = new ZipInputStream(new BufferedInputStream(fis));
         ZipEntry entry;
         while((entry = zis.getNextEntry()) != null) 
         {
            if (!entry.getName().contains(oldExtenstion))
            	continue;
        	System.out.println("Extracting: " +entry);
            int count;
            byte data[] = new byte[BUFFER];
            // write the files to the disk
            fileName = inFile.getParentFile()+"\\"+entry.getName()+newExtention;
            FileOutputStream fos = new FileOutputStream(fileName);
            dest = new BufferedOutputStream(fos, BUFFER);
            while ((count = zis.read(data, 0, BUFFER))!= -1) 
               dest.write(data, 0, count);
            
            dest.flush();
            dest.close();
         }
         zis.close();
      } catch(Exception e) {e.printStackTrace(); return null;}
      
      return fileName;
   }
   
   public static void unZip (File inFile) 
   {
      try
      {
         BufferedOutputStream dest = null;
         FileInputStream fis = new FileInputStream(inFile);
         ZipInputStream zis = new ZipInputStream(new BufferedInputStream(fis));
         ZipEntry entry;
         while((entry = zis.getNextEntry()) != null) 
         {            
        	System.out.println("Extracting: " +entry);
            int count;
            byte data[] = new byte[BUFFER];
            // write the files to the disk
            FileOutputStream fos = new FileOutputStream(inFile.getParentFile()+"\\"+entry.getName());
            dest = new BufferedOutputStream(fos, BUFFER);
            while ((count = zis.read(data, 0, BUFFER))!= -1) 
               dest.write(data, 0, count);
            
            dest.flush();
            dest.close();
         }
         zis.close();
      } catch(Exception e) {e.printStackTrace();}
   }
}

