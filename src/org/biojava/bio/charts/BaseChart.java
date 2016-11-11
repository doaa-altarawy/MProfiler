package org.biojava.bio.charts;

import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGEncodeParam;
import com.sun.image.codec.jpeg.JPEGImageEncoder;

public abstract class BaseChart 
{
	
	 public static void saveJPEG(BufferedImage thumbImage, File file, double compression) throws IOException {
			BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(file));
			saveJPEG(thumbImage, out, compression);
			out.flush();
			out.close();
		}
	 public static void saveJPEG(BufferedImage thumbImage, OutputStream out, double compression) throws IOException {
			JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(out);
			JPEGEncodeParam param = encoder.getDefaultJPEGEncodeParam(thumbImage);
			compression = Math.max(0, Math.min(compression, 100));
			param.setQuality((float)compression / 100.0f, false);
			encoder.setJPEGEncodeParam(param);
			encoder.encode(thumbImage);
		}
}
