package org.biojava.bio.charts;


import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Paint;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;

public class BarChart extends  BaseChart
{
	 int imageWidth=600;
	 int imageHeight=400;
	 
	 
	public static void main(String[] args)
	{
		BarChart bar = new BarChart();
		String methods[] = {"BioProspector", "MEME", "MotifCut"};
		String series[] = {"nSn", "nPPV"};
		
		Double[][] values = new Double[series.length][methods.length];
		
		values[0][0] = 8d;
		values[0][1] = 4d;
		values[0][2] = 4d;
		values[1][0] = 5d;
		values[1][1] = 5d;
		values[1][2] = 6d;		
		
		try
		{
			bar.createChartImage("Evaluate Tompa", ".\\charts\\new.jpg", series, methods, values, false);
			bar.createChartImage("Evaluate Tompa", ".\\charts\\new22.jpg", false);
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	 public BarChart(){}
	 
	 public BarChart(int h, int w)
	 {		
		 imageHeight = h;
		 imageWidth = w;
	 }
	 

	// for testing
	public void createChartImage(String title, String path, boolean intValues)throws IOException 
	{
		 JFreeChart jfreechart= createBarChart(title, createDataset(), intValues);
		 BufferedImage img= jfreechart.createBufferedImage(imageWidth,imageHeight);
		 File file=new File(path);
		 saveJPEG(img,file,300);		 
	}
	 
	public void createChartImage(String title, String path, String series[],String x_values[],Double y_values[][], boolean intValues)throws IOException 
	{
		 JFreeChart jfreechart= createBarChart(title, series, x_values,y_values, intValues);
		 BufferedImage img= jfreechart.createBufferedImage(imageWidth, imageHeight);
		 File file=new File(path);
		 saveJPEG(img,file,300);		 
	}
	
	public JFreeChart createBarChart(String title, String series[],String x_values[],Double y_values[][], boolean intValues)
    {		
		CategoryDataset categorydataset = createCategoryDataset(series, x_values, y_values);
		
		return createBarChart(title, categorydataset, intValues);
             
    }
	
	public JFreeChart createBarChart(String title, CategoryDataset dataset, boolean intValues) 
	{	
		
        // create the chart...
        final JFreeChart chart = ChartFactory.createBarChart(
            title,         // chart title
            "",               // domain axis label
            "",                  // range axis label
            dataset,                  // data
            PlotOrientation.VERTICAL, // orientation
            true,                     // include legend
            true,                     // tooltips?
            false                     // URLs?
        );

        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...

        // set the background color for the chart...
        chart.setBackgroundPaint(Color.white);

        // get a reference to the plot for further customization...
        final CategoryPlot plot = chart.getCategoryPlot();
        plot.setBackgroundPaint(new Color(0.85f,0.85f,1f));
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);

        // set the range axis to display integers only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        if (intValues)
        	rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        else
        	rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());

        
        // disable bar outlines...
        final BarRenderer renderer = (BarRenderer) plot.getRenderer();
        //renderer.setDrawBarOutline(false);
        
  
        // see colors: http://mc2.cchem.berkeley.edu/Java/RGB/example1.html
        // blue
        GradientPaint gp0 = new GradientPaint(
            0.0f, 0.0f, new Color(0.3f,0.5f,0.8f), 
            0.0f, 0.0f, new Color(0.3f,0.5f,0.8f)
        );        
        // dark pink
        GradientPaint gp1 = new GradientPaint(
            0.0f, 0.0f, new Color(1f,0.5f,0.5f), 
            0.0f, 0.0f, new Color(1f,0.5f,0.5f)
        );
        // yellow
        GradientPaint gp2 = new GradientPaint(
                0.0f, 0.0f, new Color(1f,1f,0.4f), 
                0.0f, 0.0f, new Color(1f,1f,0.4f)
            );
        // violet
        Paint gp3 = new GradientPaint(
            0.0f, 0.0f, new Color(0.6f,0.5f,0.7f), 
            0.0f, 0.0f, new Color(0.6f,0.5f,0.7f)
        );
        // light pink
        GradientPaint gp4 = new GradientPaint(
                0.0f, 0.0f, new Color(1f,0.95f,0.95f), 
                0.0f, 0.0f, new Color(1f,0.95f,0.95f)
            );
        // Red
        GradientPaint gp5 = new GradientPaint(
                0.0f, 0.0f, new Color(1f,0.2f,0.2f), 
                0.0f, 0.0f, new Color(1f,0.2f,0.2f)
            );
        // Green
        GradientPaint gp6 = new GradientPaint(
                0.0f, 0.0f, new Color(0.5f,0.8f,0.3f), 
                0.0f, 0.0f, new Color(0.5f,0.8f,0.3f)
            );
        
        
        renderer.setSeriesPaint(0, gp0);
        renderer.setSeriesPaint(1, gp1);
        renderer.setSeriesPaint(2, gp2);
        renderer.setSeriesPaint(3, gp3);
        renderer.setSeriesPaint(4, gp4);
        renderer.setSeriesPaint(5, gp5);
        renderer.setSeriesPaint(6, gp6);
     
        renderer.setMaximumBarWidth(0.05d);
        renderer.setItemMargin(0d);			// remove space between items
        
        
        final CategoryAxis domainAxis = plot.getDomainAxis();
        domainAxis.setCategoryLabelPositions(
            CategoryLabelPositions.createUpRotationLabelPositions(Math.PI / 6.0)
        );
        // OPTIONAL CUSTOMISATION COMPLETED.
        
        return chart;
        
    }
	
	
	
	private CategoryDataset createDataset() 
	{
        
        // row keys...
        final String series1 = "nSn";
        final String series2 = "nPPV";
        final String series3 = "nPC";
        final String series4 = "nCC";
        final String series5 = "sSn";
        final String series6 = "sPPV";
        final String series7 = "sASP";

        // column keys...
        final String category1 = "Method 1";
        final String category2 = "Method 2";
        final String category3 = "Method 3";
        final String category4 = "Method 4";
        final String category5 = "Method 5";

        // create the dataset...
        final DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        dataset.addValue(0.1, series1, category1);
        dataset.addValue(0.8, series1, category2);
        dataset.addValue(0.3, series1, category3);
        dataset.addValue(0.522, series1, category4);
        dataset.addValue(0.75, series1, category5);

        dataset.addValue(0.5, series2, category1);
        dataset.addValue(0.2, series2, category2);
        dataset.addValue(0.9, series2, category3);
        dataset.addValue(0.6, series2, category4);
        dataset.addValue(0.4, series2, category5);

        dataset.addValue(0.9, series3, category1);
        dataset.addValue(0.5, series3, category2);
        dataset.addValue(0.1, series3, category3);
        dataset.addValue(0.2, series3, category4);
        dataset.addValue(0.65, series3, category5);
        
        dataset.addValue(0.4, series4, category1);
        dataset.addValue(0.3, series4, category2);
        dataset.addValue(0.6, series4, category3);
        dataset.addValue(0.8, series4, category4);
        dataset.addValue(0.65, series4, category5);
        
        dataset.addValue(0.4, series5, category1);
        dataset.addValue(0.9, series5, category2);
        dataset.addValue(0.2, series5, category3);
        dataset.addValue(0.3, series5, category4);
        dataset.addValue(0.65, series5, category5);
        
        dataset.addValue(0.4, series6, category1);
        dataset.addValue(0.7, series6, category2);
        dataset.addValue(0.8, series6, category3);
        dataset.addValue(0.2, series6, category4);
        dataset.addValue(0.65, series6, category5);
        
        dataset.addValue(0.6, series7, category1);
        dataset.addValue(0.3, series7, category2);
        dataset.addValue(0.2, series7, category3);
        dataset.addValue(0.7, series7, category4);
        dataset.addValue(0.65, series7, category5);
        
        dataset.addValue(0.4, series4, category1);
        dataset.addValue(0.5, series4, category2);
        dataset.addValue(0.5, series4, category3);
        dataset.addValue(0.3, series4, category4);
        dataset.addValue(0.65, series4, category5);
        
        return dataset;
        
    }
	
	
	private CategoryDataset createCategoryDataset(String series[],String x_values[],Double y_values[][])
    {
     	DefaultCategoryDataset defaultcategorydataset = new DefaultCategoryDataset();
     	
     	for(int i=0; i<series.length; i++)
     		for (int j=0; j< y_values[i].length; j++)
     		defaultcategorydataset.addValue(y_values[i][j], series[i], x_values[j]);
        return defaultcategorydataset;
    }	
	
	
}
