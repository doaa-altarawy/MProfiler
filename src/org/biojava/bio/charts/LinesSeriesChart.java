package org.biojava.bio.charts;


import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.LineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.VerticalAlignment;

public class LinesSeriesChart extends  BaseChart
{ 
	 int imageWidth=600;
	 int imageHeight=400;	 
	 
	 boolean showLegend=true;
	 boolean showLines = true;
	 boolean showShapes = true;	 
	 

	Color[] colors;
	 
	 public static void main(String[] args)
		{
		 	LinesSeriesChart chart = new LinesSeriesChart();
			String methods[] = {"1", "10", "15"};
			String eval[] = {"nSn", "nPPV"};
			
			Double[][] values = new Double[eval.length][methods.length];
			
			values[0][0] = 8.5d;
			values[0][1] = 4.2d;
			values[0][2] = 4.89d;
			values[1][0] = 5.3d;
			values[1][1] = 5d;
			values[1][2] = 6d;		
			
			try
			{
				chart.createCategoryChartImage("My title", "xxx", "nSn", ".\\charts\\new.jpg", eval, methods, values);
				//chart.createXYChartImage("My title","yyy",".\\charts\\newxy.jpg",chart.createDataset());
			} 
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
	 	 
	 public LinesSeriesChart()
	 {
		colors = new Color[8];
		colors[0] = Color.BLUE;
		colors[1] = Color.RED ;
		colors[2] = Color.YELLOW ;
		colors[3] = Color.GREEN ;
		colors[4] = Color.GRAY ;
		colors[5] = Color.CYAN ;
		colors[6] = Color.PINK ;
		colors[7] = Color.ORANGE ;		
	 }
			 
	 public LinesSeriesChart(int h, int w)
	 {		
		 imageHeight = h;
		 imageWidth = w;
	 }
	
	public void createXYChartImage(String title, String xAxis, String yAxis, String path, String series[],Double x_values[],Double y_values[][]) throws IOException 
	{
		 JFreeChart jfreechart= createXYChart(title, xAxis, yAxis, series, x_values, y_values);
		 BufferedImage img= jfreechart.createBufferedImage(imageWidth,imageHeight);
		 File file=new File(path);
		 saveJPEG(img,file,300);		
	}

	
	public void createCategoryChartImage(String title, String xAxis, String yAxis, String path, String series[],String x_values[],Double y_values[][])throws IOException 
	{
		 JFreeChart jfreechart= createCategoryChart(title, xAxis, yAxis, series, x_values,y_values);
		 BufferedImage img= jfreechart.createBufferedImage(imageWidth,imageHeight);
		 File file=new File(path);
		 saveJPEG(img,file,300);		
	}
	
	public JFreeChart createCategoryChart(String title, String xAxis, String yAxis, String series[],String x_values[],Double y_values[][])
    {
		
		CategoryDataset categorydataset = createCategoryDataset(series, x_values, y_values);
		
		// create the chart...
        final JFreeChart chart = ChartFactory.createLineChart(
            title,       // chart title
            xAxis,                    // domain axis label
            yAxis,                   // range axis label
            categorydataset,           // data
            PlotOrientation.VERTICAL,  // orientation
            showLegend,                // include legend
            true,                      // tooltips
            false                      // urls
        );

        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
        chart.setBackgroundPaint(Color.white);
        
        
        CategoryPlot categoryplot = chart.getCategoryPlot();	
	    categoryplot.setBackgroundPaint(new Color(0.85f,0.85f,1f));
	    categoryplot.setRangeGridlinePaint(Color.white);
	    
	    LegendTitle legend = chart.getLegend();
	    if(legend!=null)
	    {
	       legend.setPosition(RectangleEdge.RIGHT);
	       legend.setVerticalAlignment(VerticalAlignment.TOP);
	    }       


        // customise the range axis...
        final NumberAxis rangeAxis = (NumberAxis) categoryplot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
        rangeAxis.setAutoRangeIncludesZero(true);
        
        // customise the renderer...
        final LineAndShapeRenderer renderer = (LineAndShapeRenderer)categoryplot.getRenderer();
        renderer.setShapesVisible(showShapes);
        renderer.setLinesVisible(showLines);
        
        // OPTIONAL CUSTOMISATION COMPLETED.
        
        return chart;
    }
	
	
	public JFreeChart createXYChart(String title, String xAxis, String yAxis, String series[],Double x_values[],Double y_values[][]) 
	{
		XYDataset dataset = createXYDataset(series, x_values, y_values);
		return createXYChart(title, xAxis, yAxis, dataset);
	}
	
	 private JFreeChart createXYChart(String title, String xAxis, String yAxis, XYDataset dataset) 
	 {
		   
		 // create the chart...
		final JFreeChart chart = ChartFactory.createXYLineChart(
		    title,      // chart title
		    xAxis,                    // x axis label
		    yAxis,                    // y axis label
		    dataset,                  // data
		    PlotOrientation.VERTICAL,
		    showLegend,               // include legend
		    true,                     // tooltips
		    false                     // urls
		);
		
		// NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
		chart.setBackgroundPaint(Color.white);
		
		// get a reference to the plot for further customization...
		final XYPlot plot = chart.getXYPlot();
		plot.setBackgroundPaint(new Color(0.85f,0.85f,1f));
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);
		
		final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setShapesVisible(showShapes);
		renderer.setLinesVisible(showLines);
		plot.setRenderer(renderer);
		
		LegendTitle legend = chart.getLegend();
		if(legend!=null)
		{
		   legend.setPosition(RectangleEdge.RIGHT);
		   legend.setVerticalAlignment(VerticalAlignment.TOP);
		}
		// change the auto tick unit selection to integer units only...
		final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		rangeAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
		// OPTIONAL CUSTOMISATION COMPLETED.
		        
		return chart;
	        
	 }
	
	
	private  CategoryDataset createCategoryDataset(String series[], String x_values[],Double y_values[][])
    {
     	DefaultCategoryDataset defaultcategorydataset = new DefaultCategoryDataset();
     	
     	for(int i=0; i<series.length; i++)
     		for (int j=0; j< y_values[i].length; j++)
     			defaultcategorydataset.addValue(y_values[i][j], series[i], x_values[j]);
        
     	return defaultcategorydataset;
    }
	
	private XYDataset createXYDataset(String series[], Double x_values[],Double y_values[][])
    {
		XYSeriesCollection xydataset = new XYSeriesCollection();
     	
     	for(int i=0; i<series.length; i++)
     	{
     		XYSeries seriesi = new XYSeries(series[i]);
     		for (int j=0; j< y_values[i].length; j++)
     		{
     			seriesi.add(x_values[j], y_values[i][j]);
     		}
     		xydataset.addSeries(seriesi);
     	}
     	
        return xydataset;
    }
	

	public void setShowLegend(boolean showLegend)
	{
		this.showLegend = showLegend;
	}

	public void setShowLines(boolean showLines)
	{
		this.showLines = showLines;
	}

	public void setShowShapes(boolean showShapes)
	{
		this.showShapes = showShapes;
	}
	
}
