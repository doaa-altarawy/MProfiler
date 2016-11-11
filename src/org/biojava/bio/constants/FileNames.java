package org.biojava.bio.constants;

import java.io.File;

public class FileNames
{
	private static String pathSep = File.separator;

    private static String CHARTS_FOLDER = joinPathes(new String[]{"."+pathSep,"resources","charts"});
	private static String EXCEL_FOLDER = joinPathes(new String[]{"."+pathSep,"resources","excel"});
	public static final String ANSWER_FILE =  joinPathes(new String[]{"."+pathSep,"resources"})+ "answer.tompa";
	public static final String LENGHTH_FILE = joinPathes(new String[]{"."+pathSep,"resources"})+ "length.tompa";
	public static final String TOMPA_SEQUENCES = joinPathes(new String[]{"."+pathSep,"resources","Tompa"});

	private static String[] CHARTS = {"chart1_", "chart2_", "chart3_", "chart4_", "chart5_"};
	public static String LINE_CHART = CHARTS_FOLDER + "line.jpg";
	private static String statisticExcelFileName = "statistic";

	// private constructor, class can't be instantiated
	private FileNames(){}

	private static String joinPathes(String[] pathes) {
	    StringBuilder out = new StringBuilder();

	    for (String path : pathes) {
	        out.append(path + pathSep);
	    }

	    return out.toString();
	}

	public static String[] getCharts(Object n)
	{
		String[] charts = new String[CHARTS.length];

		for (int i=0; i<CHARTS.length; i++)
			charts[i] = CHARTS_FOLDER+CHARTS[i]+n+".jpg";
		return charts;
	}

	public static String getStatisticsExcelFileName(Object n)
	{
		return EXCEL_FOLDER + statisticExcelFileName+n+".xls";
	}
}
