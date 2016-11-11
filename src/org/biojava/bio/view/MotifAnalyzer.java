
package org.biojava.bio.view;
import com.cloudgarden.layout.AnchorLayout;
import com.zfqjava.swing.JTableReadTableModelTask;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.ComponentOrientation;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.GroupLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.text.MaskFormatter;

import org.biojava.bio.BioException;
import org.biojava.bio.constants.Dataset;
import org.biojava.bio.constants.DatasetType;
import org.biojava.bio.constants.MotifFinder;
import org.biojava.bio.constants.OrganismCode;
import org.biojava.bio.motif.Motif;
import org.biojava.bio.motif.MotifFinders.Method;
import org.biojava.bio.motif.Motif.Site;
import org.biojava.bio.motif.gibbs.SimpleGibbsAlignerDemo;

import org.biojava.bio.charts.LinesSeriesChart;
import org.biojava.bio.gui.sequence.FeatureBlockSequenceRenderer;
import org.biojava.bio.gui.sequence.FilteringRenderer;
import org.biojava.bio.gui.sequence.MultiLineRenderer;
import org.biojava.bio.gui.sequence.RectangularBeadRenderer;
import org.biojava.bio.gui.sequence.RulerRenderer;
import org.biojava.bio.gui.sequence.SequencePanel;
import org.biojava.bio.gui.sequence.SymbolSequenceRenderer;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.OptimizableFilter;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.tools.FileConverterTools;
import org.biojava.bio.tools.Filter;
import org.biojava.bio.tools.MotifTools;
import org.biojava.bio.tools.TompaDataset;
import org.biojava.bio.tools.UnZip;
import org.biojava.bio.constants.FileNames;
import org.biojava.utils.ChangeVetoException;
import org.jfree.chart.ChartPanel;



import javax.swing.JTextArea;
import javax.swing.JTextPane;
import javax.swing.LayoutStyle;

/**
* This code was edited or generated using CloudGarden's Jigloo
* SWT/Swing GUI Builder, which is free for non-commercial
* use. If Jigloo is being used commercially (ie, by a corporation,
* company or business for any purpose whatever) then you
* should purchase a license for each developer using Jigloo.
* Please visit www.cloudgarden.com for details.
* Use of Jigloo implies acceptance of these licensing terms.
* A COMMERCIAL LICENSE HAS NOT BEEN PURCHASED FOR
* THIS MACHINE, SO JIGLOO OR THIS CODE CANNOT BE USED
* LEGALLY FOR ANY CORPORATE OR COMMERCIAL PURPOSE.
*/

public class MotifAnalyzer implements ActionListener
{

	{
		//Set Look & Feel
		try {
			javax.swing.UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
		} catch(Exception e) {
			e.printStackTrace();
		}
	}

	private static Logger log = Logger.getLogger(MotifAnalyzer.class.getName());
    
	private static final int INITIAL_SCALE = 25;

    private JFrame      mainFrame;
    private JPanel		inputTab;
    private JPanel		seqViewerTab;
    private JPanel		motifScoreTab;
    private JTextPane outputPane;
    private JTabbedPane tabs;
    private JScrollPane jScrollPane2;
    private JScrollPane jScrollPane1;
    private JTextArea inputTextArea;
    private static String convertMotifVoterMenuName = "Convert MotifVoter Web Option2 Result";
    private static String runPatchConverMVMenuName = "Run Patch Conversion of MV Option2(multi Dir)";
    private static String convertAllMenuName = "Convert All Motif Finders Outputs";
    private static String runPatchConverAllMenuName = "Run Patch Conversion of All Motif Finder Outputs(multi Dir)";
    private static String copyCleanInputFilesName = "Copy Clean input Files";
    private static String loadMenuName = "load Analysis";
    private static String saveMenuName = "Save last Analysis";
    private static String gibbsMenuName = "Gibbs Aligner";
    private static String compareSetsName = "Compare Sets";
    private static String mergeCompareSetsName = "Merge CompareSets";
    private static String removeRedundantSitesName = "Remove Redundant Sites";
    private static String extractSitesName = "Extract Sites";
    private static String drawKDistanceGraphName = "Draw K-Distance Graph";
    private static String computeWeightName = "Compute Weight";
    private static String infoContentMenuName = "Calculate Score";
    private static String runBioProspectorName = "BioProspector";
    private static String runKMeansName = "K-Means"; 
    private static String runDBScanName = "DBScan";
    private static String runOPTICSName = "OPTICS";
    private static String runMotifVoterName = "MotifVoter"; 
    private static String runCliqueName = "Clique"; 
    private static String evaluateTompaMotifEachName = "Evaluate Tompa Motifs";
    private static String evaluateTompaMotifAccName = "Evaluate Tompa Motif Accumulative";
    private static String evaluateTompaMotifTotalName = "Evaluate Tompa Motif Total";
    private static String evaluateMotifVoterName = "Evaluate each MotifVoter Motif in a File";
    private static String evaluateMotifVoterTotalName = "Evaluate Total Merged MV Motifs in a File";
    private static String evaluateMVDirFirstEleName = "Evaluate MV Motifs Dir (First Motif in the file)";
    private static String evaluateMVDirMergedName = "Evaluate MV Motifs Dir (Merged file Motifs)";
    private static String compareMotifVoterFirstName = "Compare two MV Motifs (First Motif in the file)";
    private static String compareMotifVoterAllName = "Compare two MV Motifs (Merge Motifs in each file)";
    private static String evaluateAllResultsName = "Evaluate All Results (dm, hm, mus, yst)";
    
    Font tabFont = new Font(Font.DIALOG, Font.BOLD, 16);
    Font menuFont = new Font(Font.DIALOG, Font.BOLD, 14);
    private JComboBox numOfClusters;
    private JLabel jLabel9;
    private JComboBox minPoints;
    private JLabel jLabel8;
    private JComboBox clusteringEpsilon;
    private JLabel jLabel7;
    private JComboBox datasetTypeComboBox;
    private JLabel jLabel6;
    private JComboBox datasetComboBox;
    private JLabel jLabel5;
    private JTextField motifWidth;
    private JLabel jLabel4;
    private JTextField numOfTrials;
    private JLabel jLabel3;
    private JComboBox organismType;
    private JLabel jLabel2;
    private JCheckBox checkReverse;
    private JTextField numOfMotifs;
    private JCheckBox multipleDataSets;
    private JLabel jLabel1;
             
	private JMenuBar jJMenuBar;

	EvaluationController controller;
	
	public MotifAnalyzer()
    {
		Logger.getLogger("").setLevel(Level.FINE);
		initComponents();
		controller = new EvaluationController();
    }

	
	public static void main(String[] args)
    {        
        new MotifAnalyzer();
    }
	
    private void initComponents()
    {
    	mainFrame = new JFrame();
    	mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	mainFrame.setBackground(Color.white);
    	mainFrame.setTitle("Motif Analayzer");
    	mainFrame.setJMenuBar(getJJMenuBar());
    	
    	tabs  = new JTabbedPane();    	
		tabs.setFont(tabFont);
    	mainFrame.add(tabs, BorderLayout.CENTER);
    	
    	    	
    	tabs.addTab("Input", getInputTab());
    	
    	seqViewerTab = new JPanel();
    	seqViewerTab.setLayout(new BorderLayout());
    	 
    	// has a problem when view change in horizontal scrolling///////////////
    	JScrollPane s = new JScrollPane(seqViewerTab);//, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED, ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
    	
    	tabs.addTab("Sequence Viewer", s);
    	tabs.addTab("Motif Scores", getMotifScoreTab());
    	
	
		mainFrame.setSize(800, 600);
    	mainFrame.setExtendedState(mainFrame.getExtendedState()|JFrame.MAXIMIZED_BOTH);
    	
		mainFrame.setVisible(true);
    }

    private JPanel getInputTab()
    {
    	inputTab = new JPanel(new GridLayout(1,1)); //////////////////
    	{
    		JPanel generalPar = new JPanel();
    		inputTab.add(generalPar);
    		GroupLayout inputTabLayout = new GroupLayout((JComponent)generalPar);
    		generalPar.setLayout(inputTabLayout);    	
    		generalPar.setBorder(BorderFactory.createTitledBorder(BorderFactory.createLineBorder(Color.black), "General Parameters", 0, 0, new Font(Font.DIALOG, Font.BOLD,16), Color.black));
    		generalPar.setPreferredSize(new java.awt.Dimension(787, 291));
    		inputTabLayout.setHorizontalGroup(inputTabLayout.createSequentialGroup()
    			.addContainerGap(12, 12)
    			.addGroup(inputTabLayout.createParallelGroup()
    			    .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			        .addComponent(getCheckReverse(), GroupLayout.PREFERRED_SIZE, 198, GroupLayout.PREFERRED_SIZE)
    			        .addGap(250))
    			    .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			        .addComponent(getMultipleDataSets(), GroupLayout.PREFERRED_SIZE, 193, GroupLayout.PREFERRED_SIZE)
    			        .addGap(255))
    			    .addGroup(inputTabLayout.createSequentialGroup()
    			        .addGroup(inputTabLayout.createParallelGroup()
    			            .addComponent(getJLabel1(), GroupLayout.Alignment.LEADING, GroupLayout.PREFERRED_SIZE, 187, GroupLayout.PREFERRED_SIZE)
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel8(), GroupLayout.PREFERRED_SIZE, 113, GroupLayout.PREFERRED_SIZE)
    			                .addGap(74))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel7(), GroupLayout.PREFERRED_SIZE, 80, GroupLayout.PREFERRED_SIZE)
    			                .addGap(107))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel2(), GroupLayout.PREFERRED_SIZE, 128, GroupLayout.PREFERRED_SIZE)
    			                .addGap(59))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel3(), GroupLayout.PREFERRED_SIZE, 148, GroupLayout.PREFERRED_SIZE)
    			                .addGap(39))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel4(), GroupLayout.PREFERRED_SIZE, 137, GroupLayout.PREFERRED_SIZE)
    			                .addGap(50))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel6(), GroupLayout.PREFERRED_SIZE, 126, GroupLayout.PREFERRED_SIZE)
    			                .addGap(61))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel5(), GroupLayout.PREFERRED_SIZE, 108, GroupLayout.PREFERRED_SIZE)
    			                .addGap(79))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getJLabel9(), GroupLayout.PREFERRED_SIZE, 155, GroupLayout.PREFERRED_SIZE)
    			                .addGap(32)))
    			        .addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
    			        .addGroup(inputTabLayout.createParallelGroup()
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getMinPoints(), GroupLayout.PREFERRED_SIZE, 122, GroupLayout.PREFERRED_SIZE)
    			                .addGap(134))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getClusteringEpsilon(), GroupLayout.PREFERRED_SIZE, 122, GroupLayout.PREFERRED_SIZE)
    			                .addGap(134))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getNumOfMotifs(), GroupLayout.PREFERRED_SIZE, 69, GroupLayout.PREFERRED_SIZE)
    			                .addGap(187))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getNumOfTrials(), GroupLayout.PREFERRED_SIZE, 69, GroupLayout.PREFERRED_SIZE)
    			                .addGap(187))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getMotifWidth(), GroupLayout.PREFERRED_SIZE, 69, GroupLayout.PREFERRED_SIZE)
    			                .addGap(187))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getDatasetTypeComboBox(), GroupLayout.PREFERRED_SIZE, 145, GroupLayout.PREFERRED_SIZE)
    			                .addGap(111))
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getDatasetComboBox(), GroupLayout.PREFERRED_SIZE, 145, GroupLayout.PREFERRED_SIZE)
    			                .addGap(111))
    			            .addComponent(getSpeciesType(), GroupLayout.Alignment.LEADING, GroupLayout.PREFERRED_SIZE, 256, GroupLayout.PREFERRED_SIZE)
    			            .addGroup(GroupLayout.Alignment.LEADING, inputTabLayout.createSequentialGroup()
    			                .addComponent(getNumOfClusters(), GroupLayout.PREFERRED_SIZE, 122, GroupLayout.PREFERRED_SIZE)
    			                .addGap(134)))))
    			.addContainerGap(317, 317));
    		inputTabLayout.setVerticalGroup(inputTabLayout.createSequentialGroup()
    			.addContainerGap()
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getNumOfMotifs(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel1(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, 20, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getSpeciesType(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel2(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, 17, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getNumOfTrials(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, 24, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel3(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getMotifWidth(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel4(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, 16, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getDatasetTypeComboBox(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel6(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getDatasetComboBox(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel5(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE))
    			.addGap(20)
    			.addComponent(getMultipleDataSets(), GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
    			.addGap(19)
    			.addComponent(getCheckReverse(), GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
    			.addGap(16)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getClusteringEpsilon(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel7(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getMinPoints(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel8(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE))
    			.addPreferredGap(LayoutStyle.ComponentPlacement.UNRELATED)
    			.addGroup(inputTabLayout.createParallelGroup(GroupLayout.Alignment.BASELINE)
    			    .addComponent(getNumOfClusters(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, 22, GroupLayout.PREFERRED_SIZE)
    			    .addComponent(getJLabel9(), GroupLayout.Alignment.BASELINE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE))
    			.addContainerGap(65, 65));
    	}
    	
    	return inputTab;
    }

    public void showMotifs(SequenceIterator seqItr, File motifFile)
    {
    	try{
    		List<Motif> motifs = FileConverterTools.readMotifVoter(motifFile, false);
    	  	showMotifs(seqItr, motifs);
    	}
    	catch (Exception e) {
			e.printStackTrace();
		}
    }
    /**
     * 
     * @param f1 sequence file
     * @param motifs list of motifs to display
     */
    public void showMotifs(SequenceIterator seqItr, List<Motif> motifs)
    {
    	seqViewerTab.removeAll();
		seqViewerTab.setLayout(new BorderLayout());
		
    	try
        {    		
    		Vector<Sequence> seqs = new Vector<Sequence>();
    		Set<MotifFinder> finders = new HashSet<MotifFinder>(); 
    		
    		while(seqItr.hasNext())
    			seqs.add(seqItr.nextSequence());  		
    		    		
    		// labels panel
    		
    		
    		JPanel labelPanel = new JPanel();    		
    		seqViewerTab.add(labelPanel, BorderLayout.NORTH);
    		
    		
    		// Sequence panels array
    		JPanel seqViewer = new JPanel(new GridLayout(seqs.size(), 1));
    		seqViewerTab.add(seqViewer, BorderLayout.CENTER);
    		
    		SequencePanel[] seqPanel = new SequencePanel[seqs.size()];
            GridLayout seqPanelLayout = new GridLayout(2, 1);
            seqPanelLayout.setColumns(1);
            seqPanelLayout.setHgap(5);
            seqPanelLayout.setVgap(5);
            seqPanelLayout.setRows(2);
            int i=0;
            
            for (Iterator<Sequence> itr=seqs.iterator(); itr.hasNext(); i++)
            {
            	Sequence sourceSeq = itr.next();
            	String name = sourceSeq.getName().substring(sourceSeq.getName().indexOf("_")+1);
            	seqPanel[i] = new SequencePanel();
            	seqPanel[i].setLayout(seqPanelLayout);
            	int num = 0;
            	String[] f = new String[motifs.size()];
            	for (Iterator<Motif> m=motifs.iterator();m.hasNext();)
	            {
		            Motif motif = m.next();		            
		            if (motif.getSites().get(name)==null) continue;
	            	for (Iterator<Site> itr2=motif.getSites().get(name).iterator(); itr2.hasNext();)
		            {
	            		Site seq = itr2.next();				        
						finders.add(motif.getFinder());
						seq.location = new RangeLocation(sourceSeq.length()+seq.location.getMin(), sourceSeq.length()+seq.location.getMax());	
						seq.type = motif.getFinder().name()+"_"+num;
						sourceSeq.createFeature(seq);			        
		            }	            	
	            	f[num] = motif.getFinder().name()+"_"+num;
	            	num++;
	            }
	            
		        seqPanel[i].setSequence(sourceSeq);
		        seqPanel[i].setRange(new RangeLocation(1, sourceSeq.length()));		       
				// Magic number from EmblViewer
		        seqPanel[i].setScale(Math.exp(-INITIAL_SCALE / 7.0) * 20.0);
		        seqPanel[i].setDirection(SequencePanel.HORIZONTAL);
		     		
				OptimizableFilter[]  optFilter = new OptimizableFilter[num];//MotifFinder.values().length];
				FeatureBlockSequenceRenderer[] renderer = new FeatureBlockSequenceRenderer[num];//MotifFinder.values().length];
				
				for (int k=0; k<num; k++)//MotifFinder.values().length; k++)
				{
					optFilter[k] = new FeatureFilter.ByType(f[k]);//MotifFinder.values()[k].name());
					renderer[k] = new FeatureBlockSequenceRenderer(new RectangularBeadRenderer(10.0f, 5.0f, Color.black, MotifFinder.valueOf(f[k].substring(0, f[k].indexOf('_'))).color, new BasicStroke()));
				}				   
			
		        MultiLineRenderer multi = new MultiLineRenderer();
		       
		        for (int k=0; k<renderer.length; k++)
		        	multi.addRenderer(new FilteringRenderer(renderer[k], optFilter[k], false));
								
		        multi.addRenderer(new SymbolSequenceRenderer());
		        multi.addRenderer(new RulerRenderer());
		       		        
		        seqPanel[i].setRenderer(multi);
		        seqPanel[i].setEnabled(false);
		        
		        JScrollPane seqScroll = new JScrollPane(seqPanel[i]);
		        
		        JPanel temp = new JPanel(new BorderLayout());
		        
		        temp.add(getControlBox(seqPanel[i], sourceSeq.getName()), BorderLayout.NORTH);
		        temp.add(seqScroll, BorderLayout.CENTER);
		        
		        seqViewer.add(temp);
            }		        
	        
            for (Iterator<MotifFinder> itr = finders.iterator(); itr.hasNext();)
            {
            	MotifFinder finder = itr.next();
            	// Add label of this MotifFinder
		        labelPanel.add(new JLabel(finder.name()));
	    		JLabel l = new JLabel("      ");
	    		l.setOpaque(true);
	    		l.setBackground(finder.color);
	    		l.setBorder(BorderFactory.createLineBorder(Color.black));
	    		labelPanel.add(l);
            }
            mainFrame.setVisible(true);		        
                       
        }
        catch (Exception e)
        {
            e.printStackTrace();           
        }
    }
    private class SliderListener implements ChangeListener
    {
    	SequencePanel seqPanel;
    	
    	public SliderListener(SequencePanel seqPanel)
    	{
    		this.seqPanel = seqPanel;
    	}
    	
    	public void stateChanged(ChangeEvent ce)
		{
		    JSlider source = (JSlider) ce.getSource();
		    if (! source.getValueIsAdjusting())
		    {
				// val is between 1 and 50
				int val = source.getValue();
				double s = Math.exp(-val / 7.0) * 20.0;
		
				seqPanel.setScale(s);
		    }
		}
    }

	/**
	 * This method initializes jJMenuBar	
	 * 	
	 * @return javax.swing.JMenuBar	
	 */
	private JMenuBar getJJMenuBar() {
		if (jJMenuBar == null) {
			jJMenuBar = new JMenuBar();
			jJMenuBar.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);			
		}
		/* File Menu */
		JMenu fileMenu = new JMenu("File");
		fileMenu.setFont(menuFont);
		jJMenuBar.add(fileMenu);
		
		JMenuItem loadItem = new JMenuItem(loadMenuName);
		loadItem.addActionListener(this);
		fileMenu.add(loadItem);
		
		JMenuItem saveItem = new JMenuItem(saveMenuName);
		saveItem.addActionListener(this);
		fileMenu.add(saveItem);
		
		JMenuItem openMotifVoterItem = new JMenuItem(convertMotifVoterMenuName);
		openMotifVoterItem.addActionListener(this);
		fileMenu.add(openMotifVoterItem);
		
		JMenuItem runPatchConverMVItem = new JMenuItem(runPatchConverMVMenuName);
		runPatchConverMVItem.addActionListener(this);
		fileMenu.add(runPatchConverMVItem);
		
		JMenuItem convertAllItem = new JMenuItem(convertAllMenuName);
		convertAllItem.addActionListener(this);
		fileMenu.add(convertAllItem);
		
		
		JMenuItem runPatchConverAllItem = new JMenuItem(runPatchConverAllMenuName);
		runPatchConverAllItem.addActionListener(this);
		fileMenu.add(runPatchConverAllItem);
					
		
		
		JMenuItem compareSetsItem = new JMenuItem(compareSetsName);
		compareSetsItem.addActionListener(this);
		fileMenu.add(compareSetsItem);
		
		JMenuItem mergeCompareSetsItem = new JMenuItem(mergeCompareSetsName);
		mergeCompareSetsItem.addActionListener(this);
		fileMenu.add(mergeCompareSetsItem);
		
		JMenuItem copyFilesWithExtensionItem = new JMenuItem(copyCleanInputFilesName);
		copyFilesWithExtensionItem.addActionListener(this);
		fileMenu.add(copyFilesWithExtensionItem);
				
		JMenuItem exitItem = new JMenuItem("Exit");
		exitItem.addActionListener(new ActionListener(){
		      public void actionPerformed(ActionEvent actionEvent) {System.exit(0);}});
		fileMenu.add(exitItem);
		
		/* Tools Menu */
		JMenu toolsMenu = new JMenu("Tools");
		toolsMenu.setFont(menuFont);
		jJMenuBar.add(toolsMenu);
		
		JMenuItem gibbsItem = new JMenuItem(gibbsMenuName);
		gibbsItem.addActionListener(this);
		toolsMenu.add(gibbsItem);
		
		JMenuItem infoContentItem = new JMenuItem(infoContentMenuName);
		infoContentItem.addActionListener(this);
		toolsMenu.add(infoContentItem);
		
		JMenuItem removeRedundantSitesItem = new JMenuItem(removeRedundantSitesName);
		removeRedundantSitesItem.addActionListener(this);
		toolsMenu.add(removeRedundantSitesItem);
		
		JMenuItem extractSitesItem = new JMenuItem(extractSitesName);
		extractSitesItem.addActionListener(this);
		toolsMenu.add(extractSitesItem);
		
		JMenuItem drawKDistanceGraphItem = new JMenuItem(drawKDistanceGraphName);
		drawKDistanceGraphItem.addActionListener(this);
		toolsMenu.add(drawKDistanceGraphItem);
		
		JMenuItem computeWeightItem = new JMenuItem(computeWeightName);
		computeWeightItem.addActionListener(this);
		toolsMenu.add(computeWeightItem);
		
		/* Run Methods Menu */
		JMenu methodMenu = new JMenu("Run Method");
		methodMenu.setFont(menuFont);
		jJMenuBar.add(methodMenu);
		
		/*
		JMenuItem bioProspectorItem = new JMenuItem(runBioProspectorName);
		bioProspectorItem.addActionListener(this);
		methodMenu.add(bioProspectorItem);
		
		JMenuItem KMeansItem = new JMenuItem(runKMeansName);
		KMeansItem.addActionListener(this);
		methodMenu.add(KMeansItem);
		
		JMenuItem DBScanItem = new JMenuItem(runDBScanName);
		DBScanItem.addActionListener(this);
		methodMenu.add(DBScanItem);
		
		JMenuItem OPTICSItem = new JMenuItem(runOPTICSName);
		OPTICSItem.addActionListener(this);
		methodMenu.add(OPTICSItem);
		
		JMenuItem motifVoteItem = new JMenuItem(runMotifVoterName);
		motifVoteItem.addActionListener(this);
		methodMenu.add(motifVoteItem);
		
		JMenuItem runCliqueItem = new JMenuItem(runCliqueName);
		runCliqueItem.addActionListener(this);
		methodMenu.add(runCliqueItem);
		*/
		
		for (int i=0; i<Method.values().length; i++)
		{
			JMenuItem item = new JMenuItem(Method.values()[i].name());
			item.addActionListener(this);
			methodMenu.add(item);
		}
		
		/* Run Analysis for Methods Menu */
		JMenu analysisMenu = new JMenu("Analysis");
		analysisMenu.setFont(menuFont);
		jJMenuBar.add(analysisMenu);
		
		JMenuItem evaluateTompaMotifEachItem = new JMenuItem(evaluateTompaMotifEachName);
		evaluateTompaMotifEachItem.addActionListener(this);
		analysisMenu.add(evaluateTompaMotifEachItem);
		
		JMenuItem evaluateTompaMotifItem = new JMenuItem(evaluateTompaMotifAccName);
		evaluateTompaMotifItem.addActionListener(this);
		analysisMenu.add(evaluateTompaMotifItem);
		
		JMenuItem evaluateTompaMotifTotalItem = new JMenuItem(evaluateTompaMotifTotalName);
		evaluateTompaMotifTotalItem.addActionListener(this);
		analysisMenu.add(evaluateTompaMotifTotalItem);
		
		JMenuItem evaluateMotifVoterItem = new JMenuItem(evaluateMotifVoterName);
		evaluateMotifVoterItem.addActionListener(this);
		analysisMenu.add(evaluateMotifVoterItem);
		
		JMenuItem evaluateMotifVoterTotalItem = new JMenuItem(evaluateMotifVoterTotalName);
		evaluateMotifVoterTotalItem.addActionListener(this);
		analysisMenu.add(evaluateMotifVoterTotalItem);
		
		JMenuItem evaluateMVDirFirstEleItem = new JMenuItem(evaluateMVDirFirstEleName);
		evaluateMVDirFirstEleItem.addActionListener(this);
		analysisMenu.add(evaluateMVDirFirstEleItem);
		
		JMenuItem evaluateMVDirMergedItem = new JMenuItem(evaluateMVDirMergedName);
		evaluateMVDirMergedItem.addActionListener(this);
		analysisMenu.add(evaluateMVDirMergedItem);
		
		
		JMenuItem compareMotifVoterFirstItem = new JMenuItem(compareMotifVoterFirstName);
		compareMotifVoterFirstItem.addActionListener(this);
		analysisMenu.add(compareMotifVoterFirstItem);
		
		JMenuItem compareMotifVoterAllItem = new JMenuItem(compareMotifVoterAllName);
		compareMotifVoterAllItem.addActionListener(this);
		analysisMenu.add(compareMotifVoterAllItem);
		
		JMenuItem EvaluateAllResultsItem = new JMenuItem(evaluateAllResultsName);
		EvaluateAllResultsItem.addActionListener(this);
		analysisMenu.add(EvaluateAllResultsItem);
		
		return jJMenuBar;
	}
	
		public void actionPerformed(ActionEvent e) 
		{
			final JFileChooser readVoter = new JFileChooser();			
			readVoter.setFileFilter(new Filter("voter"));
			readVoter.setDialogTitle("Select Motifs file (in MotifVoter format)");
			
			final JFileChooser reader = new JFileChooser();			
			reader.setDialogTitle("Select Motifs file");
			
			final JFileChooser readFASTA = new JFileChooser();
			readFASTA.setFileFilter(new Filter("fasta"));
			readFASTA.setDialogTitle("Select Sequence File (in FASTA format)");
			
			final JFileChooser readTompa = new JFileChooser();
			readTompa.setFileFilter(new Filter("tompa"));
			readTompa.setDialogTitle("Select Motifs File (in Tompa's fromat)");
			
			final JFileChooser readDirectory = new JFileChooser();
			readDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			readDirectory.setDialogTitle("Select Directory");
			
			final Dataset dataset = (Dataset) datasetComboBox.getSelectedItem();
			final DatasetType datasetType = (DatasetType) datasetTypeComboBox.getSelectedItem();
			final boolean readReverse = checkReverse.isSelected();
			final OrganismCode org = OrganismCode.getOrganism(organismType.getSelectedIndex());
			final int minPoints = (Integer)this.minPoints.getSelectedItem();
			final double epsilon  = Double.parseDouble((String)this.clusteringEpsilon.getSelectedItem());
			final int numOfClusters = (Integer)this.numOfClusters.getSelectedItem();
			final boolean multipleDataSets = this.multipleDataSets.isSelected();
			try
			{
			if (((JMenuItem)e.getSource()).getText().equals(convertMotifVoterMenuName))
			{					       	
				if (reader.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION)
				{
					String fileName = reader.getSelectedFile().getCanonicalPath();
					if (reader.getSelectedFile().getName().contains(".zip")) // if file is zipped, extract it
					{
						fileName = UnZip.unZip(reader.getSelectedFile(), ".txt", ".MVW2");
						(new File(fileName)).deleteOnExit();
					}
					
					if (fileName!=null)
						FileConverterTools.readMotifVoterWeb(new File(fileName), readReverse, true);
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(runPatchConverMVMenuName))
			{				
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{
					File[] directories = readDirectory.getSelectedFile().listFiles();	
					File out2 = new File(readDirectory.getSelectedFile().getParentFile()+"\\MVW2");
					File out = new File(readDirectory.getSelectedFile().getParentFile()+"\\MVW");
					out.mkdir();
					out2.mkdir();
					for (int i=0; i<directories.length; i++)
					{
						File[] files2 = directories[i].listFiles(new Filter("zip"));
						File[] files = directories[i].listFiles(new Filter("MVW.voter"));
						if (files2.length>0)
						{
							String fileName = UnZip.unZip(files2[0], ".txt", ".MVW2");
							if (fileName!=null)
							{
								fileName = FileConverterTools.readMotifVoterWeb(new File(fileName), readReverse, true);
								FileConverterTools.copy(new File(fileName), out2);
							}
						}
						if (files.length>0)
						{
							FileConverterTools.copy(files[0], out);							
						}
					}						
				}
			}	
			else if (((JMenuItem)e.getSource()).getText().equals(convertAllMenuName))
			{				
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{
					controller.convertFilesAndMerge(readDirectory.getSelectedFile(), dataset, readReverse);
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(runPatchConverAllMenuName))
			{				
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{
					File[] directories = readDirectory.getSelectedFile().listFiles();	
					
					for (int i=0; i<directories.length; i++)
					{
						String datasetName = directories[i].getName();
						if (datasetName.contains("_"))
							datasetName = datasetName.substring(0, datasetName.indexOf("_"));
						datasetName = datasetName.substring(0, datasetName.length()-1); // remove g,m,r
						controller.convertFilesAndMerge(directories[i], Dataset.getValue(datasetName), readReverse);
					}						
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(copyCleanInputFilesName))
			{				
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{
					File[] directories = readDirectory.getSelectedFile().listFiles();	
					File out = new File(readDirectory.getSelectedFile().getParentFile()+"\\clean");
					out.mkdir();
					for (int i=0; i<directories.length; i++)
					{
						File[] files = directories[i].listFiles(new Filter("clean.voter"));
						
						if (files.length>0)
						{
							FileConverterTools.copy(files[0], out);							
						}
					}						
				}
			}	
			else if (((JMenuItem)e.getSource()).getText().equals(loadMenuName))
			{
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{		            
		            loadAnalysis(readDirectory.getSelectedFile());		            
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(saveMenuName))
			{
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{		            
		            saveAnalysis(readDirectory.getSelectedFile());		            
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(gibbsMenuName))
			{						
				readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame);
				FileConverterTools.createLengthFile(readDirectory.getSelectedFile(), new File(readDirectory.getSelectedFile()+"//len.txt")); 
				/*
				if (readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{		            
		            Thread t = new Thread(){ public void run() 
		            {
		            	SimpleGibbsAlignerDemo demo = new SimpleGibbsAlignerDemo();
		            	int i = Integer.parseInt(numOfTrials.getText());
		            	int w = Integer.parseInt(motifWidth.getText());
		            	demo.runDemo(readFASTA.getSelectedFile(), false, w, i);}
		            };
		            t.start();		            
		           
				}
				*/
			}
			else if (((JMenuItem)e.getSource()).getText().equals(compareSetsName))
			{	
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{		            
		            File dir = readDirectory.getSelectedFile();
		            if (multipleDataSets)
		            	controller.compareSetsAllDir(dir);
		            else
		            	controller.compareSets(dir);											            
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(mergeCompareSetsName))
			{	
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{		            
		            File dir = readDirectory.getSelectedFile();
		            controller.mergeComputeSets(dir);											            
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(infoContentMenuName))
			{	
				if (reader.showOpenDialog(MotifAnalyzer.this.mainFrame) == JFileChooser.APPROVE_OPTION) 
				{		            
					List<String> sequences = null;
					if (inputTextArea.getText()!=null)
					{
						String temp = inputTextArea.getText();
						sequences = new ArrayList<String>();
						while (temp.contains("\n"))
						{
							String seq = temp.substring(0, temp.indexOf('\n'));
							temp = temp.substring(temp.indexOf('\n')+1, temp.length()-1);
							if (seq.contains(">")) continue;
							sequences.add(seq);							
						}
						sequences.add(temp);
		           	}
					
					Double[] score = MotifTools.getsequenceScore(reader.getSelectedFile(), sequences);
					StringBuffer buff = new StringBuffer();
					for (int i=0; i<score.length; i++)
						buff.append(score[i].doubleValue()+"\n");
					getOutputPane().setText(buff.toString());
					
					LinesSeriesChart linesChart = new LinesSeriesChart();
						
					String[] labels = new String[1];
					labels[0] = "label";
					String x[] = new String[score.length];
					Double[][] values = new Double[1][score.length];
					for (int i=0; i<score.length; i++)
					{
						values[0][i] = score[i];
						x[i] = i+1+"";
					}
					
					//linesChart.createChartImage("e://images//chart1.jpg", labels, x, values, false);
					
					tabs.add("Info Content", new ChartPanel(linesChart.createCategoryChart("","","",labels, x, values)));
														
					
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(removeRedundantSitesName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						SequenceIterator seqItr;						
						if (returnVal==JFileChooser.APPROVE_OPTION)
							seqItr = TompaDataset.getSequences(readFASTA.getSelectedFile());
						else
							seqItr = TompaDataset.getDatasetSequences(dataset);
																				
							List<Motif> m =controller.removeRedundantSites(readVoter.getSelectedFile().getCanonicalPath(), dataset);
							if (m!=null)
							{
								showMotifs(seqItr, m);								
								for (int i=0; i<controller.numOfEvaluations; i++)
									displayCharts(i);								
							}							
					}	
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(extractSitesName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						SequenceIterator seqItr;						
						if (returnVal==JFileChooser.APPROVE_OPTION)
							seqItr = TompaDataset.getSequences(readFASTA.getSelectedFile());
						else
							seqItr = TompaDataset.getDatasetSequences(dataset);
																				
							List<Motif> m =controller.extractSites(readVoter.getSelectedFile().getCanonicalPath(), dataset);
							if (m!=null)
							{
								showMotifs(seqItr, m);								
								for (int i=0; i<controller.numOfEvaluations; i++)
									displayCharts(i);								
							}							
					}	
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(drawKDistanceGraphName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						String fileName = controller.drawKDistanceGraph(readVoter.getSelectedFile(), dataset, minPoints);
						displayCharts(fileName);												
					}	
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(computeWeightName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						List<Double> list = controller.computeSimilarity(readVoter.getSelectedFile().getCanonicalPath(), dataset);
						double sim = list.get(0);
						double w = list.get(1);
						String out = String.format("Compactness=%.5f, Weight=%.5f", sim, w);
						System.out.println(out);
						JOptionPane.showMessageDialog(mainFrame, out, "Info", JOptionPane.INFORMATION_MESSAGE);						
					}	
				}
			}			
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateTompaMotifEachName))
			{				
				if (readTompa.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
				{
					controller.evaluateTompaMotifEach(readTompa.getSelectedFile().getCanonicalPath());
					displayCharts();						
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateTompaMotifAccName))
			{
				if (readTompa.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
				{
					controller.evaluateTompaMotifAcc(readTompa.getSelectedFile().getCanonicalPath());
					displayCharts();						
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateTompaMotifTotalName))
			{
												
				if (readTompa.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
				{					
					controller.evaluateTompaMotifTotal(readTompa.getSelectedFile().getCanonicalPath());
					displayCharts();						
				}						
				
			}
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateMotifVoterName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						SequenceIterator seqItr;						
						if (returnVal==JFileChooser.APPROVE_OPTION)
							seqItr = TompaDataset.getSequences(readFASTA.getSelectedFile());
						else
							seqItr = TompaDataset.getDatasetSequences(dataset);
						Vector<Motif> motifs = controller.evaluateMotifVoterEach(readVoter.getSelectedFile(), dataset);
						showMotifs(seqItr, motifs);	
						displayCharts();					
					}	
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateMotifVoterTotalName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						SequenceIterator seqItr;						
						if (returnVal==JFileChooser.APPROVE_OPTION)
							seqItr = TompaDataset.getSequences(readFASTA.getSelectedFile());
						else
							seqItr = TompaDataset.getDatasetSequences(dataset);
						Vector<Motif> motifs = controller.evaluateMotifVoterTotal(readVoter.getSelectedFile(), dataset);
						showMotifs(seqItr, motifs);	
						displayCharts();					
					}	
				}
			}
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateMVDirFirstEleName))
			{								
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
				{											
					File[] files = readDirectory.getSelectedFile().listFiles(new Filter("voter"));
					controller.evaluateMotifVoterDirFirstEle(files, org, datasetType, readReverse);
					displayCharts();					
				}	
				
			}
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateMVDirMergedName))
			{								
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
				{					
					File[] files = readDirectory.getSelectedFile().listFiles(new Filter("voter"));
					controller.evaluateMotifVoterDirMerged(files, org, datasetType, readReverse);
					displayCharts();					
				}	
				
			}
			else if (((JMenuItem)e.getSource()).getText().equals(compareMotifVoterFirstName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						File file1 = readVoter.getSelectedFile();
						if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)== JFileChooser.APPROVE_OPTION)
						{
							SequenceIterator seqItr;							
							if (returnVal==JFileChooser.APPROVE_OPTION)
								seqItr = TompaDataset.getSequences(readFASTA.getSelectedFile());
							else
								seqItr = TompaDataset.getDatasetSequences(dataset);
							Vector<Motif> motifs = new Vector<Motif>(2);
							motifs.add(FileConverterTools.readMotifVoter(file1, readReverse).firstElement());
							motifs.add(FileConverterTools.readMotifVoter(readVoter.getSelectedFile(), readReverse).firstElement());
							EvaluationController.evaluateMotifList(motifs, dataset);
							motifs.add(TompaDataset.getAnswer(dataset));
							showMotifs(seqItr, motifs);	
							displayCharts();			
						}
					}	
				}
			}		
			else if (((JMenuItem)e.getSource()).getText().equals(compareMotifVoterAllName))
			{
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						File file1 = readVoter.getSelectedFile();
						if (readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)== JFileChooser.APPROVE_OPTION)
						{
							SequenceIterator seqItr;							
							if (returnVal==JFileChooser.APPROVE_OPTION)
								seqItr = TompaDataset.getSequences(readFASTA.getSelectedFile());
							else
								seqItr = TompaDataset.getDatasetSequences(dataset);
							Vector<Motif> motifs = new Vector<Motif>(2);
							Vector<Motif> temp = FileConverterTools.readMotifVoter(file1, readReverse); 
							for (Iterator<Motif> i=temp.iterator(); i.hasNext();)
								temp.firstElement().merge(i.next());							
							temp.firstElement().setFinder(MotifFinder.MotifVoter);
							motifs.add(temp.firstElement());
							temp = FileConverterTools.readMotifVoter(readVoter.getSelectedFile(), readReverse);
							for (Iterator<Motif> i=temp.iterator(); i.hasNext();)
								temp.firstElement().merge(i.next());
							motifs.add(temp.firstElement());
							EvaluationController.evaluateMotifList(motifs, dataset);
							motifs.add(TompaDataset.getAnswer(dataset));
							showMotifs(seqItr, motifs);	
							displayCharts();			
						}
					}	
				}
			}	
			else if (((JMenuItem)e.getSource()).getText().equals(evaluateAllResultsName))
			{								
				if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
				{											
					File[] files = readDirectory.getSelectedFile().listFiles(new Filter("voter"));
					controller.evaluateMotifVoterDirAllDataset(files, datasetType, readReverse);
					displayCharts();					
				}	
				
			}
			else if (((JMenuItem)e.getSource()).getText().equals(runBioProspectorName))
				runMethod(MotifFinder.BP);
			/*else if (((JMenuItem)e.getSource()).getText().equals(runDBScanName)
					|| ((JMenuItem)e.getSource()).getText().equals(runKMeansName)
					|| ((JMenuItem)e.getSource()).getText().equals(runOPTICSName)
					|| ((JMenuItem)e.getSource()).getText().equals(runCliqueName)
					|| ((JMenuItem)e.getSource()).getText().equals(runMotifVoterName))
			{
				int type = -1;
				if (((JMenuItem)e.getSource()).getText().equals(runMotifVoterName))
					type = 0;
				else if (((JMenuItem)e.getSource()).getText().equals(runKMeansName))
					type = 1;
				else if (((JMenuItem)e.getSource()).getText().equals(runDBScanName))
					type = 2;
				else if (((JMenuItem)e.getSource()).getText().equals(runOPTICSName))
					type = 3;
				else if (((JMenuItem)e.getSource()).getText().equals(runCliqueName))
					type = 4;
				*/
			else if (Method.valueOf(((JMenuItem)e.getSource()).getText())!=null)
			{
				final Method methodType = Method.valueOf(((JMenuItem)e.getSource()).getText());
				int returnVal = JFileChooser.CANCEL_OPTION;
				if (!multipleDataSets && datasetComboBox.getSelectedItem().equals(Dataset.unknown))
					returnVal = readFASTA.showOpenDialog(MotifAnalyzer.this.mainFrame);
				
				if (multipleDataSets || !datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
				{											
					if (!multipleDataSets && readVoter.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						final SequenceIterator seqItr = (returnVal==JFileChooser.APPROVE_OPTION)?
							TompaDataset.getSequences(readFASTA.getSelectedFile())
							:TompaDataset.getDatasetSequences(dataset);
						Thread t = new Thread(){ public void run()
						{
							try{
								List<Motif> m =controller.runMotifVoter(readVoter.getSelectedFile().getCanonicalPath(), dataset, multipleDataSets, methodType, minPoints, epsilon, numOfClusters);
								if (m!=null)
								{
									showMotifs(seqItr, m);								
									for (int i=0; i<controller.numOfEvaluations; i++)
									displayCharts(i);								
								}
							}catch (Exception e){e.printStackTrace();}
						}};						
						t.start();
					}
					else if (multipleDataSets && readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame)==JFileChooser.APPROVE_OPTION)
					{
						Thread t2 = new Thread(){ public void run()
						{
							try{controller.runMotifVoter(readDirectory.getSelectedFile().getCanonicalPath(), dataset, multipleDataSets, methodType, minPoints, epsilon, numOfClusters);}
							catch (Exception e){e.printStackTrace();}
						}};
						t2.start();
					}
				
				}
			}
			}
			catch (Exception ex) 
			{
				ex.printStackTrace();
				JOptionPane.showMessageDialog(mainFrame, "An Error occured", "Error", JOptionPane.ERROR_MESSAGE);
			}
		}
		
		public void runMethod(MotifFinder method)
		{
			final JFileChooser fc = new JFileChooser();
			fc.setFileFilter(new Filter("fasta"));
			fc.setDialogTitle("Select Sequence File (in FASTA format)");
			int returnVal = JFileChooser.CANCEL_OPTION;
			if (datasetComboBox.getSelectedItem().equals(Dataset.unknown))
				returnVal = fc.showOpenDialog(MotifAnalyzer.this.mainFrame);
			
			if (!datasetComboBox.getSelectedItem().equals(Dataset.unknown) || (returnVal == JFileChooser.APPROVE_OPTION)) 
			{						
				try{
					final JFileChooser readDirectory = new JFileChooser();
					readDirectory.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
					if (readDirectory.showOpenDialog(MotifAnalyzer.this.mainFrame)== JFileChooser.APPROVE_OPTION)
					{						
						Dataset dataset = (Dataset) datasetComboBox.getSelectedItem();
						int type = datasetTypeComboBox.getSelectedIndex();
						String fileName = (returnVal==JFileChooser.APPROVE_OPTION)?
								fc.getSelectedFile().getCanonicalPath():
								TompaDataset.getDatasetSequenceFileName(dataset);
						int i = Integer.parseInt(numOfMotifs.getText());
						int w = Integer.parseInt(motifWidth.getText());
						int t = Integer.parseInt(numOfTrials.getText());
						OrganismCode organism = OrganismCode.getOrganism(organismType.getSelectedIndex());					
						log.info("Runninf method: "+method.name());
						log.info("FileName: "+fileName);
						log.info("Outputfile: "+readDirectory.getSelectedFile().getCanonicalPath());
						log.info("Dataset: "+dataset.name());
						log.info("DatasetType: "+type);
						controller.runAndEvaluate(method, fileName, i, w, t, organism, dataset, readDirectory.getSelectedFile().getCanonicalPath());
						displayCharts();
					}					
				}
				catch (Exception ex)
				{
					ex.printStackTrace();
					JOptionPane.showMessageDialog(mainFrame, "An Error occured", "Error", JOptionPane.ERROR_MESSAGE);
				}					
			}
		}
		
		public void displayCharts() throws IOException
		{
			displayCharts(0);
		}
		
		public void displayCharts(int chartNum) throws IOException
		{
			JFrame chartsFrame = new JFrame("Charts Analysis" + chartNum);
			chartsFrame.setSize(700,600);			
			JTabbedPane tabs = new JTabbedPane();			
			chartsFrame.add(tabs);
			String dataset = ((Dataset)datasetComboBox.getSelectedItem()).name();
			
			//tabs.addTab("Log", getOutputTextPane(log));
			tabs.addTab("Statistics", new JScrollPane(getStatisticsTab(chartNum)));
			String[] chartsNames = FileNames.getCharts(chartNum);
			for (int i=0; i<chartsNames.length; i++)
			{
				JPanel image1 = new JPanel(new BorderLayout());			
				image1.add(new JLabel(new ImageIcon(ImageIO.read(new File(chartsNames[i])))));	
				image1.setBackground(Color.white);
				tabs.addTab("Chart "+(i+1)+ "("+dataset+")", image1);
			}
			chartsFrame.setVisible(true);			
			
		}
		
		public void displayCharts(String fileName) throws IOException
		{
			JFrame chartsFrame = new JFrame("Charts Analysis");
			chartsFrame.setSize(700,600);			
			JTabbedPane tabs = new JTabbedPane();			
			chartsFrame.add(tabs);
				
			JPanel image1 = new JPanel(new BorderLayout());			
			image1.add(new JLabel(new ImageIcon(ImageIO.read(new File(fileName)))));	
			image1.setBackground(Color.white);
			tabs.addTab("Chart ", image1);
			
			chartsFrame.setVisible(true);			
			
		}
		
		/**
		 * Saves the statistics xls file and the 5 charts
		 * @param dir
		 */
		public void saveAnalysis(File dir)
		{
			
		}
		
		/**
		 * Loads the statistics xls file and the 5 charts
		 * @param dir
		 */
		public void loadAnalysis(File dir)
		{
			
		}
		
		public JTextArea getInputTextArea() {
			if(inputTextArea == null) {
				inputTextArea = new JTextArea();
				inputTextArea.setRows(15);
				inputTextArea.setFont(new java.awt.Font("Tahoma",1,14));
			}
			return inputTextArea;
		}
		
		private JScrollPane getJScrollPane1() {
			if(jScrollPane1 == null) {
				jScrollPane1 = new JScrollPane();
				jScrollPane1.setViewportView(getInputTextArea());
			}
			return jScrollPane1;
		}
		
		// create and return a read only scroll text
		private JScrollPane getOutputTextPane(String text)
		{
			JScrollPane jScrollPane = new JScrollPane();
			JTextPane outputPane = new JTextPane();
			outputPane.setEditable(false);
			outputPane.setText(text);
			jScrollPane.setViewportView(outputPane);
			
			return jScrollPane;
		}
		
		private JPanel getMotifScoreTab() {
			if(motifScoreTab == null) {
				motifScoreTab = new JPanel();
				motifScoreTab.setLayout(new GridLayout(1,2));
				motifScoreTab.add(getJScrollPane1());
				motifScoreTab.add(getJScrollPane2());
			}
			return motifScoreTab;
		}
		
		private JTextPane getOutputPane() {
			if(outputPane == null) {
				outputPane = new JTextPane();
			outputPane.setEditable(false);	
			outputPane.setFont(new java.awt.Font("Tahoma",1,14));
			}
			return outputPane;
		}
	
		private JScrollPane getJScrollPane2() {
			if(jScrollPane2 == null) {
				jScrollPane2 = new JScrollPane();
				jScrollPane2.setViewportView(getOutputPane());
			}
			return jScrollPane2;
		}		
		
		private Box getControlBox(SequencePanel seqPanel, String seqName)
		{
			JSlider     scale;
		    Box 		controlBox = null; 
			if (controlBox == null) 
			{
				controlBox = Box.createHorizontalBox();
		        scale      = new JSlider(SwingConstants.HORIZONTAL, 1, 100, INITIAL_SCALE);
		        
		        controlBox.add(new JLabel(seqName));
				controlBox.add(Box.createHorizontalGlue());
				controlBox.add(Box.createHorizontalStrut(10));
				controlBox.add(Box.createHorizontalGlue());
				controlBox.add(new JLabel("Scale"));
				controlBox.add(Box.createHorizontalStrut(5));
				controlBox.add(scale);
				controlBox.add(Box.createHorizontalGlue());			
			
				scale.addChangeListener(new SliderListener(seqPanel));
			}
			return controlBox;
		}

		private JLabel getJLabel1() {
			if(jLabel1 == null) {
				jLabel1 = new JLabel();
				AnchorLayout jLabel1Layout = new AnchorLayout();
				jLabel1.setLayout(jLabel1Layout);
				jLabel1.setText("    Number of Top Motifs:");
				jLabel1.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel1;
		}
		
		private JCheckBox getMultipleDataSets() {
			if(multipleDataSets == null) {
				multipleDataSets = new JCheckBox();
				AnchorLayout singleDataSetLayout = new AnchorLayout();
				multipleDataSets.setLayout(singleDataSetLayout);
				multipleDataSets.setText("Multiple Datasets");
				multipleDataSets.setSelected(false);
				multipleDataSets.setFont(new java.awt.Font("Arial",1,16));
			}
			return multipleDataSets;
		}
		
		private JTextField getNumOfMotifs() {
			if(numOfMotifs == null) {
				try
				{
					numOfMotifs = new JFormattedTextField(new MaskFormatter("##"));
				} 
				catch (ParseException e){e.printStackTrace();}
				numOfMotifs.setText("01");
			}
			return numOfMotifs;
		}
		
		private JCheckBox getCheckReverse() {
			if(checkReverse == null) {
				checkReverse = new JCheckBox();
				checkReverse.setText("CheckReverse Strand");
				checkReverse.setFont(new java.awt.Font("Arial",1,16));
			}
			return checkReverse;
		}
		
		private JLabel getJLabel2() {
			if(jLabel2 == null) {
				jLabel2 = new JLabel();
				jLabel2.setText("    Species Type:");
				jLabel2.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel2;
		}
		
		private JComboBox getSpeciesType() {
			if(organismType == null) {
				ComboBoxModel speciesTypeModel = 
					new DefaultComboBoxModel(OrganismCode.getNames());
				organismType = new JComboBox();
				organismType.setModel(speciesTypeModel);				
			}
			return organismType;
		}
		
		private JLabel getJLabel3() {
			if(jLabel3 == null) {
				jLabel3 = new JLabel();
				jLabel3.setText("    Number of trials:");
				jLabel3.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel3;
		}
		
		private JTextField getNumOfTrials() {
			if(numOfTrials == null) {
				try
				{
					numOfTrials = new JFormattedTextField(new MaskFormatter("##"));
				} catch (ParseException e){e.printStackTrace();}
				numOfTrials.setText("01");				
			}
			return numOfTrials;
		}
		
		private JLabel getJLabel4() {
			if(jLabel4 == null) {
				jLabel4 = new JLabel();
				jLabel4.setText("    Motif Width:");
				jLabel4.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel4;
		}
		
		private JTextField getMotifWidth() {
			if(motifWidth == null) {
				try
				{
					motifWidth = new JFormattedTextField(new MaskFormatter("##"));
				} catch (ParseException e){e.printStackTrace();}
				motifWidth.setText("12");
			}
			return motifWidth;
		}
		
		private JPanel getStatisticsTab() throws IOException 
		{
			return getStatisticsTab(0);
		}
		
		private JPanel getStatisticsTab(int excelNum) throws IOException
		{
			JPanel panel = new JPanel();
			GroupLayout panelLayout = new GroupLayout((JComponent)panel);
			panel.setLayout(panelLayout);
			panel.setPreferredSize(new java.awt.Dimension(676, 557));

			Object[][] tableData = FileConverterTools.readExcel(FileNames.getStatisticsExcelFileName(excelNum));
			Object columnName[] = new Object[tableData[0].length];
			for (int i=0; i< columnName.length; i++)
				columnName[i] = i;
			tableData = null; //free memory
			
			final JTable table = new JTable();
			JScrollPane jScrollPane = new JScrollPane(table);
			table.setFont(new Font(Font.DIALOG, Font.BOLD, 15));
			table.setEnabled(false);
			File file = new File(FileNames.getStatisticsExcelFileName(excelNum));			
			JTableReadTableModelTask task = new JTableReadTableModelTask(file, null, null, table);
			task.execute();
			
			JButton button = new JButton("Save to Excel");
			panelLayout.setVerticalGroup(panelLayout.createSequentialGroup()
				.addComponent(jScrollPane, GroupLayout.PREFERRED_SIZE, 485, GroupLayout.PREFERRED_SIZE)
				.addPreferredGap(LayoutStyle.ComponentPlacement.RELATED)
				.addComponent(button, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE, GroupLayout.PREFERRED_SIZE)
				.addContainerGap(40, 40));
			panelLayout.setHorizontalGroup(panelLayout.createParallelGroup()
				.addComponent(jScrollPane, GroupLayout.Alignment.LEADING, 0, 676, Short.MAX_VALUE)
				.addGroup(GroupLayout.Alignment.LEADING, panelLayout.createSequentialGroup()
				    .addGap(0, 493, Short.MAX_VALUE)
				    .addComponent(button, GroupLayout.PREFERRED_SIZE, 143, GroupLayout.PREFERRED_SIZE)
				    .addContainerGap(40, 40)));
			button.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent evt) {
					JFileChooser fc = new JFileChooser();
					int returnVal = fc.showOpenDialog(MotifAnalyzer.this.mainFrame);
					if (returnVal == JFileChooser.APPROVE_OPTION) 
					{											
						try{
						//File file = new File(fc.getSelectedFile().getCanonicalPath());
						//JTableWriteTableModelTask task = new JTableWriteTableModelTask(file, null, null, table);
						//task.execute();
						}
						catch (Exception e) {
							e.printStackTrace();						
							}
						}
					}
				});

			return panel;
		}		
		
		
		private JLabel getJLabel5() {
			if(jLabel5 == null) {
				jLabel5 = new JLabel();
				jLabel5.setText("    DataSet");
				jLabel5.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel5;
		}
		
		private JComboBox getDatasetComboBox() {
			if(datasetComboBox == null) {
				ComboBoxModel datasetComboBoxModel = 
					new DefaultComboBoxModel(Dataset.values());
				datasetComboBox = new JComboBox();
				datasetComboBox.setModel(datasetComboBoxModel);
			}
			return datasetComboBox;
		}
		
		private JLabel getJLabel6() {
			if(jLabel6 == null) {
				jLabel6 = new JLabel();
				jLabel6.setText("    Dataset Type");
				jLabel6.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel6;
		}
		
		private JComboBox getDatasetTypeComboBox() {
			if(datasetTypeComboBox == null) {
				ComboBoxModel datasetTypeModel = 
					new DefaultComboBoxModel(DatasetType.values());
				datasetTypeComboBox = new JComboBox();
				datasetTypeComboBox.setModel(datasetTypeModel);
			}
			return datasetTypeComboBox;
		}
		
		private JLabel getJLabel7() {
			if(jLabel7 == null) {
				jLabel7 = new JLabel();
				jLabel7.setText("    Epsilon");
				jLabel7.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel7;
		}
		
		private JComboBox getClusteringEpsilon() {
			if(clusteringEpsilon == null) {
				String[] range = new String[20];
				double c = 0;
				for (int i=0; i<range.length; i++, c+=0.05)
					range[i] = String.format("%.2f", c);
				ComboBoxModel clusteringEpsilonModel = 
					new DefaultComboBoxModel(range);
				clusteringEpsilon = new JComboBox();
				clusteringEpsilon.setModel(clusteringEpsilonModel);
				clusteringEpsilon.setSelectedIndex(16);
			}
			return clusteringEpsilon;
		}
		
		private JLabel getJLabel8() {
			if(jLabel8 == null) {
				jLabel8 = new JLabel();
				jLabel8.setText("    MinPoints");
				jLabel8.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel8;
		}
		
		private JComboBox getMinPoints() {
			if(minPoints == null) {
				ComboBoxModel minPointsModel = 
					new DefaultComboBoxModel(
							new Integer[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 });
				minPoints = new JComboBox();
				minPoints.setModel(minPointsModel);
				minPoints.setSelectedIndex(1);
			}
			return minPoints;
		}
		
		private JLabel getJLabel9() {
			if(jLabel9 == null) {
				jLabel9 = new JLabel();
				jLabel9.setText("    Num Of Clusters");
				jLabel9.setFont(new java.awt.Font("Arial",1,16));
			}
			return jLabel9;
		}
		
		private JComboBox getNumOfClusters() {
			if(numOfClusters == null) {
				ComboBoxModel numOfClustersModel = 
					new DefaultComboBoxModel(
							new Integer[] { -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 });
				numOfClusters = new JComboBox();
				numOfClusters.setModel(numOfClustersModel);
				numOfClusters.setSelectedIndex(0);
			}
			return numOfClusters;
		}

}
