/*
 * NestedMICA Motif Inference Toolkit
 *
 * Copyright (c) 2004-2007: Genome Research Ltd.
 *
 * NestedMICA is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * or see the on-line version at http://www.gnu.org/copyleft/lgpl.txt
 *
 */
 
package org.biojava.bio.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.RenderingHints;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.WeightMatrix;
import org.biojava.bio.gui.DNAStyle;
import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.gui.TextLogoPainter;

/**
 * Simple weight-matrix viewer
 *
 * @author Thomas Down
 */

@SuppressWarnings("serial")
public class WMPanel extends JPanel {
    private WeightMatrix wm;
    private DistributionLogo[] logos;
    
    public WMPanel(WeightMatrix wm) {
        super();
        this.wm = wm;
        setBackground(Color.white);
        
        RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        
        try {
            setLayout(new GridLayout(1, wm.columns()));
            logos = new DistributionLogo[wm.columns()];
            for (int pos = 0; pos < wm.columns(); ++pos) {
                Distribution dist = wm.getColumn(pos);
                DistributionLogo dl = new DistributionLogo();
                dl.setScaleByInformation(true);
                dl.setRenderingHints(hints);
                dl.setBackground(Color.white);
                dl.setOpaque(true);
                dl.setDistribution(dist);
                dl.setPreferredSize(new Dimension(40, 50));
                dl.setLogoPainter(new TextLogoPainter());
                dl.setStyle(new DNAStyle());
                add(dl);
                logos[pos] = dl;
            }
        } catch (BioException ex) {
            throw new BioError(ex);
        }
    }
    
    public void setMatrix(WeightMatrix wm) {
        this.wm = wm;
        try {
            for (int pos = 0; pos < wm.columns(); ++pos) {
                logos[pos].setDistribution(wm.getColumn(pos));
            }
        } catch (BioException ex) {
            throw new BioError(ex);
        }
        repaint();
    }
    
    public WeightMatrix getMatrix() {
        return wm;
    }
    
    public static void wmViewer(WeightMatrix wm, String message) 
    {
        WMPanel wmv = new WMPanel(wm);
        JFrame frame = new JFrame("Weight matrix viewer" + ((message == null) ? "" : (" (" + message + ")")));
        frame.getContentPane().add(wmv);
        frame.pack();
        frame.setVisible(true);
    }
}
