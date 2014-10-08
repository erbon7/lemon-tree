/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.utils;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.font.FontRenderContext;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import javax.swing.JComponent;

import net.sf.epsgraphics.ColorMode;
import net.sf.epsgraphics.EpsGraphics;

import lemontree.modulenetwork.Gene;
import lemontree.modulenetwork.Module;
import lemontree.modulenetwork.ModuleNetwork;
import lemontree.modulenetwork.TreeNode;
import lemontree.modulenetwork.Regulator;

public class DrawModules extends JComponent {

	private static final long serialVersionUID = 1L;

	private boolean drawExperimentNames = false;
	
	private String prefix = "";
	public ArrayList<Double> regClassMean = new ArrayList<Double>();
	public ArrayList<Double> regClassSd = new ArrayList<Double>();
	public HashMap<String, Integer> extraInfo = new HashMap<String, Integer>();
	
	private ModuleNetwork modNet;
	
	/**
	 * Constructor
	 * @param m module network object.
	 */
	public DrawModules(ModuleNetwork m) {
		this.modNet = m;
	}
	
	/**
	 * Directory prefix for figure files
	 * 
	 * @param str
	 */
	public void setPrefix(String str) {
		this.prefix = str;
	}
	
	public void drawAllModules() {

		for (Module mod : modNet.moduleSet) {
			if (mod.genes.size() > 0) {
				/*
				 * Some of the hierarchical trees might not be defined, so loop until
				 * we find one that is not null to be able to draw it
				 */
				for (int nb=0;nb < mod.hierarchicalTrees.size(); nb++) {
					TreeNode t = mod.hierarchicalTrees.get(nb);
					mod.hierarchicalTree = t;
					if (mod.hierarchicalTree.nodeStatus.equals("internal")) {
						if (mod.topRegClasses.size()==0)
							System.out.println("!!! no regulators !!!");
						//String file = this.prefix + "module_" + mod.number + "_tree_" + nb;
						String file = this.prefix + "module_" + mod.number;
						System.out.println("Drawing " + file + " ...");
						drawEps(mod, file);
						break;
					}
				}
			}
		}
		System.out.println("[ok]");
	}
	
	public void drawOneModule(int modNumber) {
		Module mod = modNet.moduleSet.get(modNumber);
		if (mod.genes.size() > 0) {
			/*
			 * Some of the hierarchical trees might not be defined, so loop until
			 * we find one that is not null to be able to draw it
			 */
			for (int nb=0;nb < mod.hierarchicalTrees.size(); nb++) {
				TreeNode t = mod.hierarchicalTrees.get(nb);
				mod.hierarchicalTree = t;
				if (mod.hierarchicalTree.nodeStatus.equals("internal")) {
					if (mod.topRegClasses.size()==0)
						System.out.println("!!! no regulators !!!");
					//String file = this.prefix + "module_" + mod.number + "_tree_" + nb;
					String file = this.prefix + "module_" + mod.number;
					System.out.println("Drawing " + file + " ...");
					drawEps(mod, file);
					break;
				}
			}
		}
		System.out.println("[ok]");
	}
	
	
	public void customDraw() {

		for (Module mod : modNet.moduleSet) {
			if (mod.genes.size() > 0) {
				/*
				 * Some of the hierarchical trees might not be defined, so loop until
				 * we find one that is not null to be able to draw it
				 */
				for (int nb=0;nb < mod.hierarchicalTrees.size(); nb++) {
					TreeNode t = mod.hierarchicalTrees.get(nb);
					mod.hierarchicalTree = t;
					if (mod.hierarchicalTree.nodeStatus.equals("internal")) {
						if (mod.topRegClasses.size()==0)
							System.out.println("!!! no regulators !!!");
						//String file = this.prefix + "module_" + mod.number + "_tree_" + nb;
						String file = this.prefix + "module_" + mod.number;
						System.out.println("Drawing " + file + " ...");
						customDrawEps(mod, file);
						break;
					}
				}
			}
		}
		System.out.println("[ok]");
	}

	/**
	 * Draw a list of modules
	 * 
	 * @param moduleList
	 */
	public void customDraw(ArrayList<Integer> moduleList) {

		for (int mod_number : moduleList) {
			
			Module mod = modNet.moduleSet.get(mod_number);
			if (mod.genes.size() > 0) {
				/*
				 * Some of the hierarchical trees might not be defined, so loop until
				 * we find one that is not null to be able to draw it
				 */
				for (int nb=0;nb < mod.hierarchicalTrees.size(); nb++) {
					TreeNode t = mod.hierarchicalTrees.get(nb);
					mod.hierarchicalTree = t;
					if (mod.hierarchicalTree.nodeStatus.equals("internal")) {
						if (mod.topRegClasses.size()==0)
							System.out.println("!!! no regulators !!!");
						//String file = this.prefix + "module_" + mod.number + "_tree_" + nb;
						String file = this.prefix + "module_" + mod.number;
						System.out.println("Drawing " + file + " ...");
						customDrawEps(mod, file);
						break;
					}
				}
			}
		}
		System.out.println("[ok]");
	}

	
	/**
	 * Converts an expression value to a color. Intermediate values are a gradient between blue and yellow.
	 * Mean and sigma are passed are parameters
	 * 
	 * @param value
	 * @param mean
	 * @param sigma
	 * @return Color object
	 * 
	 * @author erbon
	 */
	private Color data2Color(double value, double mean, double sigma) {
		float yellowPerc=0.5F, bluePerc=0.5F;
		double x;

		if (value <= mean)
			x = Math.max((value - mean)/(1.5 * sigma), -0.5);
		else
			x = Math.min((value - mean)/(1.5 * sigma), 0.5);
		
		yellowPerc = new Float(0.5 + x);
		bluePerc = 1.0F - yellowPerc;	
	
		return new Color(yellowPerc, yellowPerc, bluePerc);
	}

	private Color data2ColorMir(double value, double mean, double sigma) {
		
		double x;

		if (value <= mean)
			x = Math.max((value - mean)/(1.5 * sigma), -0.5);
		else
			x = Math.min((value - mean)/(1.5 * sigma), 0.5);
		
		Color c1 = Color.blue;
        Color c2 = Color.orange;

        float ratio = (float)x + 0.5f;
        
        int red = (int)(c2.getRed() * ratio + c1.getRed() * (1 - ratio));
        int green = (int)(c2.getGreen() * ratio + c1.getGreen() * (1 - ratio));
        int blue = (int)(c2.getBlue() * ratio + c1.getBlue() * (1 - ratio));
        
        //System.out.println("col: "+value+" "+ratio+" "+red+" "+green+" "+blue);
        
        Color c = new Color(red, green, blue);
        
        return(c);
	}
	
	private Color data2ColorCNV(double value) {
		Color col = null;
		
		double lower_bound = -1.3f;
		double upper_bound = 0.7f;
		
		if (value < lower_bound)
			col = new Color(0,255,127); // green
		else if (value > upper_bound)
			col = new Color(255,0,255);  // magenta
		else {
			col = new Color(155,155,155); // grey
			
			// gradient of grey
//			double x;
//
//			if (value <= mean)
//				x = Math.max((value - mean)/(1.5 * sigma), -0.5);
//			else
//				x = Math.min((value - mean)/(1.5 * sigma), 0.5);
//			
//			Color c1 = new Color(56,56,56); // dark grey
//	        Color c2 = new Color(189,189,189); // light grey
//
//	        float ratio = (float)x + 0.5f;
//	        
//	        int red = (int)(c2.getRed() * ratio + c1.getRed() * (1 - ratio));
//	        int green = (int)(c2.getGreen() * ratio + c1.getGreen() * (1 - ratio));
//	        int blue = (int)(c2.getBlue() * ratio + c1.getBlue() * (1 - ratio));
//	        
//	        col = new Color(red, green, blue);
		}
		
		return(col);
	}
	
	/**
	 * Draw module figure.
	 * 
	 * @param filename Output file name.
	 * 
	 * @author eric
	 */
	private void drawEps(Module mod, String filename) {

		int sc = 8;         //scale factor
		int x0 = 10;		// coordinates of origin
		int y0 = 10;
		int h = 1;			// tree step height
		int width;			// figure width
		int height;			// figure height
		int x = x0;
		int ybase = y0;

		/*
		 *  get max text size for condition names
		 *  create dummy image to estimate text size
		 */
		BufferedImage tmp_img = new BufferedImage(200,200, BufferedImage.TYPE_INT_BGR);
		Graphics2D tmp_g;
		tmp_g = tmp_img.createGraphics();
		tmp_g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		tmp_g.scale(sc, sc);
		FontRenderContext frc = tmp_g.getFontRenderContext();
		float MaxTextLen = 15;
		float CondTextLen = 0;
		
		// get max condition name length
		for (int i=0;i<mod.hierarchicalTree.leafDistribution.condSet.size();i++) {
			int cond_num = mod.hierarchicalTree.leafDistribution.condSet.get(i);
			String cond_name = mod.moduleNetwork.conditionSet.get(cond_num).name;
			CondTextLen = (float) (new Font("SansSerif",Font.PLAIN,1)).getStringBounds(cond_name, frc).getWidth();
			if (CondTextLen > MaxTextLen)
				MaxTextLen = CondTextLen;
		}
		
		// set width according to the number of conditions
		width = 3 * x0 + this.modNet.data[0].length + 10;

		x = x0 + mod.hierarchicalTree.leftChild.leafDistribution.condSet.size();
		ybase = y0 + h * (mod.hierarchicalTree.treeDepth() + 1) + 2;
		
		int deltaReg = 0;
		for (int i=0;i<mod.topRegClasses.size();i++) {
			for (int j=0;j<mod.topRegClasses.get(i).size();j++)
				deltaReg += 2;
			deltaReg += 2;
		}
		
		height = ybase + deltaReg + 4 + mod.genes.size() * 2 + (int) MaxTextLen;
		
		// create graphics2D object
		try{
			EpsGraphics g2 = new EpsGraphics(mod.name, 
					new FileOutputStream(new File(filename + ".eps")),
					0,
					0,
					sc*width,
					sc*height,
					ColorMode.COLOR_CMYK);
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.scale(sc, sc);

			if (mod.hierarchicalTree.nodeStatus.equals("leaf")) {
				// do something
			} else {
				g2.setPaint(Color.white);
				g2.fillRect(0, 0, width, height);
				drawNode(mod, mod.hierarchicalTree, x, y0, h, ybase, g2);
			}
			
			g2.flush();
			g2.close();
		}
		catch (Exception e){
			e.printStackTrace();
			System.out.println(e);
			System.exit(1);
		}
	}
	
	/**
	 * Draws regulatory tree leaves with top regulators classes for NCI data
	 * 
	 * @param node TreeNode to be drawn.
	 * @param x x-coordinate of start position for vertical node line.
	 * @param y y-coordinate of start position for vertical node line.
	 * @param h height of vertical node line
	 * @param ybase y-coordinate of baseline where leaves begin.
	 * @param g graphics in which node is drawn.
	 * 
	 * @author eric
	 */
	private void drawNode(Module mod, TreeNode node, int x, int y, int h, int ybase, Graphics2D g) {
		double[][] data = mod.moduleNetwork.data; 
		if (node.nodeStatus.equals("internal")) {

			// draw vertical split line
			g.setPaint(Color.black);
			g.setStroke(new BasicStroke(1));
			g.drawLine(x, y, x, y + h);
			
			// leftChild ******************************************************
			int xl;
			if (node.leftChild.nodeStatus.equals("internal")) {
				// draw tree structure
				// if child is also internal, start should be in the middle of its own children
				xl = x - node.leftChild.rightChild.leafDistribution.condSet.size();
			} 
			else {
				// paint leaves
				int cursorY = 0;
				// if child is leaf, we draw a vertical to the middle of the leaf as well as the leaf data point (vertical line of reg tree)
				xl = x - node.leftChild.leafDistribution.condSet.size() / 2;
				g.setPaint(Color.black);
				g.setStroke(new BasicStroke(1));
				g.drawLine(xl, y + h, xl, ybase);
				int xl1 = x - node.leftChild.leafDistribution.condSet.size();

				// draw top regulators classes
				if (mod.topRegClasses.size() > 0) {
					for (int i=0;i<mod.topRegClasses.size();i++) {
						if (mod.topRegClasses.get(i).size() > 0) {
							for (Regulator reg : mod.topRegClasses.get(i)) {
								for (int m = 0; m < node.leftChild.leafDistribution.condSet.size(); m++) {
									int exp = node.leftChild.leafDistribution.condSet.get(m);
									if (!Double.isNaN(data[reg.getGene().number][exp])) {
										// convert data to color
										Color col = null;
										if (mod.moduleNetwork.getRegulatorMeanForFigures() == true)
											col = data2Color(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
										else 
											col = data2Color(data[reg.getGene().number][exp], mod.mean, mod.sigma);
										g.setPaint(col);
									} else
										// missing value
										g.setPaint(Color.white);
									plotRectangle(g, xl1+m, ybase + cursorY);
								}
								cursorY += 2;
							}
							// blank line between classes
							cursorY += 2;
						}
					}
				}
				
				//spacer 
				cursorY += 2;

				// draw the data points
				Collections.sort(node.leftChild.leafDistribution.condSet);
				for (Gene gene : mod.genes) {
					for (int m = 0; m < node.leftChild.leafDistribution.condSet.size(); m++) {
						int exp = node.leftChild.leafDistribution.condSet.get(m);
						if (!Double.isNaN(data[gene.number][exp])) {
							Color col = data2Color(data[gene.number][exp], mod.mean, mod.sigma);
							g.setPaint(col);
						} else
							// missing value
							g.setPaint(Color.white);

						plotRectangle(g, xl1+m, ybase+cursorY);
					}
					cursorY += 2;
				}

				// spacer
				cursorY += 2;

				// draw condition labels
				if (this.drawExperimentNames)
					for (int c=0;c<node.leftChild.leafDistribution.condSet.size();c++) {

						// sample color
						int num = node.leftChild.leafDistribution.condSet.get(c);
						String name = mod.moduleNetwork.conditionSet.get(num).name;
						Color col = mod.moduleNetwork.conditionSet.get(num).col;
						g.setPaint(col);
						plotRectangle(g, xl1 + c, ybase + cursorY);

						// sample name
						int deltaText = 3;
						g.setPaint(Color.black);
						g.setFont(new Font("SansSerif", Font.PLAIN, 1));
						g.rotate(Math.PI/2, xl1 + c, ybase + cursorY + deltaText);
						g.drawString(name, xl1 + c, ybase + cursorY + deltaText);
						g.rotate(-Math.PI/2, xl1 + c, ybase + cursorY + deltaText);
					}
			} // end leftchild **************************************************

			// rightChild *****************************************************
			int xr;
			float MaxGeneNameLen = 0;
			float GeneNameLen=0;
			g.setFont(new Font("SansSerif", Font.BOLD, 2));
			FontRenderContext frc = g.getFontRenderContext();

			if (node.rightChild.nodeStatus.equals("internal"))
				// draw tree structure
				// if child is also internal, start should be in the middle of its own children
				xr = x + node.rightChild.leftChild.leafDistribution.condSet.size();
			else {
				// paint leaves
				int cursorY = 0;
				
				// if child is leaf, we draw a vertical to the middle of the leaf as well as the leaf data point (vertical line of the reg tree)
				xr = x + node.rightChild.leafDistribution.condSet.size() / 2;
				g.setPaint(Color.black);
				g.setStroke(new BasicStroke(1));
				g.drawLine(xr, y + h, xr, ybase);

				// draw top regulators expression and names
				for (int i=0; i<mod.topRegClasses.size();i++) {
					if (mod.topRegClasses.get(i).size() > 0) {
						for (Regulator reg : mod.topRegClasses.get(i)) {
							GeneNameLen =  (float) g.getFont().getStringBounds(reg.getGene().name, frc).getWidth();
							if (GeneNameLen > MaxGeneNameLen)
								MaxGeneNameLen = GeneNameLen;

							for (int m = 0; m < node.rightChild.leafDistribution.condSet.size(); m++) {
								int exp = node.rightChild.leafDistribution.condSet.get(m);
								if (!Double.isNaN(data[reg.getGene().number][exp])) {
									Color col = null;
									if (mod.moduleNetwork.getRegulatorMeanForFigures() == true)
										col = data2Color(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
									else 
										col = data2Color(data[reg.getGene().number][exp], mod.mean, mod.sigma);
									g.setPaint(col);
								} else
									// missing value
									g.setPaint(Color.white);
								plotRectangle(g, x + m, ybase + cursorY);
							}
							cursorY += 2;
						}
						cursorY += 2;
					}
				}
				
				// spacer
				cursorY += 2;

				// draw the data points
				Collections.sort(node.rightChild.leafDistribution.condSet);
				for (Gene gene : mod.genes) {
					GeneNameLen =  (float) g.getFont().getStringBounds(gene.name, frc).getWidth();
					if (GeneNameLen > MaxGeneNameLen)
						MaxGeneNameLen = GeneNameLen;
					for (int m = 0; m < node.rightChild.leafDistribution.condSet.size(); m++) {
						int exp = node.rightChild.leafDistribution.condSet.get(m);
						if (!Double.isNaN(data[gene.number][exp])) {
							// convert data to color
							Color col = data2Color(data[gene.number][exp], mod.mean, mod.sigma);
							g.setPaint(col);
						} else
							// missing value
							g.setPaint(Color.white);

						plotRectangle(g, x + m, ybase + cursorY);
					}
					cursorY += 2;
				}

				cursorY += 2;

				// draw condition labels
				if (this.drawExperimentNames)
					for (int c=0;c<node.rightChild.leafDistribution.condSet.size();c++) {
						int num = node.rightChild.leafDistribution.condSet.get(c);
						String name = mod.moduleNetwork.conditionSet.get(num).name;
						Color col = mod.moduleNetwork.conditionSet.get(num).col;
						// draw rectangle for tissue type
						g.setPaint(col);
						plotRectangle(g, x+c, ybase + cursorY);

						// rotate canvas and draw string 45/90 deg.
						int deltaText = 3;
						g.setPaint(Color.black);
						g.setFont(new Font("SansSerif", Font.PLAIN, 1));
						g.rotate(Math.PI/2, x + c, ybase + cursorY + deltaText);
						g.drawString(name, x + c, ybase + cursorY + deltaText);
						g.rotate(-Math.PI/2, x + c, ybase + cursorY + deltaText);
					}

				// draw extra information, like gene names, etc.
				if (node == mod.hierarchicalTree.rightMostNode()) {

					// define extra information x coord. base
					int xtra_xbase = x + node.rightChild.leafDistribution.condSet.size() + (int) MaxGeneNameLen + 2;

					// xy is the baseline for drawing strings, while xy is the upper right corner for rectangles
					// so we have to increase y to align names with expression data
					cursorY = 2;
					
					for (int i=0;i<mod.topRegClasses.size();i++) {
						for (Regulator reg : mod.topRegClasses.get(i)) {
							
							// draw gene name
							g.setPaint(Color.black);
							g.setFont(new Font("SansSerif", Font.BOLD, 2));
							g.drawString(reg.getGene().name, new Float(x + node.rightChild.leafDistribution.condSet.size() + 1), new Float(ybase + cursorY));
							cursorY += 2;
						}
						if (mod.topRegClasses.get(i).size() > 0)
							cursorY += 2;
					}
					
					// spacer
					cursorY += 2;

					for (Gene gene : mod.genes) {
						// draw gene name
						g.setPaint(Color.black);
						g.setFont(new Font("SansSerif", Font.BOLD, 2));
						g.drawString(gene.name, new Float(x + node.rightChild.leafDistribution.condSet.size() +1), new Float(ybase + cursorY));
						
						if (this.extraInfo.get(gene.name) != null){
							Color col = null;
							if (this.extraInfo.get(gene.name)<0)
								col = new Color(0,102,255); // light blue
							else
								col = new Color(255,255,153); // light yellow
							g.setPaint(col);
							plotRectangle(g, xtra_xbase, ybase + cursorY-2);
						}
						
						cursorY += 2;
					}

				} // end extra information
			} 
			// end rightchild *************************************************

			// draw horizontal line between start coordinates for the children
			// = horizontal lines for the tree
			g.setColor(Color.black);
			g.setStroke(new BasicStroke(1));
			g.drawLine(xl, y + h, xr, y + h);

			// draw the childrens
			drawNode(mod, node.leftChild, xl, y + h, h, ybase, g);
			drawNode(mod, node.rightChild, xr, y + h, h, ybase, g);

			// draw line between leaves, first for classes of regulators, then for module genes
			int cursorY = 0;
			g.setStroke(new BasicStroke(0.1F));
			g.setPaint(Color.magenta);
			
			int startY = 0;
			for (int i=0; i<mod.topRegClasses.size();i++) {
				startY = cursorY;
				if (mod.topRegClasses.get(i).size() > 0) {
					for (Regulator reg : mod.topRegClasses.get(i)) {
						cursorY += 2;
					}
					g.drawLine(x, ybase + startY, x, ybase + cursorY);
					cursorY += 2;
				}
			}
			
			cursorY += 2;
			
			g.drawLine(x, ybase + cursorY, x, ybase + cursorY + mod.genes.size()*2);
			
		}
	}

	
	private void plotRectangle(Graphics2D g, int x, int y) {
		// fill rectangle at given coordinates on object
		g.fillRect(x, y , 1, 2);
		// draw thin black rectangle around color block
		g.setPaint(Color.black);
		g.setStroke(new BasicStroke(0.1F));
		g.drawRect(x, y , 1, 2);
	}

	
	private void plotCustomRectangle(Graphics2D g, int x, int y, int height) {
		// fill rectangle at given coordinates on object
		g.fillRect(x, y , 1, height);
		// draw thin black rectangle around color block
		//g.setPaint(Color.black);
		//g.setStroke(new BasicStroke(0.1F));
		g.drawRect(x, y , 1, height);
	}
	
	
	private void customDrawEps(Module mod, String filename) {

		int sc = 8;         //scale factor for output file
		int x0 = 10;		// coordinates of origin
		int y0 = 10;
		int h = 1;			// tree step height
		int width;			// figure width
		int height;			// figure height
		int x = x0;
		int ybase = y0;
		int deltaGeneY = 2;  // height for module genes rectangles
		

		// get max text size for condition names
		// create dummy image to estimate text size
		BufferedImage tmp_img = new BufferedImage(200,200, BufferedImage.TYPE_INT_BGR);
		Graphics2D tmp_g;
		tmp_g = tmp_img.createGraphics();
		tmp_g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		tmp_g.scale(sc, sc);
		FontRenderContext frc = tmp_g.getFontRenderContext();
		float MaxTextLen = 15;
		float CondTextLen = 0;
		
		// get max condition name length
		for (int i=0;i<mod.hierarchicalTree.leafDistribution.condSet.size();i++) {
			int cond_num = mod.hierarchicalTree.leafDistribution.condSet.get(i);
			String cond_name = mod.moduleNetwork.conditionSet.get(cond_num).name;
			CondTextLen = (float) (new Font("SansSerif",Font.PLAIN,1)).getStringBounds(cond_name, frc).getWidth();
			if (CondTextLen > MaxTextLen)
				MaxTextLen = CondTextLen;
		}
		
		// set width according to the number of conditions
		width = 3 * x0 + this.modNet.data[0].length + 10;

		// set height if there is no reg. tree / there is a reg. tree
//		if (mod.hierarchicalTree.nodeStatus.equals("leaf")){
//			height = 2 * y0 + mod.genes.size()*2;
//		}
//		else {
		
		x = x0 + mod.hierarchicalTree.leftChild.leafDistribution.condSet.size();
		ybase = y0 + h * (mod.hierarchicalTree.treeDepth() + 1) + 2;
		
		int deltaReg = 0;
		for (int i=0;i<mod.topRegClasses.size();i++) {
			for (int j=0;j<mod.topRegClasses.get(i).size();j++)
				deltaReg += 2;
			deltaReg += 2;
		}
		
		height = ybase + deltaReg + 4 + mod.genes.size() * deltaGeneY + (int) MaxTextLen;
		
//		}

		// create graphics2D object
		try{
			EpsGraphics g2 = new EpsGraphics(mod.name, 
					new FileOutputStream(new File(filename + ".eps")),
					0,
					0,
					sc*width,
					sc*height,
					ColorMode.COLOR_CMYK);
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			g2.scale(sc, sc);

			if (mod.hierarchicalTree.nodeStatus.equals("leaf")) {
				// do something
			} else {
				g2.setPaint(Color.white);
				g2.fillRect(0, 0, width, height);
				customDrawNode(mod, mod.hierarchicalTree, x, y0, h, ybase, g2, deltaGeneY);
				
				// file name
				g2.setPaint(Color.black);
				g2.setFont(new Font("SansSerif", Font.PLAIN, 1));
				//g2.drawString(filename, 10, 3);
			}
			
			g2.flush();
			g2.close();
		}
		catch (Exception e){
			e.printStackTrace();
			System.out.println(e);
			System.exit(1);
		}
	}

	private void customDrawNode(Module mod, TreeNode node, int x, int y, int h, int ybase, Graphics2D g, int deltaGeneY) {
		double[][] data = mod.moduleNetwork.data;
		
		if (node.nodeStatus.equals("internal")) {

			// draw vertical split line
			g.setPaint(Color.black);
			g.setStroke(new BasicStroke(1));
			g.drawLine(x, y, x, y + h);
			
			// leftChild ******************************************************
			int xl;
			if (node.leftChild.nodeStatus.equals("internal")) {
				// draw tree structure
				// if child is also internal, start should be in the middle of its own children
				xl = x - node.leftChild.rightChild.leafDistribution.condSet.size();
			} 
			else {
				// paint leaves
				int cursorY = 0;
				// if child is leaf, we draw a vertical to the middle of the leaf as well as the leaf data point (vertical line of reg tree)
				xl = x - node.leftChild.leafDistribution.condSet.size() / 2;
				g.setPaint(Color.black);
				g.setStroke(new BasicStroke(1));
				g.drawLine(xl, y + h, xl, ybase);
				int xl1 = x - node.leftChild.leafDistribution.condSet.size();

				// draw top regulators classes
				if (mod.topRegClasses.size() > 0) {
					for (int i=0;i<mod.topRegClasses.size();i++) {
						if (mod.topRegClasses.get(i).size() > 0) {
							for (Regulator reg : mod.topRegClasses.get(i)) {

								for (int m = 0; m < node.leftChild.leafDistribution.condSet.size(); m++) {
									int exp = node.leftChild.leafDistribution.condSet.get(m);
									if (!Double.isNaN(data[reg.getGene().number][exp])) {
										//convert data to color
										//Color col = data2Color(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
										Color col = data2ColorCNV(data[reg.getGene().number][exp]);
//										Color col = null;
//										if (i==0)
//											col = data2Color(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
//										else if (i==1)
//											col = data2ColorMir(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
//										else if (i>=2)
//											col = data2ColorCNV(data[reg.getGene().number][exp],this.regClassMean.get(i), this.regClassSd.get(i));
										//else 
										//	col = data2Color(data[reg.getGene().number][exp], this.regClassMean.get(i), this.regClassSd.get(i));
										g.setPaint(col);
									} else
										// missing value
										g.setPaint(Color.white);
									plotCustomRectangle(g, xl1+m, ybase + cursorY, 2);
								}
								cursorY += 2;
							}
							// blank line between classes
							cursorY += 2;
						}
					}
				}
				
//				if (mod.topRegulators.size() > 0) {
//					for (Gene reg : mod.topRegulators) {
//						for (int m = 0; m < node.leftChild.leafDistribution.condSet.size(); m++) {
//							int exp = node.leftChild.leafDistribution.condSet.get(m);
//							if (!Double.isNaN(data[reg.number][exp])) {
//								// convert data to color
//								Color col;
//								if (mod.moduleNetwork.getRegulatorMeanForFigures() == true)
//									col = data2Color(data[reg.number][exp], reg.mean, reg.sigma);
//								else 
//									col = data2Color(data[reg.number][exp], mod.mean, mod.sigma);
//								g.setPaint(col);
//							}
//							else
//								// missing value
//								g.setPaint(Color.white);
//
//							plotRectangle(g,xl1+m, ybase + cursorY);
//						}
//						cursorY += 2;
//					}
//				}
				
				//spacer 
				cursorY += 2;

				// draw the data points
				Collections.sort(node.leftChild.leafDistribution.condSet);
				for (Gene gene : mod.genes) {
					for (int m = 0; m < node.leftChild.leafDistribution.condSet.size(); m++) {
						int exp = node.leftChild.leafDistribution.condSet.get(m);
						if (!Double.isNaN(data[gene.number][exp])) {
							Color col = data2Color(data[gene.number][exp], mod.mean, mod.sigma);
							g.setPaint(col);
						} else
							// missing value
							g.setPaint(Color.white);

						plotCustomRectangle(g, xl1+m, ybase+cursorY, deltaGeneY);
					}
					cursorY += deltaGeneY;
				}

				// spacer
				cursorY += deltaGeneY;

				// draw condition labels
				if (this.drawExperimentNames)
					for (int c=0;c<node.leftChild.leafDistribution.condSet.size();c++) {

						// sample color
						int num = node.leftChild.leafDistribution.condSet.get(c);
						String name = mod.moduleNetwork.conditionSet.get(num).name;
						Color col = mod.moduleNetwork.conditionSet.get(num).col;
						g.setPaint(col);
						plotCustomRectangle(g, xl1 + c, ybase + cursorY, deltaGeneY);

						// sample name
						int deltaText = 3;
						g.setPaint(Color.black);
						g.setFont(new Font("SansSerif", Font.PLAIN, 1));
						g.rotate(Math.PI/2, xl1 + c, ybase + cursorY + deltaText);
						g.drawString(name, xl1 + c, ybase + cursorY + deltaText);
						g.rotate(-Math.PI/2, xl1 + c, ybase + cursorY + deltaText);
					}
			} // end leftchild **************************************************

			// rightChild *****************************************************
			int xr;
			float MaxGeneNameLen = 0;
			float GeneNameLen=0;
			g.setFont(new Font("SansSerif", Font.BOLD, 2));
			FontRenderContext frc = g.getFontRenderContext();

			if (node.rightChild.nodeStatus.equals("internal"))
				// draw tree structure
				// if child is also internal, start should be in the middle of its own children
				xr = x + node.rightChild.leftChild.leafDistribution.condSet.size();
			else {
				// paint leaves
				int cursorY = 0;
				
				// if child is leaf, we draw a vertical to the middle of the leaf as well as the leaf data point (vertical line of the reg tree)
				xr = x + node.rightChild.leafDistribution.condSet.size() / 2;
				g.setPaint(Color.black);
				g.setStroke(new BasicStroke(1));
				g.drawLine(xr, y + h, xr, ybase);

				// draw top regulators expression and names
				for (int i=0; i<mod.topRegClasses.size();i++) {
					if (mod.topRegClasses.get(i).size() > 0) {
						for (Regulator reg : mod.topRegClasses.get(i)) {
							GeneNameLen =  (float) g.getFont().getStringBounds(reg.getGene().name, frc).getWidth();
							if (GeneNameLen > MaxGeneNameLen)
								MaxGeneNameLen = GeneNameLen;

							for (int m = 0; m < node.rightChild.leafDistribution.condSet.size(); m++) {
								int exp = node.rightChild.leafDistribution.condSet.get(m);
								if (!Double.isNaN(data[reg.getGene().number][exp])) {
									//Color col = data2Color(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
									Color col = data2ColorCNV(data[reg.getGene().number][exp]);
//									Color col = null;
//									if (i==0)
//										col = data2Color(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
//									else if (i==1)
//										col = data2ColorMir(data[reg.getGene().number][exp], reg.getGene().mean, reg.getGene().sigma);
//									else if (i>=2)
//										col = data2ColorCNV(data[reg.getGene().number][exp], this.regClassMean.get(i), this.regClassSd.get(i));
									//else 
									//	col = data2Color(data[reg.getGene().number][exp], this.regClassMean.get(i), this.regClassSd.get(i));
									g.setPaint(col);
								} else
									// missing value
									g.setPaint(Color.white);
								plotCustomRectangle(g, x + m, ybase + cursorY, 2);
							}
							cursorY += 2;
						}
						cursorY += 2;
					}
				}
				
//				if (mod.topRegulators.size() > 0) {
//					for (Gene reg : mod.topRegulators) {
//						for (int m = 0; m < node.rightChild.leafDistribution.condSet.size(); m++) {
//							int exp = node.rightChild.leafDistribution.condSet.get(m);
//							if (!Double.isNaN(data[reg.number][exp])) {
//								// convert data to color
//								Color col;
//								if (mod.moduleNetwork.getRegulatorMeanForFigures() == true)
//									col = data2Color(data[reg.number][exp], reg.mean, reg.sigma);
//								else
//									col = data2Color(data[reg.number][exp], mod.mean, mod.sigma);
//								g.setPaint(col);
//							}
//							else
//								// missing value
//								g.setPaint(Color.white);
//
//							plotRectangle(g, x+m, ybase + cursorY);
//						}
//						cursorY += 2;
//					}
//				}
				
				// spacer
				cursorY += 2;

				// draw the data points
				Collections.sort(node.rightChild.leafDistribution.condSet);
				for (Gene gene : mod.genes) {
					GeneNameLen =  (float) g.getFont().getStringBounds(gene.name, frc).getWidth();
					if (GeneNameLen > MaxGeneNameLen)
						MaxGeneNameLen = GeneNameLen;
					for (int m = 0; m < node.rightChild.leafDistribution.condSet.size(); m++) {
						int exp = node.rightChild.leafDistribution.condSet.get(m);
						if (!Double.isNaN(data[gene.number][exp])) {
							// convert data to color
							Color col = data2Color(data[gene.number][exp], mod.mean, mod.sigma);
							g.setPaint(col);
						} else
							// missing value
							g.setPaint(Color.white);

						plotCustomRectangle(g, x + m, ybase + cursorY, deltaGeneY);
					}
					cursorY += deltaGeneY;
				}

				cursorY += deltaGeneY;

				// draw condition labels
				if (this.drawExperimentNames)
					for (int c=0;c<node.rightChild.leafDistribution.condSet.size();c++) {
						int num = node.rightChild.leafDistribution.condSet.get(c);
						String name = mod.moduleNetwork.conditionSet.get(num).name;
						Color col = mod.moduleNetwork.conditionSet.get(num).col;
						// draw rectangle for tissue type
						g.setPaint(col);
						//plotCustomRectangle(g, x+c, ybase + cursorY);

						// rotate canvas and draw string 45/90 deg.
						int deltaText = 3;
						g.setPaint(Color.black);
						g.setFont(new Font("SansSerif", Font.PLAIN, 1));
						g.rotate(Math.PI/2, x + c, ybase + cursorY + deltaText);
						g.drawString(name, x + c, ybase + cursorY + deltaText);
						g.rotate(-Math.PI/2, x + c, ybase + cursorY + deltaText);
					}

				// draw extra information, like gene names, etc.
				if (node == mod.hierarchicalTree.rightMostNode()) {

					// define extra information x coord. base
					int xtra_xbase = x + node.rightChild.leafDistribution.condSet.size() + (int) MaxGeneNameLen + 2;

					// xy is the baseline for drawing strings, while xy is the upper right corner for rectangles
					// so we have to increase y to align names with expression data
					cursorY = 2;
					
					for (int i=0;i<mod.topRegClasses.size();i++) {
						for (Regulator reg : mod.topRegClasses.get(i)) {
							
							// draw gene name
							g.setPaint(Color.black);
							//g.setFont(new Font("SansSerif", Font.BOLD, 2));
							//g.drawString(reg.getGene().name, new Float(x + node.rightChild.leafDistribution.condSet.size() + 1), new Float(ybase + cursorY));
							
							g.setFont(new Font("SansSerif", Font.BOLD, 2));
							g.drawString(reg.getGene().name, new Float(x + node.rightChild.leafDistribution.condSet.size() + 1), new Float(ybase + cursorY));
							
							if (this.extraInfo.get(reg.getGene().name) != null){
								Color col = null;
								if (this.extraInfo.get(reg.getGene().name)<0)
									col = new Color(0,102,255); // light blue
								else
									col = new Color(255,255,153); // light yellow
								g.setPaint(col);
								//plotCustomRectangle(g, xtra_xbase, ybase + cursorY-2, 2);
							}
							// draw score ratio indicator
//							double scoreRatio = reg.getScore() / mod.maxTopRegScore;
//							// color like magenta, very light
//							Color col = new Color(255, 179, 209);
//							g.setPaint(col);
//							g.fillRect(xtra_xbase, ybase + cursorY-2, 5, 2);
//							// compute score as percentage of the width of the colored box
//							int w = (int) Math.round(scoreRatio * 5);
//							// dark magenta like for score percentage
//							col = new Color(255, 0, 12);
//							g.setPaint(col);
//							g.fillRect(xtra_xbase, ybase+cursorY-2, w, 2);
//							g.drawRect(xtra_xbase, ybase + cursorY - 2, 5, 2);
							
							cursorY += 2;
						}
						if (mod.topRegClasses.get(i).size() > 0)
							cursorY += 2;
							//cursorY += 4;
					}
					
//					for (Gene reg : mod.topRegulators) {
//						// draw gene name
//						g.setPaint(Color.black);
//						g.setFont(new Font("SansSerif", Font.BOLD, 2));
//						g.drawString(reg.name, new Float(x + node.rightChild.leafDistribution.condSet.size() + 1), new Float(ybase + cursorY));
//						
//						cursorY += 2;
//					}
					
					// spacer
					cursorY += 2;

					for (Gene gene : mod.genes) {
						// draw gene name
						g.setPaint(Color.black);
						g.setFont(new Font("SansSerif", Font.BOLD, 2));
						g.drawString(gene.name, new Float(x + node.rightChild.leafDistribution.condSet.size() +1), new Float(ybase + cursorY));
						
//						if (this.extraInfo.get(gene.name) != null){
//							Color col = null;
//							if (this.extraInfo.get(gene.name)<0)
//								col = new Color(0,102,255); // light blue
//							else
//								col = new Color(255,255,153); // light yellow
//							g.setPaint(col);
//							plotCustomRectangle(g, xtra_xbase, ybase + cursorY-2);
//						}
						
						cursorY += deltaGeneY;
					}

				} // end extra information
			} 
			// end rightchild *************************************************

			// draw horizontal line between start coordinates for the children
			// = horizontal lines for the tree
			g.setColor(Color.black);
			g.setStroke(new BasicStroke(1));
			g.drawLine(xl, y + h, xr, y + h);

			// draw the childrens
			customDrawNode(mod, node.leftChild, xl, y + h, h, ybase, g, deltaGeneY);
			customDrawNode(mod, node.rightChild, xr, y + h, h, ybase, g, deltaGeneY);

			// draw line between leaves, first for classes of regulators, then for module genes
			int cursorY = 0;
			g.setStroke(new BasicStroke(0.1F));
			g.setPaint(Color.magenta);
			
			int startY = 0;
			for (int i=0; i<mod.topRegClasses.size();i++) {
				startY = cursorY;
				if (mod.topRegClasses.get(i).size() > 0) {
					for (Regulator reg : mod.topRegClasses.get(i)) {
						cursorY += 2;
					}
					g.drawLine(x, ybase + startY, x, ybase + cursorY);
					cursorY += 2;
				}
			}
			
//			cursorY += mod.topRegulators.size() * 2;
//			g.drawLine(x, ybase, x, ybase + cursorY);
			
			cursorY += 2;
			
			g.drawLine(x, ybase + cursorY, x, ybase + cursorY + mod.genes.size() * deltaGeneY);
			
		}
	}



}
