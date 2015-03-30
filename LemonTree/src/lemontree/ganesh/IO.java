/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.ganesh;

import java.io.*;

public class IO {

	/**
	 * List of genes
	 */
	public String[][] geneSet;
	
	/**
	 * List of conditions
	 */
	public String[] conditionSet;
	
	/**
	 * Number of genes
	 */
	public int numGenes;
	
	/**
	 * Number of conditions for expression data
	 */
	public int numCond;

	/**
	 * Expression data matrix
	 */
	public double[][] data;
	
	public GibbsSampler Cinitial;
	public double[] normalGammaPrior = new double[4];
	Boolean outputMatrix;

	public IO(int initClust, String dataDir, String datasetFile, String outputfile,
			int numRuns, int numSteps, int burnIn, int sampleStep,
			Boolean oneway, Boolean outputMatrix) {
		this.readExpressionMatrix(dataDir, datasetFile);
		this.outputMatrix = outputMatrix;
		if (oneway)
			this.cluster1way(numRuns, numSteps, burnIn, sampleStep, outputfile);
		else
			this.cluster2way(initClust, numRuns, numSteps, burnIn, sampleStep, outputfile);
	}

	/**
	 * reads data matrix
	 * @param dataDir path to the directory containing datafile
	 * @param datasetfile filename containing data
	 */
	public void readExpressionMatrix(String dataDir, String datasetFile) {
		BufferedReader bufferedReader;
		try {
			bufferedReader = new BufferedReader(new FileReader(new File(dataDir, datasetFile)));
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}

		// parse first line and get number of conditions
		String firstLine="";
		try {
			firstLine = bufferedReader.readLine();
		} catch (IOException e) {
			System.out.println("error processing first line");
			e.printStackTrace();
			System.exit(1);
		}
		String[] tokens = firstLine.trim().split("\\t");
		//Columns UID and DESCRIPTION don't count
		this.numCond = tokens.length - 2;
		// the condition set
		this.conditionSet = new String[numCond];
		for (int i = 2; i < tokens.length; i++) {
			this.conditionSet[i - 2] = tokens[i];
		}

		// count number of genes and append each line to an expandable list	
		StringBuffer strbuf = new StringBuffer();
		String newLineOfText;
		int numGenes = 0;
		try {
			while ((newLineOfText = bufferedReader.readLine()) != null) {
				strbuf.append(newLineOfText + "\n");
				numGenes++;
			}
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		this.numGenes = numGenes;

		// the data array
		this.data = new double[numGenes][numCond];
		// the gene set; stores name and description of genes 
		this.geneSet = new String[2][numGenes];

		// fill the arrays
		int numMissing = 0;
		String text = new String(strbuf);
		String[] lines = text.trim().split("\\n");
		for (int i = 0; i < numGenes; i++) {
			String[] values = new String[numCond + 2];
			String[] tmp = lines[i].split("\\t");
			for (int i2 = 0; i2 < tmp.length; i2++) {
				values[i2] = tmp[i2];
			}
			for (int i2 = tmp.length; i2 < values.length; i2++) {
				values[i2] = "";
			}
			this.geneSet[0][i] = values[0].trim();
			this.geneSet[1][i] = values[1].trim();
			for (int j = 2; j < numCond + 2; j++) {
				try {
					this.data[i][j - 2] = new Double(values[j].trim())
							.doubleValue();
				} catch (NumberFormatException e) {
					data[i][j - 2] = Double.NaN;
					numMissing++;
				}
			}
		}

		try {
			bufferedReader.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("... found " + numGenes
				+ " genes, with expression values for " + numCond
				+ " conditions (but " + numMissing + " out of " + numGenes
				* numCond + " expression values were missing.)\n");

	}

	/**
	 * clusters data in both directions (genes and conditions)
	 * 
	 * @param initial_clusters initial number of clusters
	 * @param numRuns number of run of the Gibbs sampler
	 * @param numSteps number of steps for every run
	 * @param burnIn number of burnin steps
	 * @param sampleStep number of Gibbs sampler iterations between two samples
	 * @param outputfile output file name
	 */
	public void cluster2way(int initial_clusters, int numRuns, int numSteps,
			int burnIn, int sampleStep, String outputfile) {
		int counter1 = 1;
		for (int n = 0; n < numRuns; n++) {
			System.out.println("First run");
			int counter2 = 1;
			this.Cinitial = new GibbsSampler(data);
			Cinitial.RandomAssignTwoway(initial_clusters);
			Cinitial.Cluster2way();
			for (int k = 0; k < numSteps; k++) {
				if (k < burnIn) {
					Cinitial.Cluster2way();
				} else {
					Cinitial.Cluster2way();
					if (k % sampleStep < 1) {
						String s = outputfile + "_run" + counter1 + "_sample"
								+ counter2;
						if (outputMatrix)
							this.outputMatrix(s);
						else
							this.output(s);
						counter2++;
					}
				}
			}
			counter1++;
		}
	}

	/**
	 * cluster data in one direction
	 * @param numRuns number of run of the Gibbs sampler
	 * @param numSteps number of steps for every run
	 * @param burnIn number of burnin steps
	 * @param sampleStep number of Gibbs sampler iterations between two samples
	 * @param outputfile output file name
	 */
	public void cluster1way(int numRuns, int numSteps, int burnIn,
			int sampleStep, String outputfile) {
		int counter1 = 1;
		for (int n = 0; n < numRuns; n++) {
			System.out.println("First run");
			int counter2 = 1;
			this.Cinitial = new GibbsSampler(data);
			Cinitial.AssignRandomrowEachcolumn();
			Cinitial.Cluster();
			for (int k = 0; k < numSteps; k++) {
				if (k < burnIn) {
					Cinitial.Cluster();
				} else {
					Cinitial.Cluster();
					if (k % sampleStep < 1) {
						String s = outputfile + "_run" + counter1 + "_sample"
								+ counter2;
						if (outputMatrix)
							this.outputMatrix(s);
						else
							this.output(s);
						counter2++;
					}
				}
			}
			counter1++;
		}
	}

	/**
	 * Writes output in output file
	 * @param outputfile name of the file to which output is written
	 */
	public void output(String outputFile) {
		File f = new File(outputFile);
		try {
			FileWriter fw = new FileWriter(f);
			PrintWriter pw = new PrintWriter(fw);

			// Print header
			for (int i = 0; i < Cinitial.ClusterSet.size(); i++)
				for (Integer j : Cinitial.ClusterSet.get(i).RowSet)
					pw.println(geneSet[0][j] + "\t" + i);

			fw.close();
		} catch (IOException e) {
			System.out.println("IOException: " + e);
			e.printStackTrace();
		}
	}

	public void outputMatrix(String outputFile) {
		File f = new File(outputFile);
		try {
			FileWriter fw = new FileWriter(f);
			PrintWriter pw = new PrintWriter(fw);

			// Print header
			for (int i = 0; i < Cinitial.ClusterSet.size(); i++)
				for (Integer j : Cinitial.ClusterSet.get(i).RowSet)
					for (Integer k : Cinitial.ClusterSet.get(i).RowSet)
						pw.println(j + "\t" + k + "\t" + i);
			//           			pw.println(geneSet[0][j]+"\t"+geneSet[0][k]+"\t"+i);  

			fw.close();
		} catch (IOException e) {
			System.out.println("IOException: " + e);
			e.printStackTrace();
		}
	}

}
