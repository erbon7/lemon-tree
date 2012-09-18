/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.modulenetwork;

import java.util.*;
import cern.colt.matrix.impl.*;
import cytoscape.data.annotation.OntologyTerm;

/**
 * Class to store name, description and ID number for a gene. The ID number
 * refers to the row index in the data matrix.
 * 
 * @author tomic
 *
 */
public class Gene {

	/**
	 * Name of the gene.
	 */
	public String name;

	/**
	 * Description of the gene.
	 */
	public String description;

	/**
	 * ID number of the gene, refers to the row index in the data matrix.
	 */
	public int number;

	/**
	 * ORFname
	 */
	public String ORFname;

	/**
	 * Data values, to allow different normalization for regulators.
	 */
	public double[] normData;
	
	/**
	 * Expression data values
	 */
	public double[] realData;
	
	/**
	 *  HashSet to store GO codes (BiNGO ontologyterm objects) for the gene.
	 */
	public HashSet<OntologyTerm> GOcat;

	/**
	 * for microRNA genes: seed target genes.
	 */
	public HashSet<Gene> seedTargets;
	
	/**
	 * for microRNA genes: targets predicted by various programs (targetscan, pictar, etc.)
	 */
	public HashSet<Gene> predTargets;

	/**
	 * Average expression data for this gene.
	 */
	public double mean;
	
	/**
	 * standard deviation expression data for this gene.
	 */
	public double sigma;
	
	/**
	 * The possible split values of a regulator. Used for assigning regulators to
	 * tree nodes, see for instance {@link TreeNode#findRegulator(HashSet)}.
	 * For moderately large data sets, this can simply be all expression values of 
	 * the regulator in the data set.
	 */
	public ArrayList<Double> possibleSplitValues;
	
	/**
	 * For "genes" whose expression profiles are constructed by averaging over a number of genes:
	 * a map of each gene in the average and its weight (sum of all weights = 1)
	 */
	public HashMap<Gene, Double> geneWeights;
	
	/**
	 * Boolean flag to indicate if values are discrete or continuous
	 */
	public boolean discrete;

	/**
	 * Empty constructor
	 */
	public Gene() {
		this.name = "DUMMY";
		this.number = 0;
		this.GOcat = new HashSet<OntologyTerm>();
		this.seedTargets = new HashSet<Gene>();
		this.predTargets = new HashSet<Gene>();
		this.discrete = false;
	}

	/**
	 * Construct a gene from its name and number.
	 * 
	 * @param name
	 * @param number
	 */
	public Gene(String name, int number) {
		this.name = name;
		this.number = number;
		this.GOcat = new HashSet<OntologyTerm>();
		this.seedTargets = new HashSet<Gene>();
		this.predTargets = new HashSet<Gene>();
		this.discrete = false;
	}

	/**
	 * Construct a gene from its name, description and number.
	 * 
	 * @param name
	 * @param description
	 * @param number
	 */
	public Gene(String name, String description, int number) {
		this.name = name;
		this.description = description;
		this.number = number;
		this.GOcat = new HashSet<OntologyTerm>();
		this.seedTargets = new HashSet<Gene>();
		this.predTargets = new HashSet<Gene>();
		this.discrete = false;
	}

	/**
	 * Construct a regulator gene from its name, description, number and
	 * possible split values.
	 * 
	 * @param name
	 * @param description
	 * @param number
	 * @param possibleSplitValues
	 */
	public Gene(String name, String description, int number, ArrayList<Double> possibleSplitValues) {
		this.name = name;
		this.description = description;
		this.number = number;
		this.possibleSplitValues = possibleSplitValues;
		this.GOcat = new HashSet<OntologyTerm>();
		this.seedTargets = new HashSet<Gene>();
		this.predTargets = new HashSet<Gene>();
		this.discrete = false;
	}
	
	/**
	 * Creates a gene object as a weighted average of a list of genes. This list can
	 * contain the same gene more than once, then the weights are summed. 
	 * The field {@link Gene#normData} contains the averaged expression profile.
	 * 
	 * @param splitList
	 * @param data
	 */
	public Gene (ArrayList<Split> splitList, double[][] data){
		// find total weight sum for normalization
		double norm = 0.0;
		for (int k=0; k<splitList.size(); k++)
			norm += Math.exp(splitList.get(k).regulatorScore);
		this.geneWeights = new HashMap<Gene, Double>();
		for (int k=0; k<splitList.size(); k++){
			Gene reg = splitList.get(k).regulator;
			double w = Math.exp(splitList.get(k).regulatorScore)/norm;
			if (!this.geneWeights.containsKey(reg))
				this.geneWeights.put(reg, w);
			else {
				double wold = this.geneWeights.get(reg);
				this.geneWeights.put(reg, wold + w);
			}
		}
		// set data
		this.normData = new double[data[0].length];
		for (int m=0; m<data[0].length; m++){
			for (Gene gene : this.geneWeights.keySet()){
				this.normData[m] += this.geneWeights.get(gene) * data[gene.number][m];
			}
		}
		// make name from gene with highest weight
		Collections.sort(splitList);
		this.name = String.format(splitList.get(0).regulator.name + "_%.3g", this.geneWeights.get(splitList.get(0).regulator));
	}

	// METHODS

	/**
	 * Set values for normData and realData properties.
	 */
	public void setNormData(double[][] data){
		//old code = this line only: this.normData = data[this.number];
		
		
		// pointer to all gene expression values 
		this.realData = data[this.number];
		
		// set normData arrays
		if (this.discrete == true) {
			// for discrete values, add 0.5 to discrete states
			this.normData = new double[data[this.number].length];
			for (int i=0; i < this.normData.length; i++) {
				this.normData[i] = this.realData[i] + 0.5;
			}
		}
		else {
			// continuous type, point to all expression values
			this.normData = data[this.number];
		}

		
		
		//		double sum0=0, sum1=0, sum2=0;
//		for (double x : this.normData)
//		if (!Double.isNaN(x)){
//		sum0 += 1;
//		sum1 += x;
//		sum2 += Math.pow(x, 2);
//		}
//		double mean = sum1/sum0;
//		double sigma = Math.sqrt(sum2/sum0 - Math.pow(mean, 2));
//		for (int m=0; m<this.normData.length; m++){
//		if (!Double.isNaN(this.normData[m])){
//		this.normData[m] = (this.normData[m] - mean)/sigma;
//		}
//		}
	}


	/**
	 * Computes the probability for a gene to belong to its module, given the assignment 
	 * of all other genes. The probability is the geometric average of the probabilities
	 * for each tree.
	 * 
	 * @param mod Module to which the gene belongs
	 * @param modSet set of modules to which mod belongs
	 * 
	 * @author tomic
	 */
	public double computeAssignProb (ArrayList<Module> modSet, int idx){
		if (!modSet.get(idx).genes.contains(this))
			System.out.println("Gene should be in module " + idx + " of module set!");

		DenseDoubleMatrix1D allProbs = new DenseDoubleMatrix1D(modSet.size());

		// score of module to which gene belongs
		double baseScore = 0.0;
		for (TreeNode root : modSet.get(idx).hierarchicalTrees){
			baseScore += root.bayesianScore();
		}
		// normalize
		baseScore = baseScore/(double)modSet.get(idx).hierarchicalTrees.size();

		// set weights for all modules
		allProbs.set(idx, 1.0);
		for (int k=0; k<modSet.size(); k++)
			if (k != idx){
				double score = 0.0;
				for (TreeNode root : modSet.get(k).hierarchicalTrees){
					// add gene and compute score using tmp statistics
					root.statisticsAddGeneTmp(this);
					score += root.bayesianScoreTmp();
				}
				score = score/(double)modSet.get(k).hierarchicalTrees.size();
				allProbs.set(k, Math.exp(score-baseScore));
			}
		// normalization constant
		double sum = allProbs.zSum();

		return allProbs.get(idx)/sum;
	}


	/**
	 * Divide a list of experiments into over or under expressed with respect 
	 * to a split value.
	 * 
	 * @param splitValue
	 * @param condSet list of experiments.
	 * @param data matrix where the gene's expression values can be found.
	 * 
	 * @author tomic
	 * 
	 */
	public ArrayList<ArrayList<Integer>> returnMinPlus (double splitValue, ArrayList<Integer> condSet,
			double[][] data){
		ArrayList<ArrayList<Integer>> listMinPlus = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> minList = new ArrayList<Integer>();
		ArrayList<Integer> plusList = new ArrayList<Integer>();

		for (int m : condSet) {
			if (Double.isNaN(data[this.number][m])){
				// add to plus or min with equal prob
				double x =  Math.random();
				if (x <= 0.5)
					minList.add(m);
				else
					plusList.add(m);
			} else {
				if (data[this.number][m] <= splitValue)
					minList.add(m);
				else
					plusList.add(m);
			}
		}

		listMinPlus.add(minList);
		listMinPlus.add(plusList);

		return listMinPlus;
	}

	/**
	 * Probability to be upregulated given an expression value. This assumes
	 * a Gaussian decay for the probability to have basal expression even if the 
	 * expression value is above the basal level.
	 * 
	 * @param x expression value.
	 * @param mean basal expression level (usually 0)
	 * @param sigma standard deviation for Gaussian decay of basal expression 
	 * level probability.
	 * @return number between 0 and 1.
	 * 
	 * @author tomic
	 * 
	 */
	public double upProb (double x, double mean, double sigma){
		double prob=0.0;
		if (x >= mean )
			prob = 1-Math.exp(-Math.pow((x-mean)/sigma, 2));
		return prob;
	}

	/**
	 * Probability to be downregulated given an expression value. This assumes
	 * a Gaussian decay for the probability to have basal expression even if the 
	 * expression value is below the basal level.
	 * 
	 * @param x expression value.
	 * @param mean basal expression level (usually 0)
	 * @param sigma standard deviation for Gaussian decay of basal expression 
	 * level probability.
	 * @return number between 0 and 1.
	 * 
	 * @author tomic
	 * 
	 */
	public double downProb (double x, double mean, double sigma){
		double prob=0.0;
		if (x <= mean )
			prob = 1-Math.exp(-Math.pow((x-mean)/sigma, 2));
		return prob;
	}
	
	/**
	 * Set discrete property to false or true
	 * 
	 *  @author erbon
	 */
	public void setDiscrete (boolean flag) {
		this.discrete = flag;
	}
	
	/**
	 * get function for boolean discrete property 
	 * 
	 * @return false or true
	 * 
	 * @author erbon
	 */
	public boolean getDiscrete () {
		return this.discrete;
	}
}




