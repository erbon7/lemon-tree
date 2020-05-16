/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/



package lemontree.modulenetwork;

import lemontree.modulenetwork.Globals; //revamp
import java.util.*;
import cern.jet.stat.*;

/**
 * Class to store experiments below a {@link TreeNode} and the their statistics.
 * 
 * @author tomic
 *
 */
public class LeafDistribution {

	/**
	 * Statistics of a {@link TreeNode}. The array order is 
	 * {number of experiments, sum of expression values, 
	 * sum of the square of expression values}.
	 */
	public double[] statistics = new double[4]; //revamp

	/**
	 * Place to store temporary statistics. For instance when testing if a 
	 * gene should be moved to another module (see also 
	 * {@link ModuleNetwork#moduleReassign(double)}).
	 */
	public double[] statisticsTmp = new double[4]; //revamp

	/**
	 * Bayesian score corresponding to the statistics.
	 */
	public double score;

	/**
	 * Place to store score corresponding to temporary statistics. 
	 */
	public double scoreTmp;

	/**
	 * List of experiments.
	 * For a leaf node: experiments grouped at this node.
	 * For an internal node: union of the experiment set of all leaves below this node.
	 */
	public ArrayList<Integer> condSet = new ArrayList<Integer>(); 

	/**
	 * Normal-gamma prior parameters. See also {@link ModuleNetwork#normalGammaPrior}.
	 */
	public double[] normalGammaPrior = new double[4]; // {lambda0, mu0, alpha0, beta0}
	
	/**
	 * Posterior parameters for the normal distribution at this node. 
	 * Array order is {mean, standard deviation}.
	 */
	public double[] gaussianParam = new double[2];


	/**
	 * Empty constructor
	 */
	public LeafDistribution(){
		
	}
	
	/**
	 * Construct empty leaf distribution but set normal-gamma prior parameters.
	 * 
	 * @param normalGammaPrior
	 */
	public LeafDistribution(double[] normalGammaPrior) {
		score = 0.0;
		condSet = new ArrayList<Integer>();
		this.normalGammaPrior = normalGammaPrior;
	}

	/**
	 * Construct leaf distribution with given statistics and normal-gamma prior 
	 * parameters.
	 * 
	 * @param statistics
	 * @param normalGammaPrior
	 */
	public LeafDistribution(double[] statistics, double[] normalGammaPrior) {
		this.statistics = statistics;
		score = 0.0;
		condSet = new ArrayList<Integer>();
		this.normalGammaPrior = normalGammaPrior;
	}


	/**
	 * Compute the Bayesian score and posterior normal distribution parameters 
	 * for this leaf distribution.
	 * 
	 * @author tomic
	 * 
	 */
	public double bayesianScore() {
		double lambda1 = normalGammaPrior[0] + statistics[0];
		double alpha1 = normalGammaPrior[2] + 0.5 * statistics[0];
		double beta1 = normalGammaPrior[3] + 0.5 * (statistics[2] - Math.pow(statistics[1], 2) / statistics[0])
				+ normalGammaPrior[0] * Math.pow(statistics[1] - normalGammaPrior[1] * statistics[0], 2)
				/ (2 * lambda1 * statistics[0]);

		score = -0.5 * statistics[0] * Math.log(2 * Math.PI) + 0.5
				* Math.log(normalGammaPrior[0]) + normalGammaPrior[2] * Math.log(normalGammaPrior[3])
				- Gamma.logGamma(normalGammaPrior[2]) + Gamma.logGamma(alpha1) - alpha1
				* Math.log(beta1) - 0.5 * Math.log(lambda1);

		gaussianParam[0] = statistics[1]/statistics[0]; // mu1; // mean
		gaussianParam[1] = Math.sqrt(statistics[2] - statistics[0]*Math.pow(gaussianParam[0],2))/Math.sqrt(statistics[0]);  // Math.sqrt(beta1 / (alpha1 - 0.5)); // std.dev.

		return score;
	}

	/**
	 * Compute the Bayesian score and posterior normal distribution parameters using
	 * {@link LeafDistribution#statisticsTmp}.
	 * 
	 * @author tomic
	 *
	 */
	public double bayesianScoreTmp() {
		double lambda1 = normalGammaPrior[0] + statisticsTmp[0];
		double alpha1 = normalGammaPrior[2] + 0.5 * statisticsTmp[0];
		double beta1 = normalGammaPrior[3] + 0.5
				* (statisticsTmp[2] - Math.pow(statisticsTmp[1], 2) / statisticsTmp[0])
				+ normalGammaPrior[0] 
				         * Math.pow(statisticsTmp[1] - normalGammaPrior[1] * statisticsTmp[0], 2)
				/ (2 * lambda1 * statisticsTmp[0]);

		scoreTmp = -0.5 * statisticsTmp[0] * Math.log(2 * Math.PI) + 0.5
				* Math.log(normalGammaPrior[0]) + normalGammaPrior[2] * Math.log(normalGammaPrior[3])
				- Gamma.logGamma(normalGammaPrior[2]) + Gamma.logGamma(alpha1) - alpha1
				* Math.log(beta1) - 0.5 * Math.log(lambda1);
		
		// debug System.err.println("SCOREtmp " + score  + " * " + statistics[3]);

		return score;
	}
	
	/**
	 * Evaluates the normal distribution
	 * @param x
	 * @return
	 * 
	 * @author tomic
	 */
	public double evaluate (double x){
		return Math.pow(Math.sqrt(2*Math.PI)*this.gaussianParam[1], -1) * 
			Math.exp(-0.5*Math.pow((x-this.gaussianParam[0])/this.gaussianParam[1], 2));
	}
}
