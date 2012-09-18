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

import lemontree.utils.*;
import cern.colt.matrix.impl.*;

/**
 * Class to store all information to explain a split on a {@link TreeNode}. 
 * Contains several fields that are not (yet) used.
 * 
 * @author tomic
 *
 */
public class Split implements RealFunction, Comparable {

	/**
	 * Node to which split belongs.
	 */
	public TreeNode node;

    /**
 	 * Regulator to explain split.
 	 */
	public Gene regulator;
	
	/**
	 * Split value of regulator.
	 */
	public double splitValue;
	
	/**
	 * Sign of split.
	 */
	public int sign;
	
	/**
	 * Beta parameter of split.
	 */
	public double beta;
	
	/**
	 * Sum of regulator data in left child conditions
	 */
	public double sumDataLeft;
	
	/**
	 * Sum of regulator data in right child conditions
	 */
	public double sumDataRight;
	
	/**
	 * Number of non-NAN values for regulator in left child conditions
	 */
	public int numDataLeft;
	
	/**
	 * Number of non-NAN values for regulator in right child conditions
	 */
	public int numDataRight;
	
	/**
	 * Score for regulator - split value pair explaining this split. Score can be 
	 * conditional entropy or Bayesian assignment score.
	 */
	public double regulatorScore; 
	
	/**
	 * Score gain by making this split.
	 */
	public double scoreGain;

	
	/**
	 * Empty constructor.
	 */
	public Split() {
		node = new TreeNode();
		regulator = new Gene();
		splitValue = 0.0;
		regulatorScore = Double.NEGATIVE_INFINITY;
	}
	
	/**
	 * Construct empty split for given node.
	 */
	public Split(TreeNode node) {
		regulator = new Gene();
		splitValue = 0.0;
		regulatorScore = Double.NEGATIVE_INFINITY;
	}
	
	/**
	 * Construct split for a node with a regulator
	 * 
	 * @param lf TreeNode.
	 * @param reg regulator.
	 */
	public Split(TreeNode lf, Gene reg) {
		node = lf;
		regulator = reg;
	}
	
	/**
	 * Construct split for a node with a regulator and split value.
	 * 
	 * @param lf TreeNode.
	 * @param reg regulator.
	 * @param sv split value.
	 */
	public Split(TreeNode lf, Gene reg, double sv) {
		node = lf;
		regulator = reg;
		splitValue = sv;
		this.setSumData();
		this.computeScoreGain();
		this.computeSign();
	}

	/**
	 * Construct split for a node with a regulator, split value and entropy.
	 * 
	 * @param lf TreeNode.
	 * @param reg regulator.
	 * @param sv split value.
	 * @param ent entropy.
	 */
	public Split(TreeNode lf, Gene reg, double sv, double ent) {
		node = lf;
		regulator = reg;
		splitValue = sv;
		regulatorScore = ent;
		this.setSumData();
		this.computeScoreGain();
		this.computeSign();
	}
	
	/**
	 * Construct split for a node with a regulator, split value, entropy and score gain.
	 * 
	 * @param lf TreeNode.
	 * @param reg regulator.
	 * @param sv split value.
	 * @param ent entropy.
	 * @param sg score gain
	 */
	public Split(TreeNode lf, Gene reg, double sv, double ent, double sg) {
		node = lf;
		regulator = reg;
		splitValue = sv;
		regulatorScore = ent;
		scoreGain = sg;
		this.setSumData();
		this.computeSign();
	}
	
	
	/**
	 * Computes the data sums
	 */
	public void setSumData(){
		this.sumDataLeft = 0.0;
		this.sumDataRight= 0.0;
		this.numDataLeft = 0;
		this.numDataRight = 0;
		for (int m : this.node.leftChild.leafDistribution.condSet){
			//double x = this.regulator.normData[m];
			/*
			 * realData corresponds to unchanged expression values for condition m
			 */
			double x = this.regulator.realData[m];
			if (!Double.isNaN(x)){
				this.sumDataLeft += x;
				this.numDataLeft += 1;
			}
		}
		for (int m : this.node.rightChild.leafDistribution.condSet){
			//double x = this.regulator.normData[m];
			double x = this.regulator.realData[m];
			if (!Double.isNaN(x)){
				this.sumDataRight += x;
				this.numDataRight += 1;
			}
		}
	}
	
    /**
     * Compares 2 Splits by minus their regulatorScore (such that highest 
     * scores are ranked first after sorting)
     * 
     * @author tomic
     */
     public int compareTo(Object o) {
        Double p1 = -this.regulatorScore;
        Double p2 = -((Split)o).regulatorScore;
        return p1.compareTo(p2);
    }
	
	/**
	 * Computes the sign for this split.
	 * 
	 * @author tomic
	 */
	public void computeSign(){
		int signLeft=-1, signRight=-1;
		// notice in sumLeft, sumRight, NAN's are effectively set to the split value, such that
		// p(y|NAN)=0.5 which makes sense
		double sumLeft = this.sumDataLeft - this.splitValue*this.numDataLeft;
		double sumRight = this.sumDataRight - this.splitValue*this.numDataRight;
		if (sumLeft < 0)
			signLeft = 1;
		if (sumRight > 0)
			signRight = 1;
		if (signLeft == signRight)
			this.sign = signLeft;
		else
			this.sign = 0;
	}
	
	public void computeBeta(){
		if (this.sign == 0.0)
			this.beta = 0.0;
		else {
			// find optimal beta
			RootFinder rf = new RootFinder(this, 0.0, 20, 1E-5);
			// check if betaMax is large enough, then find root
			try {
				rf.expandXMax();
				try{
					rf.bisect();
					this.beta = rf.xRoot;
				} catch (Exception e){
					System.out.println(e);
					this.beta = 0.0;
				}
			} catch (Exception f){
				// these are some exceptional cases where split.evaluate()
				// is always zero, don't consider them
				this.beta = 0.0;
			}
		}
	}

	/**
     * Evaluates the function whose zero is the optimal value of beta
     * in the Bayesian score for this TreeNode, regulator and splitValue.
     * 
     * @author tomic
     */
    public double evaluate (double beta){
    	double sumLeft=0.0, sumRight=0.0, x;
    	for (int m : this.node.leftChild.leafDistribution.condSet){
    		x = this.regulator.normData[m]-this.splitValue;
    		if (!Double.isNaN(x)) // NAN's effectively have x=0, so don't sum
    			sumLeft += x/(1+Math.exp(-this.sign*beta*x)); 
    	}
    	for (int m : this.node.rightChild.leafDistribution.condSet){
    		x = this.regulator.normData[m]-this.splitValue;
    		if (!Double.isNaN(x)) // NAN's effectively have x=0, so don't sum
    			sumRight += x/(1+Math.exp(this.sign*beta*x)); 
    	}
    	return sumRight-sumLeft;
    }
    
    /**
     * Computes the regulator assignment score
     * 
     * @author tomic
     */
    public void computeRegulatorScore (){
    	double leftSum=0.0, rightSum=0.0;
    	// associate left conditions with -1, right with +1
    	for (int m: this.node.leftChild.leafDistribution.condSet){
    		double x = this.regulator.normData[m] - this.splitValue;
    		if (!Double.isNaN(x))
    			leftSum -= Math.log(1.0+Math.exp(this.sign*this.beta*x));
    		else // NAN's should contribute complete uncertainty (p=0.5)
    			leftSum -= Math.log(2.0);
    	}
    	for (int m: this.node.rightChild.leafDistribution.condSet){
    		double x = this.regulator.normData[m] - this.splitValue;
    		if (!Double.isNaN(x))
    			rightSum -= Math.log(1.0+Math.exp(-this.sign*this.beta*x));
    		else // NAN's should contribute complete uncertainty (p=0.5)
    			rightSum -= Math.log(2.0);
    	}
    	this.regulatorScore = leftSum+rightSum;
    }
	
    
    /**
     * Given an expression value, return left or right child of node
     * with probability defined by beta-function
     * @param xReg
     * @return
     */
    public TreeNode sampleDirection (double xReg){
    	if (this.node.nodeStatus.equals("leaf"))
    		return this.node;
    	else {
    		double pRight = 1.0/(1.0 + Math.exp(-this.sign*this.beta*(xReg-this.splitValue)));
    		Random rand = new Random();
    		if (rand.nextDouble() < pRight)
    			return this.node.rightChild;
    		else
    			return this.node.leftChild;
    	}
    }
    
    
    /**
     * Computes the gain in Bayesian score of this split.
     * 
     * @author tomic
     *
     */
    public void computeScoreGain(){
    	this.scoreGain = -this.node.leafDistribution.score + this.node.leftChild.leafDistribution.score
    		+ this.node.rightChild.leafDistribution.score;
    }
    
}
