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

import lemontree.modulenetwork.Globals; //revamp
import lemontree.utils.*;
import nu.xom.Element;

/**
 * Class for nodes in the hierarchical trees representing regulation programs
 * of modules. 
 * 
 * @author tomic
 *
 */
public class TreeNode implements Comparable {
    
    private final String INTERNAL = "internal";
    
    private final String LEAF = "leaf";
    
    /**
     * This node belongs to a tree for a certain module.
     */
    public Module module;
    
    /**
     * Pointer to the parent of this node in the tree.
     */
    public TreeNode parentNode;
    
    /**
     * INTERNAL or LEAF
     */
    public String nodeStatus;
    
     /**
     * Regulator information for this node {@link TreeNode#testGene}, 
     * {@link TreeNode#splitValue}, {@link TreeNode#entropy}).
     */
    public Split regulationSplit;
    
    /**
     * Score for merging the children of this node in the hierarchical tree.
     * See also {@link Module#hierarchicalClustering(List, boolean)}.
     */
    public double mergeScore;

    /**
     * Highest scoring partition of the experiments below this node.
     */
    public boolean partition;
    
    /**
     * Statistical properties of the data below this node.
     */
    public LeafDistribution leafDistribution;
    
    /**
     * List of splits chosen by 
     * {@link TreeNode#findRegulatorStoch(int, HashSet, double)}
     */
    public ArrayList<Split> testSplits;

    /**
     * Lists of random splits of the same size as {@link TreeNode#testSplits}
     * for statistical comparisons.
     */
    public ArrayList<Split> testSplitsRandom;
    
	/**
	 * Map of all regulators to their weight.
	 */
	public HashMap<Gene, Double> regulatorWeights;

	/**
	 * Map of all random regulators to their weight.
	 */
	public HashMap<Gene, Double> regulatorWeightsRandom;

	
    /**
     * Pointer to the left child node of this node in the tree.
     */
    public TreeNode leftChild;
    
    /**
     * Pointer to the right child node of this node in the tree.
     */
    public TreeNode rightChild;
    
    
    // CONSTRUCTOR
    
    public TreeNode() {
        this.nodeStatus = LEAF;
        this.leafDistribution = new LeafDistribution();
    }
    
    public TreeNode(String nodeStatus, double[] normalGammaPrior) {
        this.nodeStatus = nodeStatus;
        this.leafDistribution = new LeafDistribution(normalGammaPrior);
    }
    
    public TreeNode(Module mod, double[] normalGammaPrior) {
    	this.module = mod;
        this.nodeStatus = LEAF;
        this.leafDistribution = new LeafDistribution(normalGammaPrior);
    }
    
    public TreeNode(Module mod, String nodeStatus, double[] normalGammaPrior) {
    	this.module = mod;
        this.nodeStatus = nodeStatus;
        this.leafDistribution = new LeafDistribution(normalGammaPrior);
    }
    
    /**
     * Merges conditions of a list of nodes into one leaf node
     * @param mod
     * @param nodeList
     * @param normalGammaPrior
     */
    public TreeNode(Module mod, List<TreeNode> nodeList, double[] normalGammaPrior) {
    	this.module = mod;
    	this.nodeStatus = LEAF;
        this.leafDistribution = new LeafDistribution(normalGammaPrior);
        for (TreeNode node : nodeList)
        	for (int m : node.leafDistribution.condSet)
        		this.leafDistribution.condSet.add(m);
        this.statistics();
		this.leafDistribution.bayesianScore();
    }
    
    /**
     * Creates internal node by joining two children
     * @param mod
     * @param leftChild
     * @param rightChild
     * @param normalGammaPrior
     */
    public TreeNode(Module mod, TreeNode leftChild, TreeNode rightChild, TreeNode parentNode, double[] normalGammaPrior){
    	this.module = mod;
    	this.nodeStatus = LEAF;
        this.leafDistribution = new LeafDistribution(normalGammaPrior);
        for (int m : leftChild.leafDistribution.condSet)
        	this.leafDistribution.condSet.add(m);
        for (int m : rightChild.leafDistribution.condSet)
        	this.leafDistribution.condSet.add(m);
        this.statistics();
		this.leafDistribution.bayesianScore();
		this.nodeStatus = INTERNAL;
		this.parentNode = parentNode;
        this.leftChild = leftChild;
        this.leftChild.parentNode = this;
        this.rightChild = rightChild;
        this.rightChild.parentNode = this;
    }
    
    public TreeNode(TreeNode parentNode) {
    	this.module = parentNode.module;
        this.parentNode = parentNode;
        this.nodeStatus = LEAF;
        this.leafDistribution = new LeafDistribution(parentNode.leafDistribution.normalGammaPrior);
    }
    
    public TreeNode(Module mod, Gene testGene, double splitValue, String nodeStatus,
            TreeNode leftChild, TreeNode rightChild, double[] normalGammaPrior) {
    	this.module = mod;
    	this.regulationSplit = new Split(this, testGene, splitValue);
        this.nodeStatus = nodeStatus;
        this.leftChild = leftChild;
        this.rightChild = rightChild;
        if (nodeStatus.equals(LEAF))
            this.leafDistribution = new LeafDistribution(normalGammaPrior);
    }
    
    
    /**
     * Compares 2 TreeNodes by their mergeScore, used in 
     * {@link Module#hierarchicalClustering(List, boolean)}.
     * 
     * @author tomic
     */
     public int compareTo(Object o) {
        Double p1 = this.mergeScore;
        Double p2 = ((TreeNode)o).mergeScore;
        return p1.compareTo(p2);
    }
    
    
    /**
     * Depth is the largest distance to a leaf node.
     *   
     * @return depth.
     * @author tomic
     */
     public int treeDepth() {
        int depth;
        if (this.nodeStatus.equals(LEAF))
            depth = 0;
        else
            depth = 1 + Math.max(this.leftChild.treeDepth(), this.rightChild.treeDepth());
        return depth;
    }
    
    /**
     * Level of TreeNode is the distance to the root. 
     * 
     * @return level.
     * @author tomic
     */
     public int treeLevel(){
        int level;
        if (this.parentNode instanceof TreeNode)
            level = 1 + this.parentNode.treeLevel();
        else
            level = 0;
        return level;
     }	
    
      
     /**
      * Returns the weight of the path to a leaf in the mixture distribution
      * 
      * @return
      * 
      * @author tomic
      */
     public double pathWeight(){
    	 double weight;
    	 if (this.parentNode instanceof TreeNode){
    		 if (this.equals(this.parentNode.rightChild)){
     			 Split splt = this.parentNode.regulationSplit;
     			 double p = Math.exp(splt.regulatorScore/(double)splt.node.leafDistribution.condSet.size());
    			 weight = p * this.parentNode.pathWeight();
    		 } else {
     			 Split splt = this.parentNode.regulationSplit;
     			 double p = Math.exp(splt.regulatorScore/(double)splt.node.leafDistribution.condSet.size());
    			 weight = p * this.parentNode.pathWeight();
    		 }
    	 } else {
    		 weight = 1.0;
    	 }
    	 return weight;
     }
     
     /**
      * Returns the weight of the path to a leaf in the mixture distribution for
      * a specific experiment
      * 
      * @return
      * 
      * @author tomic
      */
     public double pathWeight(double[] xExp){
    	 double weight;
    	 if (this.parentNode instanceof TreeNode){
    		 if (this.equals(this.parentNode.rightChild)){
     			 Split splt = this.parentNode.regulationSplit;
     			 double p = 1.0/(1.0 + Math.exp(-splt.sign*splt.beta*(xExp[splt.regulator.number]-splt.splitValue)));
    			 weight = p * this.parentNode.pathWeight(xExp);
    		 } else {
     			 Split splt = this.parentNode.regulationSplit;
     			 double p = 1.0/(1.0 + Math.exp(splt.sign*splt.beta*(xExp[splt.regulator.number]-splt.splitValue)));
    			 weight = p * this.parentNode.pathWeight(xExp);
    		 }
    	 } else {
    		 weight = 1.0;
    	 }
    	 return weight;
     }
     
     /**
      * Compute regulator weights
      * 
      * @author tomic	
      */
     public void setRegulatorWeights (){
    	 this.regulatorWeights = new HashMap<Gene, Double>();
    	 double fact = (double)this.leafDistribution.condSet.size()/(double)this.module.moduleNetwork.data[0].length;
    	 for (Split splt : this.testSplits){
				if (!this.regulatorWeights.containsKey(splt.regulator))
					this.regulatorWeights.put(splt.regulator, 0.0);
				double w = this.regulatorWeights.get(splt.regulator)
					+ Math.exp(splt.regulatorScore/this.leafDistribution.condSet.size())*fact;
				this.regulatorWeights.put(splt.regulator, w);
			}
     }
     
     /**
      * Compute regulator weights for randomly assigned regulators
      * 
      * @author eric  
      */
     public void setRegulatorWeightsRandom (){
       this.regulatorWeightsRandom = new HashMap<Gene, Double>();
       double fact = (double)this.leafDistribution.condSet.size()/(double)this.module.moduleNetwork.data[0].length;
       for (Split splt : this.testSplitsRandom){
        if (!this.regulatorWeightsRandom.containsKey(splt.regulator))
          this.regulatorWeightsRandom.put(splt.regulator, 0.0);
        double w = this.regulatorWeightsRandom.get(splt.regulator)
          + Math.exp(splt.regulatorScore/this.leafDistribution.condSet.size())*fact;
        this.regulatorWeightsRandom.put(splt.regulator, w);
      }
     }

     
     /**
     * Finds rightmost leaf node below this one. Used to attach gene labels in 
     * {@link Module#draw(String)}.
     * 
     * @author stmae
     */
     public TreeNode rightMostNode () {
		if (this.nodeStatus.equals("leaf")){ // no tree
			return this;
		}	
		else {
			TreeNode tempNode = this;
		 	while(tempNode.rightChild.nodeStatus.equals("internal"))
		 		tempNode = tempNode.rightChild;
			return tempNode;
		}
    }
    
    /**
     * Prints the statistics of all leaves below this node.
     * 
     * @author tomic
     */
     public void printStatistics(){
        if (this.nodeStatus.equals(LEAF)) {
            System.out.println(this.leafDistribution.statistics[0] + "\t"
                    + this.leafDistribution.statistics[1] + "\t"
                    + this.leafDistribution.statistics[2] + "\t"
                    + this.leafDistribution.condSet.size());
        } else {
            this.leftChild.printStatistics();
            this.rightChild.printStatistics();
        }
    }
    
    /**
     * Print the temporary statistics of all leaves below this node.
     * 
     * @author tomic
     */
     public void printStatisticsTmp(){
        if (this.nodeStatus.equals(LEAF)) {
            System.out.println(this.leafDistribution.statisticsTmp[0] + "\t"
                    + this.leafDistribution.statisticsTmp[1] + "\t"
                    + this.leafDistribution.statisticsTmp[2] );
        } else {
            this.leftChild.printStatisticsTmp();
            this.rightChild.printStatisticsTmp();
        }
    }
    
    /**
     * Check if this node is higher in the tree than a given node.
     * 
     * @param node
     * @return <code>true</code> upon success.
     */
    public boolean isParent(TreeNode node) {
        if (this.nodeStatus == LEAF)
            return false;
        else if (this.leftChild == node | this.rightChild == node)
            return true;
        else if (this.leftChild.isParent(node) | this.rightChild.isParent(node))
            return true;
        return false;
    }
    

    /**
     * Gather all leaves below this node.
     * 
     * @param leaves should be empty initially, will be filled recursively.
     */
    public void gatherLeaves(HashSet<TreeNode> leaves) {
        if (this.nodeStatus.equals(LEAF)) {
            leaves.add(this);
        } else {
            this.leftChild.gatherLeaves(leaves);
            this.rightChild.gatherLeaves(leaves);
        }
    }
    

    /**
     * Gather all leaves below this node.
     * 
     * @param leafList should be empty initially, will be filled recursively.
     */
    public void gatherLeafList(ArrayList<TreeNode> leafList) {
        if (this.nodeStatus.equals(LEAF)) {
            leafList.add(this);
        } else {
            this.leftChild.gatherLeafList(leafList);
            this.rightChild.gatherLeafList(leafList);
        }
    }
    
    /**
     * Gather all internal nodes below and including this node.
     * 
     * @param internals should be empty initially, will be filled recursively.
     */
    public void gatherInternals(ArrayList<TreeNode> internals) {
        if (this.nodeStatus.equals(INTERNAL)) {
            internals.add(this);
            this.leftChild.gatherInternals(internals);
            this.rightChild.gatherInternals(internals);
        }
    }
    
    /**
     * Return list of leaf nodes below and including this node.
     * 
     * @author tomic
     */
    public ArrayList<TreeNode> getLeafNodes (){
    	ArrayList<TreeNode> nodes = new ArrayList<TreeNode>();
    	this.gatherLeafList(nodes);
    	return nodes;
    }
    
    /**
     * Return list of internal nodes below and including this node.
     * 
     * @author tomic
     */
    public ArrayList<TreeNode> getInternalNodes (){
    	ArrayList<TreeNode> nodes = new ArrayList<TreeNode>();
    	this.gatherInternals(nodes);
    	return nodes;
    }
    
    /**
     * Converts this node to a leaf
     *
     */
    public void toLeaf(){
    	if (this.nodeStatus.equals(INTERNAL)){
    		this.nodeStatus = LEAF;
    		this.leftChild = null;
    		this.rightChild = null;
    	}
    }
    
    /**
     * Returns a partition of the conditions from the condition sets
     * of the leaf nodes of this node.
     * @return
     * 
     * @author tomic
     */
    public ArrayList<ArrayList<Integer>> condSets (){
    	ArrayList<ArrayList<Integer>> list = new ArrayList<ArrayList<Integer>>();
    	for (TreeNode leaf : this.getLeafNodes())
    		list.add(leaf.leafDistribution.condSet);
    	return list;
    }
    
    /**
     * For each leaf below this node: copy statistics to statisticsTmp
     * 
     *  @author tomic
     */
    public void copyStatisticsTmp() {
        if (this.nodeStatus.equals(LEAF)) {
            for (int i = 0; i < 4; i++) //revamp
                this.leafDistribution.statisticsTmp[i] = this.leafDistribution.statistics[i];
        } else {
            this.leftChild.copyStatisticsTmp();
            this.rightChild.copyStatisticsTmp();
        }
    }
    
    /**
     * For each leaf below this node: copy statisticsTmp to statistics
     * 
     *  @author tomic
     */
    public void copyStatistics() {
        if (this.nodeStatus.equals(LEAF)) {
            for (int i = 0; i < 4; i++) //revamp
                this.leafDistribution.statistics[i] = this.leafDistribution.statisticsTmp[i];
        } else {
            this.leftChild.copyStatistics();
            this.rightChild.copyStatistics();
        }
    }
    
    /**
     * Computes the Bayesian score sum of all leaves below this node. 
     * Does not perform any score computation, but recursively sums the content of the 
     * {@link LeafDistribution#score} fields.
     * 
     * @return sum of scores.
     * @author tomic
     */
    public double bayesianScore() {
        double score;
        if (this.nodeStatus.equals(LEAF)) {
            this.leafDistribution.bayesianScore();
            //revamp debugging
            //System.err.println("BayesianScore mod num " + this.leafDistribution.score  +
            //		" * " + this.leafDistribution.condSet.get(0) +
            //		"/" + this.leafDistribution.condSet.size() +
            //		" : " + this.module.moduleNetwork.condition_weight[ this.leafDistribution.condSet.get(0) ] +  "   ");
            if (Globals.unweighted) {
                score = this.leafDistribution.score;
            } else {
                score = this.leafDistribution.score * this.module.moduleNetwork.condition_weight[ this.leafDistribution.condSet.get(0) ]; //revamp
            }
        } else {
            score = this.leftChild.bayesianScore() + this.rightChild.bayesianScore();
        }
        return score;
    }
    
    /**
     * Computes the Bayesian score sum of all leaves below this node based on the
     * temporary statistics. Does not perform any score computation, but recursively 
     * sums the content of the {@link LeafDistribution#scoreTmp} fields.
     * 
     * @return sum of scores.
     * @author tomic
     */
    public double bayesianScoreTmp() {
        double score;
        if (this.nodeStatus.equals(LEAF)) {
            this.leafDistribution.bayesianScoreTmp();
            if (Globals.unweighted) {
                score = this.leafDistribution.scoreTmp;
            } else {
                score = this.leafDistribution.scoreTmp * this.module.moduleNetwork.condition_weight[ this.leafDistribution.condSet.get(0) ]; //revamp
            }
        } else {
            score = this.leftChild.bayesianScoreTmp() + this.rightChild.bayesianScoreTmp();
        }
        return score;
    }
    
    
    /**
     * <p>Computes the log of the partition sum. The partition sum is the sum over all
     * partitions consistent with the tree of the exponential of the score of the partition. 
     * Or, if a node is a leaf, the partition sum is defined as the exponential of the 
     * leaf score, <em>Z=e<sup>score</sup></em>. 
     * And if a node is an internal node, it is defined recursively as <em>Z=e<sup>score</sup> + 
     * Z<sub>1</sub> Z<sub>2</sub></em>, with <em>score</em> the score of all experiments
     * below this node clustered into 1 set, and <em>Z<sub>1</sub></em> and 
     * <em>Z<sub>2</sub></em> the partition sums of the children of the node.</p>
     * <p>This is used in {@link Module#hierarchicalClustering(List, boolean)} if the merge
     * score based on Bayesian Hierarchical Clustering is chosen.
     * 
     * @return log of partition sum.
     * @author tomic
     * 
     */
    public double logPartSum(){
        double lp;
        if (this.nodeStatus.equals(LEAF)){
            lp = this.leafDistribution.score;
        } else {
            lp = this.leafDistribution.score +
            	Math.log(1 + Math.exp(-this.leafDistribution.score +
            			this.leftChild.logPartSum() + this.rightChild.logPartSum()));
        }
        return lp;
    }
    
 
    /**
     * Sets the merge score for {@link Module#hierarchicalClustering(List, boolean)}.
     * This can be the difference between the score of all experiments below this node 
     * clustered into 1 set, and the sum of scores of all experiments below this node 
     * clustered into 2 sets corresponding to the condition sets of the 2 children of
     * this node. It can also be based on Bayesian Hierarchical Clustering. Then it is 
     * defined as the difference between the score of all experiments below this node 
     * clustered into 1 set and the sum of the log of the partition sums of the 2 children 
     * of this node. The results is stored in the field {@link TreeNode#mergeScore}.
     * 
     * @param useBHCscore which type of merge score to compute.
     * @author tomic
     */
    public void setMergeScore(boolean useBHCscore){
    	if (useBHCscore){
    		this.mergeScore =  this.leafDistribution.score
        		- this.leftChild.logPartSum() - this.rightChild.logPartSum();
    	} else {
    		this.mergeScore = this.leafDistribution.score
    			- this.leftChild.leafDistribution.score - this.rightChild.leafDistribution.score;
    	}
    }
    
    /**
     * Trims a hierarchical tree at a certain level.
     * 
     * @author tomic
     */
    public void testLevel(int level){
        if (this.nodeStatus.equals("internal")) {
            if (this.treeLevel() <= level) {
                this.leftChild.testLevel(level);
                this.rightChild.testLevel(level);
            } else {
            	this.nodeStatus = "leaf";
                this.leftChild = null;
                this.rightChild = null;
            }
        }
    }
    
    /**
     * Trims a hierarchical tree keeping only nodes with positive merge score.
     * 
     * @author tomic
     */
    public void testMergeScore(){
        if (this.nodeStatus.equals("internal")) {
            if (this.mergeScore < 0) {
                this.leftChild.testMergeScore();
                this.rightChild.testMergeScore();
            } else {
                 this.nodeStatus = "leaf";
                this.leftChild = null;
                this.rightChild = null;
            }
        }
    }
    
     
    /**
     * Trims a hierarchical tree keeping only nodes with a score difference higher 
     * than some value.
     * 
     * @param scoregain threshold value.
     * @author tomic
     */
     public void testScore(double scoregain){
        double diff;
        //System.out.println("testScore of module " + this.leafDistribution.statistics[3] );
        if (this.nodeStatus.equals("internal")) {
            diff = (-this.leafDistribution.score
                    + this.leftChild.leafDistribution.score + this.rightChild.leafDistribution.score)/this.module.moduleNetwork.geneSet.size();
             if (diff > scoregain) {
                this.leftChild.testScore(scoregain);
                this.rightChild.testScore(scoregain);
            } else {
                this.nodeStatus = "leaf";
                this.leftChild = null;
                this.rightChild = null;
            }
        }
    }
    
     /**
      * Trims a hierarchical tree keeping only nodes whose regulator
      * has high enough average prediction probability.
      * 
      * @param minScore
      */
    public void testRegulatorScore(double minProb){
    	if (this.nodeStatus.equals("internal")){
    		double prob = 0.0;
    		if (this.testSplits.isEmpty()){
    			prob = Math.exp(this.regulationSplit.regulatorScore/(double)this.leafDistribution.condSet.size());
    		} else {
    			for (Split splt : this.testSplits){
    				prob += Math.exp(splt.regulatorScore/(double)this.leafDistribution.condSet.size());
    			}
    			prob = prob/(double)this.testSplits.size();
    		}
    		if (prob > minProb){
    			this.leftChild.testRegulatorScore(minProb);
    			this.rightChild.testRegulatorScore(minProb);
    		} else {
                this.nodeStatus = "leaf";
                this.leftChild = null;
                this.rightChild = null;
            }
    	}
    }
     
    /**
     * Trims a hierarchical tree such that the partition of the experiments has the
     * highest possible score. See also {@link TreeNode#setPartitions(double)}.
     * 
     * @author tomic
     */
    public void testPartition(){
        if (this.nodeStatus.equals("internal")){
            if (this.partition == true){
                this.leftChild.testPartition();
                this.rightChild.testPartition();
            } else {
                this.nodeStatus = "leaf";
                this.leftChild = null;
                this.rightChild = null;
            }
        }
    }
    
    /**
     * Gets all nodes at a certain depth below this node.
     * 
     * @param depth
     * @param trees initially empty, will be filled recursively.
     */
    public void getTreesAtDepth(int depth, ArrayList<TreeNode> trees){
        if (this.treeDepth() <= depth){
            trees.add(this);
        } else {
            this.leftChild.getTreesAtDepth(depth, trees);
            this.rightChild.getTreesAtDepth(depth, trees);
        }
    }
    

    /**
     * This sets for all internal nodes below and including this one a boolean 
     * {@link TreeNode#partition} as follows. For the internal nodes right above the 
     * leaf nodes, it is <code>true</code> if the node has a negative score 
     * difference with its leaf children, <code>false</code> otherwise. Then it moves
     * up in the tree and sets {@link TreeNode#partition} to <code>true</code> if the 
     * score of the node (=all experiments below it in 1 cluster) is higher than the
     * highest scoring partition below that node, <code>false</code> otherwise.
     * The highest scoring partition can be found because if all nodes below a node have 
     * their  {@link TreeNode#partition} field set, we only have to sum the scores of all 
     * nodes with {@link TreeNode#partition}<code>=false</code> below the node.
     *
     *  @param scoregain threshold for accepting score differences.
     *  @author tomic
     */
    public void setPartitions(double scoregain){
        HashSet<TreeNode> leaves = new HashSet<TreeNode>();
        this.gatherLeaves(leaves);
        // new depth
        int depth = this.treeDepth();
        // internal nodes
        ArrayList<TreeNode> internals = new ArrayList<TreeNode>();
        this.gatherInternals(internals);
        // set partitions at the bottom
        for (TreeNode node : internals){
            if (node.treeDepth() == 1){
                if (-node.leafDistribution.score
                        + node.leftChild.leafDistribution.score + node.rightChild.leafDistribution.score
                        > scoregain){
                    node.partition = true;
                } else {
                    node.partition = false;
                }
            }
        }
        // move up to root
        for (int l = 2; l<=depth; l++){
            // loop over nodes at level l
            for (TreeNode node : internals)
                if (node.treeDepth() == l){
                double sum = node.leftChild.bestPartitionScore()
                	+ node.rightChild.bestPartitionScore();
                if (sum - node.leafDistribution.score > scoregain)
                    node.partition = true;
                else
                    node.partition = false;
                }
        }
    }
    
    /**
     * Computes the score of the best partition below this node by summing the scores
     * of all nodes below this one with {@link TreeNode#partition}<code>=false</code>.
     * 
     * @return score.
     * @author tomic
     * 
     */
    public double bestPartitionScore(){
        double sum;
        if (this.nodeStatus.equals("leaf"))
            sum = this.leafDistribution.score;
        else {
            if (this.partition == false)
                sum = this.leafDistribution.score;
            else{
                sum = this.leftChild.bestPartitionScore() + this.rightChild.bestPartitionScore();
            }
        }
        return sum;
    }
    

    /**
     * Finds the regulator with minimal Bayesian posterior assignment probability
     * for an internal node.
     * 
     * @param regulatorSet set of regulators from which to choose.
     * @author tomic
     */
    public void findRegulatorBayes(ArrayList<Gene> regulatorSet, double beta){
    	if (this.nodeStatus.equals("internal")){
    		double [][] data = this.module.moduleNetwork.data;
    		this.regulationSplit = new Split(this);
            Split bestSplit = new Split(this);
            Split testSplit;
            double bestScore = Double.NEGATIVE_INFINITY;
            for (Gene reg : regulatorSet){
            	// limit split values to values in condition set
            	reg.setNormData(data);
            	reg.possibleSplitValues = new ArrayList<Double>();
            	for (int m : this.leafDistribution.condSet)
            		if (!Double.isNaN(reg.normData[m]))
            			reg.possibleSplitValues.add(reg.normData[m]);
                for (double sv : reg.possibleSplitValues){
                	testSplit = new Split(this, reg, sv);
                	// compute sign of split, quit if it's no good
                	testSplit.computeSign();
//                	testSplit.beta = beta;
                	if (testSplit.sign == 0){
                		continue;
                	}
                	// find optimal beta
                	RootFinder rf = new RootFinder(testSplit, 0.0, beta, 1E-5);
                	// check if betaMax is large enough, then find root
            		try {
            			rf.expandXMax();
            			try{
                    		rf.bisect();
                    		testSplit.beta = rf.xRoot;
                    	} catch (Exception e){
                    		System.out.println(e);
                    		continue;
                    	}
            		} catch (Exception f){
            			// these are some exceptional cases where split.evaluate()
            			// is always zero, don't consider them
            			continue;
            		}
                	testSplit.computeRegulatorScore();
                	if (testSplit.regulatorScore > bestScore){
                        bestScore = testSplit.regulatorScore;
                        bestSplit = testSplit;
                        bestSplit.computeScoreGain();
                	} 
                }
            }
            this.regulationSplit = bestSplit;
    	}
    }
    
    /**
     * Finds a number of regulators with low Bayesian assignment score stochastically. Regulators are
     * stored in {@link TreeNode#testSplits}. In addition a same number of random 
     * assignments (with uniform probability) are stored in {@link TreeNode#testSplitsRandom} 
     * for statistical comparisons.
     * 
     * @param numRegAssign number of regulators to assign.
     * @param regulators set of regulators from which to choose.
     * @param beta parameter in Bayesian probability of each regulator.
     * 
     * @author tomic
     */
    public void findRegulatorBayesStoch(int numRegAssign, ArrayList<Gene> regulators, double beta){
        if (this.nodeStatus.equals("internal")){
        	this.testSplits = new ArrayList<Split>(numRegAssign);
        	this.testSplitsRandom = new ArrayList<Split>(numRegAssign);
        	ArrayList<Split> regulatorSplits = this.computeRegulatorSplitsBayes(regulators, beta);
        	ArrayList<Double> regulatorSplitProbs = new ArrayList<Double>(regulatorSplits.size());
          double maxRegulatorScore = Double.NEGATIVE_INFINITY;
          for (Split splt : regulatorSplits){
            if (maxRegulatorScore < splt.regulatorScore) {
              maxRegulatorScore = splt.regulatorScore;
            }
        	}
          // fill regulatorSplitProbs with weights, compute normalization along the way
          double probnorm = 0.0;
        	for (Split splt : regulatorSplits){
        		double weight = Math.exp(splt.regulatorScore - maxRegulatorScore);
        		probnorm += weight;
        		regulatorSplitProbs.add(weight);
        	}
        	// renormalize each value to get probabilities
        	for (int k=0; k<regulatorSplitProbs.size(); k++){
        		regulatorSplitProbs.set(k, regulatorSplitProbs.get(k)/probnorm);
        	}
        	
        	// set regulator weights 
//    		for (int k=0; k<regulatorSplits.size(); k++){
//    			if (!this.module.regulatorWeights.keySet().contains(regulatorSplits.get(k).regulator))
//    				this.module.regulatorWeights.put(regulatorSplits.get(k).regulator, 0.0);
//    			double w = this.module.regulatorWeights.get(regulatorSplits.get(k).regulator);
//    			if (regulatorSplitProbs.get(k) > w)
//    				this.module.regulatorWeights.put(regulatorSplits.get(k).regulator, regulatorSplitProbs.get(k));
//    		}
        	
            for (int k=0; k<numRegAssign; k++){
            	// roll the dice and choose a split according to SplitProbs
            	Random random = new Random();
            	double rand = random.nextDouble();
            	int splitind = 0;
            	double partialsum = regulatorSplitProbs.get(0);
            	while (rand > partialsum){
            		splitind += 1;
            		partialsum += regulatorSplitProbs.get(splitind);
            	}
            	// assign split corresponding to splitind 
            	this.regulationSplit = regulatorSplits.get(splitind);
            	this.testSplits.add(regulatorSplits.get(splitind));
            	// also pick a completely random one
            	int ri = random.nextInt(regulatorSplits.size());
            	this.testSplitsRandom.add(regulatorSplits.get(ri));
            }
        }
    }
    
    public void findAverageRegulatorBayes(ArrayList<Gene> regulators, double beta){
    	if (this.nodeStatus.equals("internal")){
    		System.out.println("Computing all regulator - split value scores ...");
    		ArrayList<Split> regulatorSplits = this.computeRegulatorSplitsBayes(regulators, beta);
    		System.out.println("... done");
        	// Create "average" regulator
        	Gene avgReg = new Gene(regulatorSplits, this.module.moduleNetwork.data);
        	// compute average split value
        	double zAvg = 0.0, norm = 0.0;
        	for (Split splt : regulatorSplits){
        		zAvg += splt.splitValue * Math.exp(splt.regulatorScore);
        		norm += Math.exp(splt.regulatorScore);
        	}
        	zAvg = zAvg/norm;
        	this.regulationSplit = new Split(this, avgReg, zAvg);
        	this.regulationSplit.computeBeta();
        	this.regulationSplit.computeRegulatorScore();
    	}
    }
    
    /**
     * Computes and stores the Bayesian score for each regulator - split value pair
     * 
     * @param regulators regulators for which to compute the score.
     * @return list of splits, one for each regulator - split value pair.
     * 
     * @author tomic
     */
    public ArrayList<Split> computeRegulatorSplitsBayes (ArrayList<Gene> regulators, double beta){
    	ArrayList<Split> regulatorSplits = new ArrayList<Split>();
    	if (this.nodeStatus.equals("internal")){
    		double [][] data = this.module.moduleNetwork.data;
     		this.regulatorWeights = new HashMap<Gene,Double>();
    		for (Gene reg : regulators){
            	// limit split values to values in condition set
    			reg.setNormData(data);
            	reg.possibleSplitValues = new ArrayList<Double>();
            	for (int m : this.leafDistribution.condSet)
            		if (!Double.isNaN(reg.normData[m]))
            			reg.possibleSplitValues.add(reg.normData[m]);
     			for (double sv : reg.possibleSplitValues){
    				Split splt = new Split(this, reg, sv);
    				// compute sign of split, go to next split value if it's no good
                	splt.computeSign();
                	if (splt.sign == 0)
                		continue;
                	
                	// find optimal beta
                	RootFinder rf = new RootFinder(splt, 0.0, 20, 1E-5);
                	// check if betaMax is large enough, then find root
            		try {
            			rf.expandXMax();
            			try{
                    		rf.bisect();
                    		splt.beta = rf.xRoot;
                    	} catch (Exception e){
                    		System.out.println(e);
                    		continue;
                    	}
            		} catch (Exception f){
            			// these are some exceptional cases where split.evaluate()
            			// is always zero, don't consider them
            			continue;
            		}
            		splt.computeRegulatorScore();
            		regulatorSplits.add(splt);
    			}
    		}
    		if (regulatorSplits.isEmpty()){ // this is unusual but seems to happen sometimes ... add dummy to prevent null pointers
    			this.toLeaf();
    			regulatorSplits.add(new Split());
    		}
    	}
    	return regulatorSplits;
    }
    
    /**
     * Finds the regulator with minimal entropy for an internal node.
     * The entropy is the conditional entropy for an experiment to belong to
     * the left or right child, given that the regulator is expressed
     * above or below a split value in that experiment.
     * 
     * @param regulatorSet  set of regulators from which to choose.
     * @author tomic
     */
    public void findRegulator(HashSet<Gene> regulatorSet){
    	if (this.nodeStatus.equals("internal")){
    		double [][] data = this.module.moduleNetwork.data;
            this.regulationSplit = new Split();
            Split bestSplit = new Split();
            double bestEnt = 1.0; // upper bound for entropy
            for (Gene reg : regulatorSet)
                for (double split : reg.possibleSplitValues){
                    ArrayList<Integer> lowConds = new ArrayList<Integer>();
                    ArrayList<Integer> highConds = new ArrayList<Integer>();
                    for (int m: this.leafDistribution.condSet)
                        if (data[reg.number][m] <= split)
                            lowConds.add(m);
                        else
                            highConds.add(m);
                    double p1 = 0.0; // fraction of exp in lowConds \cap this.left.condset
                    double p2 = 0.0; // fraction of exp in lowConds \cap this.right.condset
                    for (int m: this.leftChild.leafDistribution.condSet)
                        if (lowConds.contains(m))
                            p1 += 1.0;
                        else
                            p2 += 1.0;
                    p1 = p1/lowConds.size();
                    p2 = p2/highConds.size();
                    double entropy = this.ent(p1) * lowConds.size()/this.leafDistribution.condSet.size()
                    + this.ent(p2)* highConds.size()/this.leafDistribution.condSet.size();
                    if (entropy < bestEnt){
                        bestEnt = entropy;
                        bestSplit.regulator = reg;
                        bestSplit.splitValue = split;
                        bestSplit.regulatorScore = entropy;
                    } 
                }
                this.regulationSplit = bestSplit;
        }
    }
    
    /**
     * Finds a number of low-entropy regulators stochastically. Regulators are chosen
     * with a probability proportional to <em>e<sup>-&beta;*entropy<sup>2</sup></sup></em> 
     * and stored in {@link TreeNode#testSplits}. In addition a same number of random 
     * assignments (with uniform probability) are stored in {@link TreeNode#testSplitsRandom} 
     * for statistical comparisons.
     * 
     * @param numRegAssign number of regulators to assign.
     * @param regulators set of regulators from which to choose.
     * @param beta parameter in probability weight of each regulator.
     * 
     * @author tomic
     */
    public void findRegulatorStoch(int numRegAssign, HashSet<Gene> regulators, double beta){
        if (this.nodeStatus.equals("internal")){
        	this.testSplits = new ArrayList<Split>(numRegAssign);
        	this.testSplitsRandom = new ArrayList<Split>(numRegAssign);
        	ArrayList<Split> regulatorSplits = this.computeRegulatorSplits(regulators);
        	ArrayList<Double> regulatorSplitProbs = new ArrayList<Double>(regulatorSplits.size());
        	// fill regulatorSplitProbs with weights, compute norm along the way
        	double probnorm = 0.0;
        	for (Split splt : regulatorSplits){
        		double weight = Math.exp(-beta*Math.pow(splt.regulatorScore,2));
        		probnorm += weight;
        		regulatorSplitProbs.add(weight);
        	}
        	// renormalize each value to get probabilities
        	for (int k=0; k<regulatorSplitProbs.size(); k++){
        		regulatorSplitProbs.set(k, regulatorSplitProbs.get(k)/probnorm);
        	}
        	
            for (int k=0; k<numRegAssign; k++){
            	// roll the dice and choose a split according to SplitProbs
            	Random random = new Random();
            	double rand = random.nextDouble();
            	int splitind = 0;
            	double partialsum = regulatorSplitProbs.get(0);
            	while (rand > partialsum){
            		splitind += 1;
            		partialsum += regulatorSplitProbs.get(splitind);
            	}
            	// assign split corresponding to splitind 
            	this.regulationSplit = regulatorSplits.get(splitind);
            	this.testSplits.add(regulatorSplits.get(splitind));
            	// also pick a completely random one
            	int ri = random.nextInt(regulatorSplits.size());
            	this.testSplitsRandom.add(regulatorSplits.get(ri));
            }
        }
    }
    
    /**
     * Computes the entropy for each (regulator, split value) pair.
     * 
     * @param regulators regulators for which to compute the entropy.
     * @return list of splits, one for each (regulator, split value) pair.
     * 
     * @author tomic
     */
    public ArrayList<Split> computeRegulatorSplits (HashSet<Gene> regulators){
    	ArrayList<Split> regulatorSplits = new ArrayList<Split>();
    	if (this.nodeStatus.equals("internal")){
    		double [][] data = this.module.moduleNetwork.data;
    		// set regulator - split value pairs
    		for (Gene reg : regulators)
    			for (double sv : reg.possibleSplitValues)
    				regulatorSplits.add(new Split(this, reg, sv));
    		// compute entropy for each regulator - split value pair
    		for (Split split : regulatorSplits) {
    			ArrayList<Integer> lowConds = new ArrayList<Integer>();
    			ArrayList<Integer> highConds = new ArrayList<Integer>();
    			for (int m: this.leafDistribution.condSet)
    				if (data[split.regulator.number][m] <= split.splitValue)
    					lowConds.add(m);
    				else
    					highConds.add(m);
    			double p1 = 0.0; // fraction of exp in lowConds \cap this.left.condset
    			double p2 = 0.0; // fraction of exp in lowConds \cap this.right.condset
    			for (int m: this.leftChild.leafDistribution.condSet)
    				if (lowConds.contains(m))
    					p1 += 1.0;
    				else
    					p2 += 1.0;
    			p1 = p1/lowConds.size();
    			p2 = p2/highConds.size();
    			double entropy = this.ent(p1) * lowConds.size()/this.leafDistribution.condSet.size()
    				+ this.ent(p2)* highConds.size()/this.leafDistribution.condSet.size();
    			split.regulatorScore = entropy;
    		}
    	}
    	return regulatorSplits;
    }
    
     
    /**
     * The binary entropy function.
     * 
     * @param p value between 0 and 1.
     * @return -p*log(p) -(1-p)*log(1-p)
     * 
     * @author tomic
     */
    private double ent(double p){
        double plogp = 0.0;
        if (p > 0.0 && p < 1.0)
            plogp = -p*Math.log(p)-(1-p)*Math.log(1-p);
        return plogp;
    }
    
    
    /**
     * Tries to split a leaf node using a certain regulator and split value. Used in
     * top down search of regulation programs.
     * 
     * @param testGene candidate regulator
     * @param splitValue
     * @return increase in Bayesian score if this split were made.
     * 
     * @author tomic
     */
    public double trySplit(Gene testGene, double splitValue) {
        if (this.nodeStatus == INTERNAL) {
            return 0.0;
        } else {
        	double [][] data = this.module.moduleNetwork.data;
            // use ficticious child nodes
            TreeNode leftChild = new TreeNode(this.module, this.leafDistribution.normalGammaPrior);
            TreeNode rightChild = new TreeNode(this.module, this.leafDistribution.normalGammaPrior);
            // compute statistics for fictitious leafs
            TreeNode node;
            // loop over all conditions that end at current node
            for (Integer m : this.leafDistribution.condSet) {
                // check new split value
                if (data[testGene.number][m] <= splitValue)
                    node = leftChild;
                else
                    node = rightChild;
                // add m to conditions of new node
                node.leafDistribution.condSet.add(m);
                // add expression values for this condition and for all genes in
                // module to present leaf statistics
                for (Gene gene : this.module.genes) {
                    if (Double.isNaN(data[gene.number][m]) == false) {
                        node.leafDistribution.statistics[0]++;
                        node.leafDistribution.statistics[1] = node.leafDistribution.statistics[1]
                                + data[gene.number][m];
                        node.leafDistribution.statistics[2] = node.leafDistribution.statistics[2]
                                + Math.pow(data[gene.number][m], 2);
                    }
                }
            }
            // scores for new leafs
            leftChild.leafDistribution.bayesianScore();
            rightChild.leafDistribution.bayesianScore();
            // calculate projected score : subtract old score and add 2 new scores
            double scoreGain = leftChild.leafDistribution.score
                    + rightChild.leafDistribution.score
                    - this.leafDistribution.score;
            
            return scoreGain;
        }
    }
    
    /**
     * Splits a leaf. This turns a leaf into an internal node, creating 2 child nodes and 
     * set their statistics and score. Used in top down search of regulation programs.
     * 
     * @param testGene candidate regulator
     * @param splitValue
     * 
     * @author tomic
     */
    public void split(Gene testGene, double splitValue) {
    	this.regulationSplit = new Split(this, testGene, splitValue);
        this.leftChild = new TreeNode(this);
        this.rightChild = new TreeNode(this);
        this.nodeStatus = INTERNAL;
        double[][] data = this.module.moduleNetwork.data;
        
        // compute statistics for new leafs and update score
        TreeNode node;
        // loop over all conditions that end at this
        for (Integer m : this.leafDistribution.condSet) {
            // check new split value
            if (data[testGene.number][m] <= splitValue)
                node = this.leftChild;
            else
                node = this.rightChild;
            // add m to conditions of new node
            node.leafDistribution.condSet.add(m);
            // add expression values for this condition and for all genes in
            // Genes to present leaf statistics
            for (Gene gene : this.module.genes) {
                if (Double.isNaN(data[gene.number][m]) == false) {
                    node.leafDistribution.statistics[0]++;
                    node.leafDistribution.statistics[1] += data[gene.number][m];
                    node.leafDistribution.statistics[2] += Math.pow(data[gene.number][m], 2);
                }
            }
        }
        // scores for new leafs
        this.leftChild.leafDistribution.bayesianScore();
        this.rightChild.leafDistribution.bayesianScore();
    }
    
    
    /**
     * Computes the statistics for all leaves below this node.
     * 
     * @author tomic
     */
    public void statistics() {
    	double [][] data = this.module.moduleNetwork.data;
        ArrayList<TreeNode> leaves = new ArrayList<TreeNode>();
        this.gatherLeafList(leaves);
        for (TreeNode leaf : leaves)
            for (int m : leaf.leafDistribution.condSet) {
                for (Gene gene : this.module.genes)
                    if (Double.isNaN(data[gene.number][m]) == false) {
                    	leaf.leafDistribution.statistics[0]++;
                    	leaf.leafDistribution.statistics[1] = leaf.leafDistribution.statistics[1]
                            + data[gene.number][m];
                    	leaf.leafDistribution.statistics[2] = leaf.leafDistribution.statistics[2]
                            + Math.pow(data[gene.number][m], 2);
                    }
            //System.err.println("statistics  " + leaf.leafDistribution.statistics[3]);    
            leaf.leafDistribution.statistics[3] = this.module.moduleNetwork.condition_weight[m];
        }
    }

    
    /**
     * Computes the statistics as if a certain gene is added to the module of this node. 
     * The result is stored in {@link LeafDistribution#statisticsTmp}.
     * 
     * @param gene gene to be added.
     * 
     * @author tomic
     */
    public void statisticsAddGeneTmp(Gene gene) {
    	double [][] data = this.module.moduleNetwork.data;
        ArrayList<TreeNode> leaves = new ArrayList<TreeNode>();
        this.gatherLeafList(leaves);
        for (TreeNode leaf : leaves)
            for (int m : leaf.leafDistribution.condSet) {
                if (Double.isNaN(data[gene.number][m]) == false) {
                	leaf.leafDistribution.statisticsTmp[0]++;
                	leaf.leafDistribution.statisticsTmp[1] = leaf.leafDistribution.statisticsTmp[1] + data[gene.number][m];
                	leaf.leafDistribution.statisticsTmp[2] = leaf.leafDistribution.statisticsTmp[2] + Math.pow(data[gene.number][m], 2);
                }
                leaf.leafDistribution.statistics[3] = this.module.moduleNetwork.condition_weight[m]; //revamp
            }    
                
    }
    

    /**
     * Computes the statistics as if a certain gene is removed from the module of this node. 
     * The result is stored in {@link LeafDistribution#statisticsTmp}.
     * 
     * @param gene gene to be removed.
     * 
     * @author tomic
     */
    public void statisticsRemoveGeneTmp(Gene gene) {
    	double[][] data = this.module.moduleNetwork.data;
        ArrayList<TreeNode> leaves = new ArrayList<TreeNode>();
        this.gatherLeafList(leaves);
        for (TreeNode leaf : leaves)
            for (int m : leaf.leafDistribution.condSet) {
                if (Double.isNaN(data[gene.number][m]) == false) {
                	leaf.leafDistribution.statisticsTmp[0]--;
                	leaf.leafDistribution.statisticsTmp[1] = leaf.leafDistribution.statisticsTmp[1] - data[gene.number][m];
                	leaf.leafDistribution.statisticsTmp[2] = leaf.leafDistribution.statisticsTmp[2] - Math.pow(data[gene.number][m], 2);
                }
                leaf.leafDistribution.statistics[3] = this.module.moduleNetwork.condition_weight[m]; //revamp
            }    
    }
    
    
    
    
    /** compute the statistics for all leaves below this node from scratch
     * 
     * @param data
     * @param Genes
     * @param Experiments
     * 
     * @author tomic
     */
    public void statistics(double[][] data, HashSet<Gene> Genes, HashSet<Experiment> Experiments) {
        for (Experiment exp : Experiments) {
            TreeNode pos = this; // start at this node
            while (pos.nodeStatus.equals(INTERNAL)) { // go down in tree
                if (data[pos.regulationSplit.regulator.number][exp.number] <= pos.regulationSplit.splitValue)
                    pos = pos.leftChild;
                else
                    pos = pos.rightChild;
            }
            // add m to the set of conditions for this leaf
            pos.leafDistribution.condSet.add(exp.number);
            // add expression values for this condition and for all gene-exps in
            // module to present leaf statistics
            for (Gene gene : Genes) {
                if (Double.isNaN(data[gene.number][exp.number]) == false) {
                    pos.leafDistribution.statistics[0]++;
                    pos.leafDistribution.statistics[1] = pos.leafDistribution.statistics[1]
                            + data[gene.number][exp.number];
                    pos.leafDistribution.statistics[2] = pos.leafDistribution.statistics[2]
                            + Math.pow(data[gene.number][exp.number], 2);
                }
            }
        }
    }
    
    /**
     * update the statistics when gene is added to the module, and stores the
     * result in the "temporary" leaf statistics variable.
     * adding a gene affects all leaf statistics, so bayesianScore has to be
     * recomputed from scratch afterwards.
     * 
     * @author tomic
     */ 
    public void statisticsAddGeneTmp(double[][] data, HashSet<Experiment> Experiments, Gene gene) {
        for (Experiment exp : Experiments) {
            TreeNode pos = this; // start at this
            while (pos.nodeStatus.equals(INTERNAL)) { // go down in tree
                if (data[pos.regulationSplit.regulator.number][exp.number] <= pos.regulationSplit.splitValue)
                    pos = pos.leftChild;
                else
                    pos = pos.rightChild;
            }
            // add expression values for this condition and this gene to present
            // leaf statistics
            if (Double.isNaN(data[gene.number][exp.number]) == false) {
                pos.leafDistribution.statisticsTmp[0]++;
                pos.leafDistribution.statisticsTmp[1] = pos.leafDistribution.statisticsTmp[1]
                        + data[gene.number][exp.number];
                pos.leafDistribution.statisticsTmp[2] = pos.leafDistribution.statisticsTmp[2]
                        + Math.pow(data[gene.number][exp.number], 2);
            }
        }
    }
    
    /** update the statistics when gene is removed from the module
     * removing a gene affects all leaf statistics, so bayesianScore has to be
     * recomputed from scratch afterwards.
     * 
     * @param data
     * @param Experiments
     * @param gene
     * 
     * @author tomic
     */
    public void statisticsRemoveGeneTmp(double[][] data, HashSet<Experiment> Experiments, Gene gene) {
        for (Experiment exp : Experiments) {
            TreeNode pos = this; // start at this
            while (pos.nodeStatus.equals(INTERNAL)) { // go down in tree
                if (data[pos.regulationSplit.regulator.number][exp.number] <= pos.regulationSplit.splitValue)
                    pos = pos.leftChild;
                else
                    pos = pos.rightChild;
            }
            // remove expression values for this condition and this gene from
            // present leaf statistics
            if (Double.isNaN(data[gene.number][exp.number]) == false) {
                pos.leafDistribution.statisticsTmp[0]--;
                pos.leafDistribution.statisticsTmp[1] = pos.leafDistribution.statisticsTmp[1]
                        - data[gene.number][exp.number];
                pos.leafDistribution.statisticsTmp[2] = pos.leafDistribution.statisticsTmp[2]
                        - Math.pow(data[gene.number][exp.number], 2);
            }
        }
    }
 
    /**
     * Writes the information of this node to XML.
     * 
     * @param sbn StringBuffer to write XML output.
     * 
     * @author anjos, tomic
     */
	public void toXML(StringBuffer sbn){
		
		String newline = System.getProperty("line.separator");
		
		if (this.testSplits == null)
			this.testSplits = new ArrayList<Split>();
			
		if (this.testSplitsRandom == null)
			this.testSplitsRandom = new ArrayList<Split>();
		
		if (this.nodeStatus.equals(LEAF)) {
			sbn.append("<TreeNode numChildren=\"" + 0 
					+ "\" score=\"" + this.leafDistribution.score 
					+ "\" condSet=\"" 
					+ this.exptsToString() + "\">");
			sbn.append("</TreeNode>"+newline);
		} else {
			sbn.append("<TreeNode numChildren=\"" + 2 
					+ "\" score=\"" + this.leafDistribution.score 
					+ "\" condSet=\"" + this.exptsToString() + "\">"+newline);
			if (this.testSplits.isEmpty()){
				if (this.regulationSplit == null) {
					this.regulationSplit = new Split();
					this.testSplits.add(this.regulationSplit);
				}
				sbn.append("<Regulator name=\"" + this.regulationSplit.regulator.name
						+ "\" entropy=\"" + this.regulationSplit.regulatorScore
						+ "\" splitValue=\"" + this.regulationSplit.splitValue
						+ "\">");
				sbn.append("</Regulator>"+newline);
			}
			else {
				for (Split splt : this.testSplits){
					sbn.append("<Regulator name=\"" + splt.regulator.name 
							+ "\" entropy=\"" + splt.regulatorScore
							+ "\" splitValue=\"" + splt.splitValue
							+ "\">");
					sbn.append("</Regulator>"+newline);
				}
				for (Split splt : this.testSplitsRandom){
					sbn.append("<RandomRegulator name=\"" + splt.regulator.name 
							+ "\" entropy=\"" + splt.regulatorScore
							+ "\" splitValue=\"" + splt.splitValue
							+ "\">");
					sbn.append("</RandomRegulator>"+newline);
				}
			}
			this.leftChild.toXML(sbn);
			this.rightChild.toXML(sbn);
			sbn.append("</TreeNode>");
		}
	}
	
	/**
	 * Fills the information of this node by reading an XML element.
     * 
     * @param child XML element to read.
     * 
     * @author anjos, tomic
     */
	public void fromXML(Element child){
		this.leafDistribution.score = 
			Double.parseDouble(child.getAttributeValue("score"));
		String expts = child.getAttributeValue("condSet");
		this.leafDistribution.condSet=StringtoArray(expts);
		this.nodeStatus = LEAF;
		// compute statistics and score as if node was leaf
		this.statistics();
		this.leafDistribution.bayesianScore();
		// for internal nodes, read more information
		if (Integer.parseInt(child.getAttributeValue("numChildren")) == 2){
			this.nodeStatus = INTERNAL;
			this.leftChild = new TreeNode(this);
			this.rightChild = new TreeNode(this);
			this.testSplits = new ArrayList<Split>();
			this.testSplitsRandom = new ArrayList<Split>();
			// children in XML are either Regulator, RandomRegulator or TreeNode
			// loop over "Regulators"
			for (int i=0; i<child.getChildElements("Regulator").size(); i++){
				Element e = child.getChildElements("Regulator").get(i);
				String regname = e.getAttributeValue("name");
				double entropy = Double.parseDouble(e.getAttributeValue("entropy"));
				double splitvalue = Double.parseDouble(e.getAttributeValue("splitValue"));
				// look up regulator name
				for (Gene reg : this.module.moduleNetwork.regulatorSet){
					if (reg.name.equals(regname)){
						Split splt = new Split(this, reg, splitvalue, entropy);
						this.testSplits.add(splt);
						this.regulationSplit = splt;
						break;
					}
				}
			}
			// loop over "RandomRegulators"
			for (int i=0; i<child.getChildElements("RandomRegulator").size(); i++){
				Element e = child.getChildElements("RandomRegulator").get(i);
				String regname = e.getAttributeValue("name");
				double entropy = Double.parseDouble(e.getAttributeValue("entropy"));
				double splitvalue = Double.parseDouble(e.getAttributeValue("splitValue"));
				// look up regulator name
				for (Gene reg : this.module.moduleNetwork.regulatorSet){
					if (reg.name.equals(regname)){
						Split splt = new Split(this, reg, splitvalue, entropy);
						this.testSplitsRandom.add(splt);
						break;
					}
				}
			}

			// set children
			this.leftChild.fromXML(child.getChildElements("TreeNode").get(0));
			this.rightChild.fromXML(child.getChildElements("TreeNode").get(1));
		}
	}
		 		 
	/**
	 * Converts a string of integers separated by ";" into an ArrayList
	 * 
	 *  @param expts string of integers.
	 *  @return List of integers.
	 *  
	 *  @author anjos
	 */
	public ArrayList<Integer> StringtoArray(String expts){
		ArrayList<Integer> list = new ArrayList<Integer>();
		StringTokenizer tokens = new StringTokenizer(expts,";");
		while(tokens.hasMoreTokens()){
			list.add(Integer.parseInt(tokens.nextElement().toString()));
		}
		return list;
	}
	
	/**
	 * Converts {@link LeafDistribution#condSet} into a string of integers
	 *  separated by ";".
	 *  
	 *  @author anjos
	 *  @return string of integers.
	 */
	public String exptsToString(){
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < this.leafDistribution.condSet.size(); i++) {
			result.append(this.leafDistribution.condSet.get(i));
			if (i < this.leafDistribution.condSet.size() - 1)
				result.append(";");
		}
		return result.toString();
	}   
}
