/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.modulenetwork;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Scanner;

import lemontree.utils.SetOperation;

import lemontree.ganesh.Cluster;
import lemontree.ganesh.GibbsSampler;
import lemontree.modulenetwork.Globals; //revamp
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

import cytoscape.data.annotation.OntologyTerm;



/**
 * Collection of genes with one or more possible regulation programs. 
 * The main building block in a {@link ModuleNetwork}. It extends 
 * JComponent to allow a graphical representation of the regulation
 * program(s).
 * 
 * @author tomic
 *
 */
public class Module {

	/**
	 * ModuleNetwork the module belongs to.
	 */
	public ModuleNetwork moduleNetwork;

	/**
	 * Numeric ID of module.
	 */
	public int number;

	/**
	 * Name of module (currently not used)
	 */
	public String name;

	/**
	 * List of genes in the module.
	 */ 
	public ArrayList<Gene> genes;

	/** 
	 * List of parents of the module.
	 */
	public ArrayList<Gene> parents;

	/**
	 * Map of all regulators (not just those appearing in a tree) to their weight.
	 */
	public HashMap<Gene, Double> regulatorWeights;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their weight,
	 * computed using the maximum weight in the testSplits of each node
	 */
	public HashMap<Gene, Double> regulatorWeightsMax;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their weight,
	 * computed using the sum of all weights in the testSplits of each node
	 */
	public HashMap<Gene, Double> regulatorWeightsSum;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their expected sign,
	 * computed using the sum of all weights in the testSplits of each node.
	 */
	public HashMap<Gene, Double> regulatorSign;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their left condition set,
	 * computed using the sum of all weights in the testSplits of each node. It gives for each condition
	 * the probability to belong to the condition set.
	 */
	public HashMap<Gene, ArrayList<Double>> regConditionsLeft;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their right condition set,
	 * computed using the sum of all weights in the testSplits of each node. It gives for each condition
	 * the probability to belong to the condition set.
	 */
	public HashMap<Gene, ArrayList<Double>> regConditionsRight;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their complete condition set,
	 * computed using the sum of all weights in the testSplits of each node. It gives for each condition
	 * the probability to belong to the condition set.
	 */
	public HashMap<Gene, ArrayList<Double>> regConditions;
	
	/**
	 * Map of all parents occuring in {@link Module#hierarchicalTrees} to their sign
	 */
	public HashMap<Gene, Integer> parentSigns;
	
	/**
	 * Map of all random parents occuring in {@link Module#hierarchicalTrees} to 
	 * their weight.
	 */
	public HashMap<Gene, Double> regulatorWeightsRandom;

	/**
	 * Hierarchical tree representing the regulation program. Used in 
	 * {@link ModuleNetwork#heuristicSearchMax(double, boolean, double, String[])} 
	 * when one tree per module is inferred.
	 */
	public TreeNode hierarchicalTree; 
	
	/**
	 * LEAF node to hold non-informative conditions (= conditions not in hierarchicalTree)
	 */
	public TreeNode nonInformCond;
	
	
	/**
	 * Different equally likely hierarchical trees. Based on different equally likely 
	 * clusterings of the experiments. Used in 
	 * {@link ModuleNetwork#gibbsSamplerExpts(int, int, int, int, double, boolean, boolean, double[])}.
	 */
	public ArrayList<TreeNode> hierarchicalTrees; 
	
	/**
	 * Map of each experiment to the leaf number it belongs to in each
	 * hierarchical tree
	 */
	public HashMap<Integer, ArrayList<Integer>> inExpCluster;
	
	/**
	 * Current best split for this module. Used in top-down search for regulation programs.
	 */
	public Split bestSplit;
	
	/**
	 * Module score. Computed by acting with {@link TreeNode#bayesianScore()} on
	 * {@link Module#hierarchicalTree}.
	 */
	public double moduleScore;

	/**
	 * Keep track if the list of genes has been changed compared to a previous situation.
	 * Used in {@link ModuleNetwork#heuristicSearchMax(double, boolean, double, String[])} to avoid
	 * rebuilding {@link Module#hierarchicalTree} if no genes were added/removed during
	 * {@link ModuleNetwork#moduleReassign(double)}.
	 */
	public boolean changed;
	
	/**
	 * GO
	 */
	
//	public BingoResults bingoResults;

	/**
	 *  GOsuper set to store GO "super" categories (highest term in the hierarchy) for this module.
	 */
	public LinkedHashSet<OntologyTerm> GOsuper;

	/**
	 *  GOall set to store all GO categories for this module. 
	 */
	public LinkedHashSet<OntologyTerm> GOall;

	/**
	 * ArrayList of top regulators for this module (post-processing data)
	 */
	public ArrayList<Gene> topRegulators;

	/**
	 * Ordered classes of top regulators (TF, microRNA, CNV, etc.)
	 */
	public ArrayList<ArrayList<Regulator>> topRegClasses;
	
	/**
	 * Map to top regulators indexes.
	 */
	public HashMap<String, Integer> topRegClassesMap;
	
	/**
	 * module mean gene expression value
	 */
	public double mean;
	
	/**
	 * module gene expression standard deviation value
	 */
	public double sigma;
	
	/**
	 * Empty constructor
	 */
	public Module() {
		this.changed = true;
		this.genes = new ArrayList<Gene>();
	    this.GOsuper = new LinkedHashSet<OntologyTerm>();
	    this.GOall = new LinkedHashSet<OntologyTerm>();
	    this.topRegulators = new ArrayList<Gene>();
	    this.topRegClasses = new ArrayList<ArrayList<Regulator>>();
	}

	/**
	 * Constructs initial module without regulation program.
	 * 
	 * @param M module network the module belongs to.
	 * @param number numeric ID of module.
	 * @param normalGammaPrior normal-gamma prior parameters used in the Bayesian score.
	 * 
	 */
	public Module(ModuleNetwork M, int number, double[] normalGammaPrior) {
		this.moduleNetwork = M;
		this.number = number;
		this.hierarchicalTree = new TreeNode(this, normalGammaPrior);
		for (int m =0; m<M.data[0].length; m++)
			this.hierarchicalTree.leafDistribution.condSet.add(m);
		this.changed = true;
		this.hierarchicalTrees = new ArrayList<TreeNode>();
		this.genes = new ArrayList<Gene>();
		this.parents = new ArrayList<Gene>();
		this.regulatorWeights = new HashMap<Gene,Double>();
		this.regulatorWeightsMax = new HashMap<Gene,Double>();
		this.regulatorWeightsRandom = new HashMap<Gene,Double>();
	    this.GOsuper = new LinkedHashSet<OntologyTerm>();
	    this.GOall = new LinkedHashSet<OntologyTerm>();
	    this.topRegulators = new ArrayList<Gene>();
	    this.topRegClasses = new ArrayList<ArrayList<Regulator>>();
	}

	// METHODS
	
	public void assignRegulatorsNoAcyclStoch (int level, double betaMax, int numRegAssign) {
    	System.out.println("\nAssigning regulators stochastically down to level " + level + " ...");
    	if (!this.genes.isEmpty()){
    		// check for multiple trees
    		if (this.hierarchicalTrees.isEmpty()){
    			this.hierarchicalTrees.add(this.hierarchicalTree);
    		} 
    		for (TreeNode root : this.hierarchicalTrees) {
    			// trim root below level
    			root.testLevel(level);
    			// assign regulators to remaining nodes
    			for (TreeNode node : root.getInternalNodes()){
    				node.findRegulatorBayesStoch(numRegAssign, this.moduleNetwork.regulatorSet, betaMax);
    			}
    		}
    	}
    	System.out.println("... done.");
    }
	
	/**
	 * 
	 */
	public void setInExpCluster(){
		// initialize
		this.inExpCluster = new HashMap<Integer, ArrayList<Integer>>();
		// fill
		for (int m=0; m<this.moduleNetwork.data[0].length; m++){
			this.inExpCluster.put(m, new ArrayList<Integer>());
		}
		
		for (TreeNode root : this.hierarchicalTrees){
			ArrayList<TreeNode> leaves = root.getLeafNodes();
			for (int k=0; k<leaves.size(); k++){
				for (Integer m : leaves.get(k).leafDistribution.condSet){
					this.inExpCluster.get(m).add(k);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param file
	 */
	public void printExpClusters(String file){
		this.setInExpCluster();
		try{
    		File f = new File(file);
    		PrintWriter pw = new PrintWriter(f);
    		for (int e=0; e<this.moduleNetwork.data[0].length; e++){
    			String s = String.format("%s\t%d", "exp", e);
    			for (int m : this.inExpCluster.get(e))
    				s += "\t" + m;
    			pw.println(s);
    		}
    		pw.close();
    	} catch (IOException e){
			System.out.println(e);
		}
	}
	
	public ArrayList<ArrayList<Integer>> readExpClusters (String file, int numCluster){
	       System.out.println("Reading experiment clusters ...");
	        int assigned = 0;
	        
	        ArrayList<ArrayList<Integer>> expclust = new ArrayList<ArrayList<Integer>>(numCluster);
	        for (int k=0; k<numCluster; k++)
	        	expclust.add(new ArrayList<Integer>());
	        
	        try {
	            Scanner clusterScanner = new Scanner(new File(file));
	            clusterScanner.useDelimiter("\\n");
	            // walk through clusterfile, skip comment lines, check if counter starts with 0 or 1
	            int first = 1;
	            while (clusterScanner.hasNext()) {
	                Scanner line = new Scanner(clusterScanner.next().toString());
	                line.useDelimiter("\\s");
	                // read experiment number
	                int exp = line.nextInt();
	                int number = line.nextInt();
	                if (number == 0){
	                	first = 0;
	                }
	                expclust.get(number-first).add(exp);
	                assigned += 1;
	                line.close();
	            }
	            clusterScanner.close();
	            System.out.println("... assigned " + assigned	+ " experiments to "
	                    + expclust.size() + " clusters.");
	            System.out.println();
	            clusterScanner.close();
	        } catch (FileNotFoundException e) {
	            System.out.println("Cluster file not found.");
	        }
	        
	        return expclust;
	}
	
	
	/**
	 * 
	 * @param numRuns
	 * @param numSteps
	 * @param burnIn
	 * @param sampleStep
	 * @param scoregain
	 * @param useBHCscore
	 */
	public void gibbsSampler(int numRuns, int numSteps, int burnIn, int sampleStep,
    		double scoregain, boolean useBHCscore){
		System.out.println("Gibbs sampling experiment partitions ...");
		this.hierarchicalTrees = new ArrayList<TreeNode>();
		HashSet<ArrayList<ArrayList<Integer>>> clusterList 
			= this.gibbsSamplerClustering(numRuns, numSteps, burnIn, sampleStep);
		for (ArrayList<ArrayList<Integer>> cluster : clusterList){
			this.hierarchicalTree = this.hierarchicalClusteringOrdered(this.partition2TreeNode(cluster), useBHCscore);
			this.hierarchicalTree.testScore(scoregain);
			this.hierarchicalTrees.add(this.hierarchicalTree);
		}
		System.out.println("... done");
	}
	
	/**
	 * Gibbs sampling to construct several equally likely
	 * partitions of the experiments for this module.
	 * 
	 * @see CRC
     * @param numRuns number of Gibbs sampler runs.
     * @param numSteps number of steps in each run of the Gibbs sampler.
     * @param burnIn number of intial steps before sampling starts.
     * @param sampleStep number of steps between successive samples.
     * 
     * @return Set of partitions. Each partition is a list of lists of integers
     * (referring to the experiment ID's) forming the sets of the partitions.
     * 
     * @author anjos, tomic
	 */
	public HashSet<ArrayList<ArrayList<Integer>>> gibbsSamplerClustering
		(int numRuns, int numSteps, int burnIn, int sampleStep) {

		// set to store different experiment clusterings
		HashSet<ArrayList<ArrayList<Integer>>> expClustering = new HashSet<ArrayList<ArrayList<Integer>>>();

		// get the data submatrix for this module
		DoubleMatrix2D data = new DenseDoubleMatrix2D(this.moduleNetwork.data);
		int[] rowList = new int[genes.size()];
		int i=0;
		for (Gene gene : genes)	{
			rowList[i]=gene.number;
			i++;
		}
		double[][] moduledata = Algebra.DEFAULT.transpose(data.viewSelection(rowList,null)).toArray();

		// cluster experiments using GaneSh
		for(int j=1; j<=numRuns; j++){
			GibbsSampler C = new GibbsSampler(moduledata);
			C.AssignRandomrowSinglecolumn();
			for(int k=1; k<=numSteps; k++){
				C.Cluster();
				if(k > burnIn && (numSteps - k) % sampleStep == 0){
					// list of experiment clusters 
					ArrayList<ArrayList<Integer>> expClusters = new ArrayList<ArrayList<Integer>>();
					for(Cluster clust : C.ClusterSet){
						ArrayList<Integer> expList 
							= new ArrayList<Integer>(clust.RowSet.size());
						for(Integer in : clust.RowSet){
							expList.add(in);
						}
						expClusters.add(expList);
					}
					expClustering.add(expClusters);
				}
			}
		}
		return(expClustering);
	} 

	/**
	 * Converts a partition of the experiments to a list of {@link TreeNode}s. 
	 * For each set of experiments, creates a TreeNode with this set as its condition set. 
	 * This list of TreeNodes can be used as input for 
	 * {@link Module#hierarchicalClustering(List, boolean)}.
	 * 
	 * @param clusters partition of experiments, list of lists of integers.
	 * @return List of {@link TreeNode}s.
	 * 
	 * @author tomic
	 *   
	 */
	public List<TreeNode> partition2TreeNode (ArrayList<ArrayList<Integer>> clusters) {
		List<TreeNode> treeList = new ArrayList<TreeNode>();
		
		for (ArrayList<Integer> list : clusters){
			TreeNode node = new TreeNode(this, this.hierarchicalTree.leafDistribution.normalGammaPrior);
			node.leafDistribution.condSet = list;
			node.statistics();
			node.leafDistribution.bayesianScore();
			treeList.add(node);
		}
		
		return treeList;
	}
	
	
	/**
	 * Takes the output of {@link Module#gibbsSamplerClustering(int, int, int, int)} and creates a list
	 * of condition sets which are clustered together in each Gibbs sample and have a minimal size.
	 * 
	 * @param partitionList
	 * @param minSize
	 * @return
	 * 
	 * @author tomic
	 */
	public ArrayList<ArrayList<Integer>> makeConsensusCondSets (HashSet<ArrayList<ArrayList<Integer>>> partitionList, int minSize){
		// create list of intersections of all partitions
		Iterator<ArrayList<ArrayList<Integer>>> it = partitionList.iterator();
		ArrayList<ArrayList<Integer>> list = it.next();
		while (it.hasNext()){
			list = SetOperation.prodPart(list, it.next());
		}
		// find intersections that are too small
		ArrayList<ArrayList<Integer>> remove = new ArrayList<ArrayList<Integer>>();
		for (ArrayList<Integer> l : list)
			if (l.size() < minSize)
				remove.add(l);
		for (ArrayList<Integer> l : remove)
			list.remove(l);
		// done
		return list;
	}
	
	/**
	 * Takes the leaf nodes of all trees in {@link Module#hierarchicalTrees} and creates a list
	 * of condition sets which are clustered together in each tree and have a minimal size.
	 * 
	 * @param partitionList
	 * @param minSize
	 * @return
	 * 
	 * @author tomic
	 */
	public ArrayList<ArrayList<Integer>> makeConsensusCondSetsFromTrees (int minSize){
		// create partition list
		HashSet<ArrayList<ArrayList<Integer>>> partitionList = new HashSet<ArrayList<ArrayList<Integer>>> ();
		for (TreeNode root : this.hierarchicalTrees){
			partitionList.add(root.condSets());
		}
		// create list of intersections of all partitions
		Iterator<ArrayList<ArrayList<Integer>>> it = partitionList.iterator();
		ArrayList<ArrayList<Integer>> list = it.next();
		while (it.hasNext()){
			list = SetOperation.prodPart(list, it.next());
		}
		// find intersections that are too small
		ArrayList<ArrayList<Integer>> remove = new ArrayList<ArrayList<Integer>>();
		for (ArrayList<Integer> l : list)
			if (l.size() < minSize)
				remove.add(l);
		for (ArrayList<Integer> l : remove)
			list.remove(l);
		// done
		return list;
	}
	
	/**
	 * Makes a consensus tree from {@link Module#hierarchicalTrees}, stores the result in {@link Module#hierarchicalTree},
	 * and empties {@link Module#hierarchicalTrees}.
	 * @param minSize
	 * 
	 * @author tomic
	 */
	public void makeConsensusTree (int minSize){
		ArrayList<ArrayList<Integer>> condsets =  this.makeConsensusCondSetsFromTrees(minSize);
		if (condsets.size() > 0){
			// cluster hierarchically
			this.hierarchicalTree = this.hierarchicalClusteringOrdered(this.partition2TreeNode(condsets), true);
		} else {
			this.hierarchicalTree.nodeStatus = "leaf";
			this.hierarchicalTree.leftChild = null;
			this.hierarchicalTree.rightChild = null;
		}
		// empty trees
		this.hierarchicalTrees = new ArrayList<TreeNode>();
	}
	
	/**
	 * Hierarchical clustering of experiments (or experiment clusters).
	 * This builds a full hierarchical clustering tree of the data in the
	 * module, which can then be converted in a regulation program tree.
	 * 
	 * Usage of {@link TreeNode} and {@link LeafDistribution} is somewhat 
	 * non-standard:
	 * <ul> 
	 * <li>{@link TreeNode#testGene} and {@link TreeNode#splitValue} are not set on 
	 * internal nodes.</li>
	 * <li>Internal nodes have a {@link LeafDistribution} as well, containing 
	 * the statistics of all leaves (=experiments) below it.</li> 
 	 * </ul>
 	 * 
 	 * The hierarchical clustering proceeds by merging in each step the pair of 
 	 * nodes with the smallest {@link TreeNode#mergeScore}. This can be either
 	 * the score difference between the merged nodes and the sum of scores of each 
 	 * node, or the score based on Bayesian Hierarchical Clustering which also 
 	 * takes into account the other partitions below the pair of nodes.
	 * 
	 * @param treeList list of TreeNodes to be clustered. This can be each 
	 * experiment in an individual leaf node such as returned by {@link Module#initTreeList()} or
	 * {@link Module#initTreeList(ArrayList)}, 
	 * or a set of experiment clusters such as returned by 
	 * {@link Module#partition2TreeNode(ArrayList)}.
	 * @param useBHCscore use merge score from Bayesian Hierarchical Clustering or 
     * simple score difference when hierarchically clustering experiments.
	 * @return root node of the hierarchical tree.
	 * 
	 * @author tomic
	 * 
	 */
	public TreeNode hierarchicalClustering(List<TreeNode> treeList, boolean useBHCscore) {
		int numNode = treeList.size();
		double[] normalGammaPrior = treeList.get(0).leafDistribution.normalGammaPrior;
		List<TreeNode> pairTreeList = new ArrayList<TreeNode>(numNode*(numNode-1)/2);
		
		// initialize list of pairwise merges
		for (int i = 0; i < numNode; i++){
			for (int j = 0; j < i; j++) {
				TreeNode tree = new TreeNode(this, normalGammaPrior);
				tree.nodeStatus = "internal";
				// child nodes, order by means of normal distribution at nodes
				if (treeList.get(i).leafDistribution.gaussianParam[0] < 
						treeList.get(j).leafDistribution.gaussianParam[0]){
					tree.leftChild = treeList.get(i);
					tree.rightChild = treeList.get(j);
				} else {
					tree.leftChild = treeList.get(j);
					tree.rightChild = treeList.get(i);
				}
				// statistics of merged node
				for (int s = 0; s < 3; s++)
					tree.leafDistribution.statistics[s] = tree.leftChild.leafDistribution.statistics[s]
							+ tree.rightChild.leafDistribution.statistics[s];
				// union of condition sets
				tree.leafDistribution.condSet = SetOperation.union(tree.leftChild.leafDistribution.condSet, 
						tree.rightChild.leafDistribution.condSet);
				// compute score
				tree.leafDistribution.bayesianScore();
				tree.setMergeScore(useBHCscore);
				// add to list
				pairTreeList.add(tree);
			}
		}

		// merge trees in treeList until one tree is left
		TreeNode bestTree = new TreeNode(this, normalGammaPrior);
		bestTree.mergeScore = 1.0;
		while (treeList.size() > 1) {
			// find the pair tree with the highest posterior probability
			//Collections.sort(pairTreeList);
			bestTree = Collections.max(pairTreeList);
			// set bestTree as parent of its children
			bestTree.leftChild.parentNode = bestTree;
			bestTree.rightChild.parentNode = bestTree;
			// remove children of bestTree from treeList; bestTree will be added at the end
			treeList.remove(bestTree.leftChild);
			treeList.remove(bestTree.rightChild);
			// loop over the pair list and remove trees whose condSet intersects
			// with the condSet of the highest scoring tree, including the highest one
			// itself
			Iterator<TreeNode> it = pairTreeList.iterator();
			while (it.hasNext()) {
				TreeNode tree = it.next();
				ArrayList<Integer> intersect 
					= SetOperation.intersection(tree.leafDistribution.condSet, 
							bestTree.leafDistribution.condSet);
				if (!intersect.isEmpty())
					it.remove();
			}
			// add new trees which merge the bestTree with the other trees in
			// treeList;
			Iterator<TreeNode> it2 = treeList.iterator();
			while (it2.hasNext()) {
				TreeNode ittree = it2.next();
				TreeNode tree = new TreeNode(this, normalGammaPrior);
				tree.nodeStatus = "internal";
				// child nodes, order means of normal distributions at the nodes
				if (ittree.leafDistribution.gaussianParam[0] < 
						bestTree.leafDistribution.gaussianParam[0]){
					tree.leftChild = ittree;
					tree.rightChild = bestTree;
				} else {
					tree.leftChild = bestTree;
					tree.rightChild = ittree;
				}
				// statistics of merged node
				for (int s = 0; s < 3; s++)
					tree.leafDistribution.statistics[s] = tree.leftChild.leafDistribution.statistics[s]
							+ tree.rightChild.leafDistribution.statistics[s];
				// union of condition sets
				tree.leafDistribution.condSet = SetOperation.union(tree.leftChild.leafDistribution.condSet, 
						tree.rightChild.leafDistribution.condSet);
				// compute scores
				tree.leafDistribution.bayesianScore();
				tree.setMergeScore(useBHCscore);
				// add to list
				pairTreeList.add(tree);
			}
			// add bestTree again to treeList
			treeList.add(bestTree);
		}
		// 1 tree left, this is the one
		return treeList.get(0);
	}
	
	/**
	 * Hierarchical clustering of a set of leaf nodes. Similar to {@link Module#hierarchicalClustering(List, boolean)},
	 * but this method orders the nodes by their mean expression level and only allows merges
	 * between neighboring nodes. This is faster and prevents the 'wrong' linking that is sometimes seen in
	 * {@link Module#hierarchicalClustering(List, boolean)}.
	 * 
	 * @param treeList
	 * @param useBHCscore
	 * @return
	 * 
	 * @author tomic
	 */
	public TreeNode hierarchicalClusteringOrdered(List<TreeNode> treeList, boolean useBHCscore) {
		double[] normalGammaPrior = treeList.get(0).leafDistribution.normalGammaPrior;
		// list of possible merges
		List<TreeNode> pairList = this.pairTreeList(treeList, useBHCscore);
		// merge trees in treeList until one tree is left
		TreeNode bestTree = new TreeNode(this, normalGammaPrior);
		bestTree.mergeScore = 1.0;
		while (treeList.size() > 1) {
			// find the pair tree with the highest posterior probability
			bestTree = Collections.max(pairList);
			// set bestTree as parent of its children
			bestTree.leftChild.parentNode = bestTree;
			bestTree.rightChild.parentNode = bestTree;
			// remove children of bestTree from treeList and add bestTree
			treeList.remove(bestTree.leftChild);
			treeList.remove(bestTree.rightChild);
			treeList.add(bestTree);
			// sort and pair again
			pairList = this.pairTreeList(treeList, useBHCscore);
		}
		// 1 tree left, this is the one
		return treeList.get(0);
	}
	
	/**
	 * Creates from a list of tree nodes a new list of tree nodes by merging
	 * neighboring pairs (neighbor in the sense of data mean). 
	 * 
	 * @param treeList
	 * @param useBHCscore
	 * @return
	 * 
	 * @author tomic
	 */
	public List<TreeNode> pairTreeList(List<TreeNode> treeList, boolean useBHCscore){
		double[] normalGammaPrior = treeList.get(0).leafDistribution.normalGammaPrior;
		// sort treeList by increasing mean
		Comparator mc = new MeanComparator();
		Collections.sort(treeList, mc);
		// list of merges between neighboring pairs
		List<TreeNode> pairList = new ArrayList<TreeNode>(treeList.size()-1);
		for (int i = 0; i < treeList.size()-1; i++){
			TreeNode tree = new TreeNode(this, normalGammaPrior);
			// attach children
			tree.nodeStatus = "internal";
			// child nodes	
			tree.leftChild = treeList.get(i);
			tree.rightChild = treeList.get(i+1);
			// statistics of merged node
			for (int s = 0; s < 3; s++)
				tree.leafDistribution.statistics[s] = tree.leftChild.leafDistribution.statistics[s]
						+ tree.rightChild.leafDistribution.statistics[s];
			tree.leafDistribution.statistics[3] = tree.leftChild.leafDistribution.statistics[3]; //revamp
			// union of condition sets
			tree.leafDistribution.condSet = SetOperation.union(tree.leftChild.leafDistribution.condSet, 
					tree.rightChild.leafDistribution.condSet);
			// compute score
			tree.leafDistribution.bayesianScore();
			tree.setMergeScore(useBHCscore);
			// dummy regulator
			tree.regulationSplit = new Split();
			tree.testSplits = new ArrayList<Split>();
			// add to list
			pairList.add(tree);
		}
		return pairList;
	}
	
	/**
	 * From a list of experiment clusters (in the form of leaf nodes), find a
	 * subset which form mutually significant splits. The siginificance test
	 * is whatever is implemented in 
	 * {@link Module#compareSignificance(TreeNode, TreeNode)} 
	 * The idea is to apply this method
	 * to the list of leaves obtained from {@link Module#makeConsensusTree(int)}.
	 * 
	 * @param treeList
	 * @param fact multiplication factor in significance test
	 * @return
	 * 
	 * @author tomic
	 */
	public TreeNode findSignificantSplits (List<TreeNode> treeList, double fact){
		int numTree = treeList.size();
		// help tree list to store selected nodes from treeList
		List<TreeNode> tmpList = new ArrayList<TreeNode>(); 
		// sort treeList by increasing mean
		Comparator mc = new MeanComparator();
		Collections.sort(treeList, mc);
		
		if (treeList.size()==1)
			return treeList.get(0);
		else {
			// compare first and last
			boolean init = this.compareSignificance(treeList.get(0), treeList.get(numTree-1), fact);
			if (!init){ // put everything in 1 leaf
				TreeNode tree = this.hierarchicalClusteringOrdered(treeList, true);
				tree.testLevel(-1);
				return tree;
			} else {// we have something
				tmpList.add(treeList.get(0));
				tmpList.add(treeList.get(numTree-1));
				treeList.remove(numTree-1);
				treeList.remove(0);
				if (!treeList.isEmpty()){// recursion
					TreeNode tree2 = this.findSignificantSplits(treeList, fact);
					// check if we can add tree2
					boolean add = true;
					for (TreeNode node : tmpList){
						add = add && this.compareSignificance(node, tree2, fact);
					}
					if (add)
						tmpList.add(tree2);
				}
				// link chosen nodes into 1 tree
				TreeNode tree = this.hierarchicalClusteringOrdered(tmpList, true);
				return tree;
			}
		}
	}
	
	public void expandSignificantSplits(double fact){
		this.setNonInformCond();
		if (!this.nonInformCond.leafDistribution.condSet.isEmpty()){
			//int num = 0;
			ArrayList<TreeNode> leaves = this.hierarchicalTree.getLeafNodes();
			for (int m : this.nonInformCond.leafDistribution.condSet){
				double mean = 0;
				for (Gene gene : this.genes){
					mean += this.moduleNetwork.data[gene.number][m];
				}
				mean /= this.genes.size();
				// add m to a leaf if it is within (fact x standard deviation) of the mean;
				// because of significant splits this can happen with at most 1 leaf, if "fact"
				// has the same value as used to determine significance
				for (TreeNode leaf : leaves){
					if (Math.abs(leaf.leafDistribution.gaussianParam[0] - mean) < fact*leaf.leafDistribution.gaussianParam[1]){
						// add m
						leaf.leafDistribution.condSet.add(m);
						// update statistics
						leaf.statistics();
						// update score and normal dist param
						leaf.leafDistribution.bayesianScore();
						// done
						//num += 1;
						break;
					}
				}
			}
			this.setNonInformCond();
		}
	}
	
	/**
	 * Checks if the normal distributions for two nodes are significantly separated, i.e.,
	 * if (mu1 - mu2) > fact*(sigma1 + sigma2)
	 * 
	 * @param node1
	 * @param node2
	 * @param fact
	 * @return
	 * 
	 * @author tomic
	 */
	public boolean compareSignificance(TreeNode node1, TreeNode node2, double fact){
		Double mu1 = node1.leafDistribution.gaussianParam[0];
		Double sigma1 = node1.leafDistribution.gaussianParam[1];
		Double mu2 = node2.leafDistribution.gaussianParam[0];
		Double sigma2 = node2.leafDistribution.gaussianParam[1];
		if (Math.abs(mu1-mu2) > fact*(sigma1+sigma2))
			return true;
		else
			return false;
	}
	
	
	
	/**
	 * Putd each experiment in its own leaf node. This list of TreeNodes can be
	 * used as input for {@link Module#hierarchicalClustering(List, boolean)}.
	 * 
	 * @author tomic
	 */
	public List<TreeNode> initTreeList () {
		double[][] data = this.moduleNetwork.data;
		List<TreeNode> treeList = new ArrayList<TreeNode>(data[0].length);
		for (int i=0; i<data[0].length; i++) {
			TreeNode tree = new TreeNode(this, this.moduleNetwork.normalGammaPrior);
			tree.nodeStatus = "leaf";
			tree.leafDistribution.condSet.add(i);
			for (Gene gene : this.genes)
				if (!Double.isNaN(data[gene.number][i])) {
					tree.leafDistribution.statistics[0]++;
					tree.leafDistribution.statistics[1] += data[gene.number][i];
					tree.leafDistribution.statistics[2] += Math.pow(data[gene.number][i], 2);
				}
			tree.leafDistribution.bayesianScore();
			tree.mergeScore = 1.0;
			treeList.add(tree);
		}
		return treeList;
	}
	
	/**
	 * Puts each experiment of a selected list in its own leaf node. This list of TreeNodes 
	 * can be used as input for {@link Module#hierarchicalClustering(List, boolean)}.
	 * 
	 * @param condSet for each experiment in this list a leaf node is created.
	 * 
	 * @author tomic
	 * 
	 */
	public List<TreeNode> initTreeList (ArrayList<Integer> condSet) {
		List<TreeNode> treeList = new ArrayList<TreeNode>(condSet.size());
		double[][] data = this.moduleNetwork.data;
		for (int i : condSet) {
			TreeNode tree = new TreeNode(this, this.moduleNetwork.normalGammaPrior);
			tree.nodeStatus = "leaf";
			tree.leafDistribution.condSet.add(i);
			for (Gene gene : this.genes)
				if (!Double.isNaN(data[gene.number][i])) {
					tree.leafDistribution.statistics[0]++;
					tree.leafDistribution.statistics[1] += data[gene.number][i];
					tree.leafDistribution.statistics[2] += Math.pow(data[gene.number][i], 2);
				}
			// revamp: only works if there's one node per condition
			tree.leafDistribution.statistics[3] = this.moduleNetwork.condition_weight[i]; //revamp
			tree.leafDistribution.bayesianScore();
			tree.mergeScore = 1.0;
			treeList.add(tree);
		}
		return treeList;
	}

	
	public void setRegulatorSign (){
		// initialize
		//int numCond = this.moduleNetwork.data[0].length;
		this.regulatorSign = new HashMap<Gene,Double>();
		for (Gene reg : this.regulatorWeights.keySet())
			this.regulatorSign.put(reg, 0.0);
		// compute
		if (this.hierarchicalTrees.isEmpty())
			this.hierarchicalTrees.add(this.hierarchicalTree);
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = Math.pow((double)node.leafDistribution.condSet.size()/(double)this.moduleNetwork.data[0].length,1);
				for (Split splt : node.testSplits){
					double pIn = this.regulatorSign.get(splt.regulator);
					double pAdd = splt.sign * Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact
						/this.regulatorWeights.get(splt.regulator);
					this.regulatorSign.put(splt.regulator, pIn + pAdd);
				}
			}
		}
	}
	
	
	public void setRegConditions (){
		// initialize
		int numCond = this.moduleNetwork.data[0].length;
		this.regConditionsLeft = new HashMap<Gene, ArrayList<Double>>();
		this.regConditionsRight = new HashMap<Gene, ArrayList<Double>>();
		this.regConditions = new HashMap<Gene, ArrayList<Double>>();
		for (Gene reg : this.regulatorWeights.keySet()){
			this.regConditionsLeft.put(reg, new ArrayList<Double>(numCond));
			this.regConditionsRight.put(reg, new ArrayList<Double>(numCond));
			this.regConditions.put(reg, new ArrayList<Double>(numCond));
			for (int m=0; m<numCond; m++){
				this.regConditionsLeft.get(reg).add(0.0);
				this.regConditionsRight.get(reg).add(0.0);
				this.regConditions.get(reg).add(0.0);
			}
		}
		// compute
		if (this.hierarchicalTrees.isEmpty())
			this.hierarchicalTrees.add(this.hierarchicalTree);
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = Math.pow((double)node.leafDistribution.condSet.size()/(double)this.moduleNetwork.data[0].length,1);
				for (Split splt : node.testSplits){
					for (int m : node.leftChild.leafDistribution.condSet){
						double pIn = this.regConditionsLeft.get(splt.regulator).get(m);
						double pAdd = Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact
							/this.regulatorWeights.get(splt.regulator);
						this.regConditionsLeft.get(splt.regulator).set(m, pIn + pAdd);
						this.regConditions.get(splt.regulator).set(m, pIn + pAdd);
					}
					for (int m : node.rightChild.leafDistribution.condSet){
						double pIn = this.regConditionsRight.get(splt.regulator).get(m);
						double pAdd = Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact
							/this.regulatorWeights.get(splt.regulator);
						this.regConditionsRight.get(splt.regulator).set(m, pIn + pAdd);
						this.regConditions.get(splt.regulator).set(m, pIn + pAdd);
					}
				}
			}
		}
	}
	/**
	 * Computes the regulator weights. Sums regulator weights at
	 * all nodes.
	 *
	 * @author tomic
	 */
	public void setRegulatorWeights(){
		this.regulatorWeights = new HashMap<Gene, Double>();
		if (this.hierarchicalTrees.isEmpty())
			this.hierarchicalTrees.add(this.hierarchicalTree);
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = Math.pow((double)node.leafDistribution.condSet.size()/(double)this.moduleNetwork.data[0].length,1);
				for (Split splt : node.testSplits){
					if (!this.regulatorWeights.containsKey(splt.regulator))
						this.regulatorWeights.put(splt.regulator, 0.0);
					double w = this.regulatorWeights.get(splt.regulator)
						+ Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact;
					this.regulatorWeights.put(splt.regulator, w);
				}
			}
		}
	}
	

	/**
	 * Computes the weight of each regulator at the nodes of all
	 * hierarchical trees for this module. The weight is given by
	 * the average probability of a correct prediction, i.e.,
	 * w = exp(regulatorScore/node.numConditions), summed over all
	 * occurrences of a regulator in the set of regulation trees.
	 *
	 * @param downWeightLevel reduce weights of assignments below this tree level
	 * @param downWeightFact multiplication factor by which to reduce weights
	 * @deprecated use setRegulatorWeights instead
	 * @author tomic
	 */
	public void setRegulatorWeightsSum(int downWeightLevel, double downWeightFact) {
		this.regulatorWeightsSum = new HashMap<Gene, Double>();
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = 1.0;
				for (Split splt : node.testSplits){
					if (!this.regulatorWeightsSum.containsKey(splt.regulator))
						this.regulatorWeightsSum.put(splt.regulator, 0.0);
					double w = this.regulatorWeightsSum.get(splt.regulator)
						+ Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact;
					this.regulatorWeightsSum.put(splt.regulator, w);
				}
			}
		}
	}
	
	
	/**
	 * Computes the weight of each regulator at the nodes of all
	 * hierarchical trees for this module. The weight is given by
	 * the average probability of a correct prediction, i.e.,
	 * w = exp(regulatorScore/node.numConditions), maximized over all
	 * occurrences of a regulator in the set of regulation trees.
	 *
	 * @param downWeightLevel reduce weights of assignments below this tree level
	 * @param downWeightFact multiplication factor by which to reduce weights
	 * @deprecated use setRegulatorWeights instead
	 * @author tomic
	 */
	public void setRegulatorWeightsMax(int downWeightLevel, double downWeightFact){
		this.regulatorWeightsMax = new HashMap<Gene, Double>();
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = (double)node.leafDistribution.condSet.size()/(double)this.moduleNetwork.data[0].length;
				for (Split splt : node.testSplits){
					double score = Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact;
					if (!this.regulatorWeightsMax.containsKey(splt.regulator))
						this.regulatorWeightsMax.put(splt.regulator, score);
					else{
						double w = this.regulatorWeightsMax.get(splt.regulator);
						if (score > w)
							this.regulatorWeightsMax.put(splt.regulator, score);
					}
				}
			}
		}
	}
	
	

	
	/**
	 * Computes the weight of each random regulator at the nodes of all
	 * hierarchical trees for this module. The weight is given by
	 * the average probability of a correct prediction, i.e.,
	 * w = exp(regulatorScore/node.numConditions), summed over all
	 * occurrences of a random regulator in the set of regulation trees.
	 *
	 * @author tomic
	 */
	public void setRegulatorWeightsRandom() {
		this.regulatorWeightsRandom = new HashMap<Gene, Double>();
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = (double)node.leafDistribution.condSet.size()/(double)this.moduleNetwork.data[0].length;
				for (Split splt : node.testSplitsRandom){
					if (!this.regulatorWeightsRandom.containsKey(splt.regulator))
						this.regulatorWeightsRandom.put(splt.regulator, 0.0);
					double w = this.regulatorWeightsRandom.get(splt.regulator)
						+ Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact;
					this.regulatorWeightsRandom.put(splt.regulator, w);
				}
			}
		}
	}
	
	/**
	 * Computes the weight of each random regulator at the nodes of all
	 * hierarchical trees for this module. The weight is given by
	 * the average probability of a correct prediction, i.e.,
	 * w = exp(regulatorScore/node.numConditions), summed over all
	 * occurrences of a random regulator in the set of regulation trees.
	 *
	 * @param downWeightLevel reduce weights of assignments below this tree level
	 * @param downWeightFact multiplication factor by which to reduce weights
	 * @deprecated use setRegulatorWeightsRandom instead
	 * @author tomic
	 */
	public void setParentWeightsRandomSum(int downWeightLevel, double downWeightFact) {
		this.regulatorWeightsRandom = new HashMap<Gene, Double>();
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				double fact = 1.0;
				if (node.treeLevel() > downWeightLevel)
					fact = downWeightFact;
				for (Split splt : node.testSplitsRandom){
					if (!this.regulatorWeightsRandom.containsKey(splt.regulator))
						this.regulatorWeightsRandom.put(splt.regulator, 0.0);
					double w = this.regulatorWeightsRandom.get(splt.regulator)
						+ Math.exp(splt.regulatorScore/node.leafDistribution.condSet.size())*fact;
					this.regulatorWeightsRandom.put(splt.regulator, w);
				}
			}
		}
	}
	
	/**
	 * Sets the sign (activator/repressor) for each regulator at the nodes of all
	 * hierarchical trees for this module.
	 *
	 * @return false if a regulator occurs with different sign at different nodes
	 * 
	 * @author tomic
	 */
	public boolean setParentSigns(){
		boolean consistent = true;
		this.parentSigns = new HashMap<Gene, Integer>();
		for (TreeNode regTree: this.hierarchicalTrees) {
			for (TreeNode node : regTree.getInternalNodes()){
				for (Split splt : node.testSplits){
					splt.setSumData();
					splt.computeSign();
					if (!this.parentSigns.containsKey(splt.regulator))
						this.parentSigns.put(splt.regulator, splt.sign);
					else // check consistency
						if (splt.sign != this.parentSigns.get(splt.regulator))
							consistent = false;
				}
			}
		}
		return consistent;
	}

	/** 
	 * Sets the leaf node holding the non-informative conditions
	 * 
	 * @author tomic
	 */
	public void setNonInformCond(){
		this.nonInformCond = new TreeNode(this, this.moduleNetwork.normalGammaPrior);
		for (int m=0; m<this.moduleNetwork.data[0].length; m++)
			if (!this.hierarchicalTree.leafDistribution.condSet.contains(m))
				this.nonInformCond.leafDistribution.condSet.add(m);
		this.nonInformCond.statistics();
	}
	
	/**
	 * Validates the regulation program in hierarchicalTree on the conditions
	 * in nonInformCond
	 * 
	 * @author tomic
	 */
	public double validateRegProgram(){
		this.setNonInformCond();
		double p = 0.0;
		EnsembleCalculator calc = new EnsembleCalculator(this.moduleNetwork);
		for (int m : this.nonInformCond.leafDistribution.condSet){
			double[] xExp = calc.getDataColumn(m);
			for (Gene gene : this.genes){
				double pi = 0.0;
				for (TreeNode leaf : this.hierarchicalTree.getLeafNodes()){
					double w = leaf.pathWeight(xExp);
					pi += w * leaf.leafDistribution.evaluate(xExp[gene.number]);
				}
				p += Math.log(pi);
			}
		}
		return Math.exp(p/((double)this.nonInformCond.leafDistribution.condSet.size()*this.genes.size()));
	}
	
	

	/**
	 * Learn a regulation tree top-down using exact splits on a given set of
	 * parent for this module.
	 * 
	 * @author tomic
	 *
	 */
	public void learnTreeTopDown(ArrayList<Gene> regs, int maxParents){
		System.out.println("Learning tree for module " + this.number);
		// reset tree and recompute statistics
		this.parents = new ArrayList<Gene>();
		this.hierarchicalTrees = new ArrayList<TreeNode>();
		this.hierarchicalTree = new TreeNode(this, this.moduleNetwork.normalGammaPrior);
		for (int m = 0; m < this.moduleNetwork.data[0].length; m++)
			this.hierarchicalTree.leafDistribution.condSet.add(m);
		this.hierarchicalTree.statistics();
		// update score value
		this.moduleScore = this.hierarchicalTree.bayesianScore();
		// keep splitting as long as we have parents and can improve the score
		//ArrayList<Gene> regs = this.parents;
		boolean foundSplit = true;
		while (foundSplit && this.parents.size() < maxParents){
			// find current best split
			foundSplit = this.bestSplitTopDown(regs);
			// apply best split
			if (foundSplit){
				System.out.println("\t found split: " + this.bestSplit.regulator.name);
				this.bestSplit.node.split(this.bestSplit.regulator, this.bestSplit.splitValue);
				this.parents.add(this.bestSplit.regulator);
			}
		}
	}
	
	/**
	 * Finds the best split in a module using exact splits on a subset of regulators.
	 * If one is found it is stored in {@link Module#bestSplit}.
	 * 
	 * @return <code>true</code> if a split is found.
	 * 
	 * @author tomic
	 */
	public boolean bestSplitTopDown(ArrayList<Gene> regulators) {
		HashSet<TreeNode> leaves = new HashSet<TreeNode>();
		TreeNode optLeaf = new TreeNode();
		double optSplit = 0.0;
		Gene optReg = new Gene();
		double optModuleScore;
		double maxModuleScore;
		double scoreGain;

		maxModuleScore = this.moduleScore;

		this.hierarchicalTree.gatherLeaves(leaves);
		for (TreeNode leaf : leaves)
			for (Gene reg : regulators) {
				for (Double splt : reg.possibleSplitValues) {
					// try to split
					optModuleScore = leaf.trySplit(reg, splt);
					if (optModuleScore - maxModuleScore > 0.0) {
						optLeaf = leaf;
						optReg = reg;
						optSplit = splt;
						maxModuleScore = optModuleScore;
					}
				}
			}

		// we have found a best split if max has changed its value in the
		// previous loops
		if (maxModuleScore - this.moduleScore > 0.0) {
			scoreGain = maxModuleScore - this.moduleScore;
			this.bestSplit = new Split(optLeaf, optReg, optSplit, 0.0, scoreGain);
			return true;
		} else {
			this.bestSplit = new Split();
			return false;
		}
	}
	
	/**
	 * Finds the best split in a module which preseves acyclicity. 
	 * If one is found it is stored in {@link Module#bestSplit}.
	 * 
	 * @param maxParents maximum number of parents for a module.
	 * @return <code>true</code> if a split is found.
	 * 
	 * @author tomic
	 */
	public boolean bestSplitTopDown(int maxParents) {
		HashSet<TreeNode> leaves = new HashSet<TreeNode>();
		TreeNode optLeaf = new TreeNode();
		double optSplit = 0.0;
		Gene optReg = new Gene();
		double optModuleScore;
		double maxModuleScore;
		double scoreGain;

		System.out.println("Looking for best split in module " + this.number);
		maxModuleScore = this.moduleScore;

		// only proceed if module has less than maxParents parents
		if (this.parents.size() < maxParents) {
			for (Gene reg : this.moduleNetwork.regulatorSet) {
				// only proceed if adding reg to parents preserves acyclicity
				if (!this.moduleNetwork.moduleGraph.testCycle(reg, this)) {
					this.hierarchicalTree.gatherLeaves(leaves);
					for (TreeNode leaf : leaves)
						for (Double splt : reg.possibleSplitValues) {
							// try to split
							optModuleScore = leaf.trySplit(reg, splt);
							// System.out.println(test);
							if (optModuleScore - maxModuleScore > 0.0) {
								optLeaf = leaf;
								optReg = reg;
								optSplit = splt;
								maxModuleScore = optModuleScore;
							}
						}
				}
			}
		}

		// we have found a best split if max has changed its value in the
		// previous loops
		if (maxModuleScore - this.moduleScore > 0.0) {
			scoreGain = maxModuleScore - this.moduleScore;
			this.bestSplit = new Split(optLeaf, optReg, optSplit, 0.0, scoreGain);
			System.out.println("... best split found with score gain "
					+ scoreGain / this.moduleNetwork.numGenes);
			return true;
		} else {
			this.bestSplit = new Split();
			System.out.println("... no best split found");
			return false;
		}
	}
	

	

	/**
	 * Computes the best split of the merged leaves below a node after a
	 * reassignment.
	 * 
	 * @param node node for which split is searched.
	 * 
	 * @author tomic
	 * 
	 */
	public Split bestMergeSplit(TreeNode node) {
		if (node.nodeStatus == "leaf")
			return new Split();
		// get list of leaves below node
		ArrayList<TreeNode> leafList = new ArrayList<TreeNode>();
		node.gatherLeafList(leafList);
		// create fictitious node to hold all leaves, note that
		// statistics are additive
		TreeNode fictnode = new TreeNode(this, node.leafDistribution.normalGammaPrior);
		for (TreeNode leaf : leafList) {
			fictnode.leafDistribution.statistics[0] += leaf.leafDistribution.statistics[0];
			fictnode.leafDistribution.statistics[1] += leaf.leafDistribution.statistics[1];
			fictnode.leafDistribution.statistics[2] += leaf.leafDistribution.statistics[2];
			for (int m : leaf.leafDistribution.condSet)
				fictnode.leafDistribution.condSet.add(m);
		}
		// compute fictitious node score: this is the unsplitted score
		fictnode.leafDistribution.bayesianScore();
		// find the best split of fictnode, set initially to current
		double optSplit = node.regulationSplit.splitValue;
		Gene optReg = node.regulationSplit.regulator;
		// score of current split
		double scoreGain = fictnode.trySplit(optReg, optSplit);
		double optScoreGain;
		// first try other split values with optReg
		for (Double splt : optReg.possibleSplitValues) {
			if (splt != node.regulationSplit.splitValue) {
				optScoreGain = fictnode.trySplit(optReg, splt);
				if (optScoreGain - scoreGain > 0.0) {
					optSplit = splt;
					scoreGain = optScoreGain;
				}
			}
		}
		// try other regulators that preserve acyclicity
		for (Gene reg : this.moduleNetwork.regulatorSet) {
			if (reg.number != optReg.number
					&& !this.moduleNetwork.moduleGraph.testCycle(reg, this)) {
				for (Double splt : reg.possibleSplitValues) {
					optScoreGain = fictnode.trySplit(reg, splt);
					if (optScoreGain - scoreGain > 0.0) {
						optReg = reg;
						optSplit = splt;
						scoreGain = optScoreGain;
					}
				}
			}
		}

		// there must be a best split, at the very least the current
		// one; note this is a split for fictnode
		return new Split(fictnode, optReg, optSplit, 0.0, scoreGain);
	}
	
	/**
	 * Trims a regulation tree by removing nodes that are no longer best splits
	 * after a reassignment. Resets mod.bestSplit if applicable
	 * 
	 * @param mod Module where trimming is performed.
	 * @param node node below which trimming is performed.
	 * 
	 * @author tomic
	 * 
	 */
	public void treeTrim(TreeNode node) {
		Split splt;
		if (node.nodeStatus == "internal") {
			System.out.println("Testing node " + node.regulationSplit.regulator.name);
			// find the best split
			splt = this.bestMergeSplit(node);
			// check if it is still the same, if yes go down in tree,
			// if not prune the tree, but don't apply the new split
			if (splt.regulator.number == node.regulationSplit.regulator.number
					&& splt.splitValue == node.regulationSplit.splitValue) {
				System.out.println("... keep");
				this.treeTrim(node.leftChild);
				this.treeTrim(node.rightChild);
			} else {
				System.out.println("... trim");
				// remove regulators below and including node from parent set
				// and modulegraph
				removeRegulators(node);
				// make it a leaf, note that all the information is allready in
				// the bestMergeSplit
				node.nodeStatus = "leaf";
				node.leftChild = null;
				node.rightChild = null;
				node.leafDistribution = splt.node.leafDistribution;
				// set splt.leaf to node so later split is applied in the right place
				splt.node= node;
				// check if the split would be a best split for the current module
				if (splt.scoreGain >= this.bestSplit.scoreGain)
					this.bestSplit = splt;
			}
		}
	}
	
	
	/**
	 * Removes all regulators in a subtree starting at a node.
	 * 
	 * @param node start node.
	 * 
	 * @author tomic
	 */
	public void removeRegulators(TreeNode node) {
		if (node.nodeStatus == "internal") {
			this.moduleNetwork.moduleGraph.removeRegulator(node.regulationSplit.regulator, node.module);
			node.module.parents.remove(node.regulationSplit.regulator);
			removeRegulators(node.leftChild);
			removeRegulators(node.rightChild);
		}
	}

}
