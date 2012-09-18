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
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import lemontree.networktools.*;

/**
 * Class to perform various calculations using the probabilistic models
 * in a module networks ensemble.
 * 
 * @author tomic
 *
 */

public class EnsembleCalculator {

	/**
	 * Module network containing all necessary information
	 */
	ModuleNetwork moduleNetwork;
	
    
	public EnsembleCalculator(){
		
	}
	
	public EnsembleCalculator(ModuleNetwork M){
		this.moduleNetwork = M;
	}
	
    /**
     * Get one column of the data matrix
     * @param m
     * @return
     * 
     * @author tomic
     */
    public double[] getDataColumn(int m){
    	double[] x = new double[this.moduleNetwork.geneSet.size()];
    	for (int k=0; k<x.length; k++)
    		x[k] = this.moduleNetwork.data[k][m];
    	return x;
    }
 
    public void sampleModuleNetwork(){
    	Random rand = new Random();
    	// sample module set
		this.moduleNetwork.moduleSet = this.moduleNetwork.moduleSets.get(rand.nextInt(this.moduleNetwork.moduleSets.size()));
		for (Module mod : this.moduleNetwork.moduleSet){
			// sample tree
			mod.hierarchicalTree = mod.hierarchicalTrees.get(rand.nextInt(mod.hierarchicalTrees.size()));
			for (TreeNode node : mod.hierarchicalTree.getInternalNodes()){
				// sample split
				node.regulationSplit = node.testSplits.get(rand.nextInt(node.testSplits.size()));
			}
		}
    }
    
    public double logProbData(){
    	double p = 0.0;
    	for (int m=0; m<this.moduleNetwork.data[0].length; m++){
    		double[] col = this.getDataColumn(m);
    		for (Module mod : this.moduleNetwork.moduleSet){
    			for (Gene gene : mod.genes){
    				double pi = 0.0;
    				for (TreeNode leaf : mod.hierarchicalTree.getLeafNodes()){
    					double w = leaf.pathWeight(col);
    					pi += w * leaf.leafDistribution.evaluate(col[gene.number]);
    				}
    				p += Math.log(pi);
    			}
    		}
    	}
    	return p;
    }
    
    
    
    /**
     * Evaluates a single module network model on a new experiment 
     * @param xExp expression value for each gene
     * @return log of probability distribution evaluated at xExp
     * 
     * @author tomic
     */
    public double evaluateModel (double[] xExp){
    	if (xExp.length != this.moduleNetwork.geneSet.size())
    		return Double.NaN;
    	else {
    		double p = 0.0;
    		for (Module mod : this.moduleNetwork.moduleSet){
    			for (Gene gene : mod.genes){
					double pi = 0.0;
					for (TreeNode leaf : mod.hierarchicalTree.getLeafNodes()){
						double w = leaf.pathWeight(xExp);
						pi += w * leaf.leafDistribution.evaluate(xExp[gene.number]);
					}
					p += Math.log(pi);
				}
    		}
    		System.out.println(p/xExp.length);
    		return p/xExp.length;
    	}
    }
    

    /**
     * Evaluates the average of a given number of module network models sampled from the ensemble
     * on a new experiment. 
     * @param xExp expression value for each gene: genes not included have xExp=NaN
     * @return log of probability distribution evaluated at xExp
     * 
     * @author tomic
     */
    public double[] evaluateModelEnsemble (double[] xExp, int numSample){
    	Random rand = new Random();
    	double[] pEns = new double[numSample];
    	for (int k=0; k<numSample; k++){
    		// sample module set
    		this.moduleNetwork.moduleSet = this.moduleNetwork.moduleSets.get(rand.nextInt(this.moduleNetwork.moduleSets.size()));
    		for (Module mod : this.moduleNetwork.moduleSet){
    			// sample tree
    			mod.hierarchicalTree = mod.hierarchicalTrees.get(rand.nextInt(mod.hierarchicalTrees.size()));
    			for (TreeNode node : mod.hierarchicalTree.getInternalNodes()){
    				// sample split
    				node.regulationSplit = node.testSplits.get(rand.nextInt(node.testSplits.size()));
    			}
    		}
    		// evaluate sampled model and add to average
    		pEns[k] = this.evaluateModel(xExp);
    	}
    	return pEns;
    }
    
    /**
     * Evaluates a single module network model on a new experiment
     * for a subset of genes.
     * @param xExp expression values for a subset of genes: genes not included have xExp=NaN
     * @return log of probability distribution evaluated at xExp
     * 
     * @author tomic
     */
    public double evaluateModelSub (double[] xExp){
    	double p = 0.0;
    	int num = 0;
    	for (Module mod : this.moduleNetwork.moduleSet){
    		if (mod.hierarchicalTree.nodeStatus.equals("internal")){
    			mod.parents = new ArrayList<Gene>();
    			for (TreeNode node : mod.hierarchicalTree.getInternalNodes())
    				mod.parents.add(node.regulationSplit.regulator);
    			// all parents need to be in subset of genes
    			boolean b = true;
    			for (Gene reg : mod.parents)
    				b = b & !Double.isNaN(xExp[reg.number]);
    			// if true we can compute probabilities
    			if (b){
    				for (Gene gene : mod.genes){
    					if (!Double.isNaN(xExp[gene.number])){
    						double pi = 0.0;
    						for (TreeNode leaf : mod.hierarchicalTree.getLeafNodes()){
    							double w = leaf.pathWeight(xExp);
    							pi += w * leaf.leafDistribution.evaluate(xExp[gene.number]);
    						}
    						p += Math.log(pi);
    						num += 1;
    					}
    				}
    			}
    		}
    	}
    	return p/num;
    }
    
    /**
     * Evaluates the average of a given number of module network models sampled from the ensemble
     * on a new experiment for a subset of genes. 
     * @param xExp expression value for each gene
     * @return log of probability distribution evaluated at xExp
     * 
     * @author tomic
     */
    public double[] evaluateModelEnsembleSub (double[] xExp, int numSample){
    	Random rand = new Random();
    	double[] pEns = new double[numSample];
    	for (int k=0; k<numSample; k++){
    		// sample module set
    		this.moduleNetwork.moduleSet = this.moduleNetwork.moduleSets.get(rand.nextInt(this.moduleNetwork.moduleSets.size()));
    		for (Module mod : this.moduleNetwork.moduleSet){
    			// sample tree
    			mod.hierarchicalTree = mod.hierarchicalTrees.get(rand.nextInt(mod.hierarchicalTrees.size()));
    			for (TreeNode node : mod.hierarchicalTree.getInternalNodes()){
    				// sample split
    				node.regulationSplit = node.testSplits.get(rand.nextInt(node.testSplits.size()));
    			}
    		}
    		// evaluate sampled model and add to average
    		pEns[k] = this.evaluateModelSub(xExp);
    	}
    	return pEns;
    }
 
   
    
    public void computeEdgeProbs(int numSample){
    	Network net = new Network();
    	Random rand = new Random();
    	HashMap<String, HashMap<String,boolean[]>> edges = new HashMap<String, HashMap<String,boolean[]>>();
    	
    	double Z = 0.0;
    	
    	// first take sample to normalize score
    	this.sampleModuleNetwork();
    	double p0 = this.logProbData();

    	// next take samples to compute probabilities
    	for (int k=0; k<numSample; k++){
    		// sample module network
    		this.sampleModuleNetwork();
    		System.out.println(k);
    		for (Module mod : this.moduleNetwork.moduleSet){
    			mod.parents = new ArrayList<Gene>();
    			for (TreeNode node : mod.hierarchicalTree.getInternalNodes()){
    				mod.parents.add(node.regulationSplit.regulator);
    			}
    			for (Gene reg : mod.parents){
    				if (!net.edges.containsKey(reg.name))
    					net.edges.put(reg.name, new HashMap<String,Double>());
    				for (Gene gene : mod.genes){
    					if (!net.edges.get(reg.name).containsKey(gene.name))
    						net.edges.get(reg.name).put(gene.name, 0.0);
    					double w = net.edges.get(reg.name).get(gene.name);
    					net.edges.get(reg.name).put(gene.name, w + 1.0);
    				}
    			}
    		}
    		
    		
    		// log-probability of sampled model
    		double p = this.logProbData();
    		System.out.println((p-p0)/p0);
    		// update normalization constant (partition function)
    		Z += Math.exp(p-p0);
    	}
    	net.scaleEdgeWeights((double)numSample);
    	net.setProperties();
    	net.printEdges("test.txt", 1000);
    	
     }
	
 }
