/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/



package lemontree.networktools;

import java.util.*;

import lemontree.modulenetwork.*;

/**
 * A class for weighted edges, mainly to sort edges in a network by their weight
 * 
 * @author tomic
 *
 */

public class Edge implements Comparable {
	
	public String source;
	
	public String target;
	
	public double weight;
	
	public int[] inModule;
	
	public int[] signs;

	public int sign;
	
    public ArrayList<Integer> condSet; 
	
	/**
	 * +1 if the edge is a tp, -1 if it is a fp, 0 if undefined (one of the nodes does not
	 * belong to the reference network)
	 */
	public int valid;
	
	public Edge(){
		
	}
	
	public Edge (String src, String tgt){
		this.source = src;
		this.target = tgt;
		this.weight = 1.0;
		this.valid = 0;
		this.inModule = new int[0];
		this.signs = new int[0];
		this.condSet = new ArrayList<Integer>();
	}
	
	public Edge (String src, String tgt, double w){
		this.source = src;
		this.target = tgt;
		this.weight = w;
		this.valid = 0;
		this.inModule = new int[0];
		this.signs = new int[0];
		this.condSet = new ArrayList<Integer>();
	}
	

	/**
	 * Compares 2 edges by minus their weight (we want them sorted from high to low)
	 * 
	 * @author tomic
	 */
	public int compareTo(Object o) {
        Double w1 = -this.weight;
        Double w2 = -((Edge)o).weight;
        return w1.compareTo(w2);
    }
	
	/**
	 * Edges are equal if they have the same source and target (for recall - precision computations)
	 * 
	 * @author tomic
	 */
	public boolean equals(Edge edge) {
		if ( edge.source.equals(this.source) && edge.target.equals(this.target)){
			return true;
		} else
			return false;
	}
	
	/**
	 * Creates string.
	 * 
	 * @author tomic
	 */
	public String toString(){
		StringBuilder sb = new StringBuilder();
		for (int k=0; k<this.inModule.length; k++){
			String sgn = "0";
			if (this.signs[k] == 1)
				sgn = "+";
			else if (this.signs[k] == -1)
				sgn = "-";
			String s = String.format("\t %3d(" + sgn + ")",	this.inModule[k], this.signs[k]);
			sb.append(s);
		}

		StringBuilder sb2 = new StringBuilder();
		for (int m : this.condSet){
		    sb2.append(String.format("\t%3d", m));
		}

		return String.format(this.source + "\t" + this.target + "\t %f \t %d" + sb.toString() + sb2.toString(), this.weight, this.valid);
	}
	
	/**
	 * Validates edge.
	 * 
	 * @author tomic
	 */
	public void validate (Network net){
		if (net.containsEdge(this))
			this.valid = 1;
		else if (net.edges.keySet().contains(this.source) && net.nodeList.contains(this.target))
			this.valid = -1;
		else 
			this.valid = 0;
	}

	/**
	 * Finds the rank of an edge in a network, store in field {@link valid}.
	 * @param net
	 * 
	 * @author tomic
	 */
	public void setRank(Network net){
		net.sortEdges();
		for (int k=0; k<net.edgeList.size(); k++){
			if (this.equals(net.edgeList.get(k))){
				this.valid = k;
				break;
			}
		}
	}
	
	
	/**
	 * Find modules in a module network to which the target belongs
	 * 
	 * @author tomic
	 */
	public void setInModule (ModuleNetwork M){
		this.inModule = new int[M.moduleSets.size()];
		this.signs = new int[M.moduleSets.size()];
		ArrayList<ArrayList<Integer>> allConds = new ArrayList<ArrayList<Integer>>();
		Gene tgtGene = M.geneMap.get(this.target);
		Gene srcGene = M.geneMap.get(this.source);
		for (int k=0; k<M.moduleSets.size(); k++){
			for (Module mod : M.moduleSets.get(k)){
				if (mod.genes.contains(tgtGene)){
					this.inModule[k] = mod.number;
					// find regulator sign
					if (mod.parentSigns.containsKey(srcGene)){
						this.signs[k] = mod.parentSigns.get(srcGene);
						for (TreeNode root : mod.hierarchicalTrees){
						    for (TreeNode node : root.getInternalNodes()){
							for (Split splt : node.testSplits){
							    if (splt.regulator.name.equals(this.source)){
								allConds.add(node.leafDistribution.condSet);
								break;
							    }
							}
						    }
						}
					}
					else
						this.signs[k] = 0;
					break;
				}
			}
		}
		this.setOverallSign();
	}

	/**
	 * Sets the overall sign.
	 * 
	 * @author tomic
	 */
	public void setOverallSign(){
		int sum = 0;
		for (int s : this.signs){
			sum += s;
		}
		if (sum > 0)
			this.sign = 1;
		else if (sum < 0)
			this.sign = -1;
		else
			this.sign = 0;
	}
	
}
