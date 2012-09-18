/*
 * PFClustering - Java package for clustering weighted networks based on the Perron-Frobenius Theorem
 * 
 * Copyright (c) 2011 Tom Michoel
 *
 * This is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 * 
 */


package lemontree.pfclustering;



import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class PFClustering {

	/**
	 * Provides the EdgeClust algorithm
	 */
	
	public Network network;
	
	public ArrayList<HashSet<String>> clustersX;
	// public ArrayList<HashSet<String>> clustersY; // for asymmetric networks, not yet supported
		
	public double tolerance; // convergence tolerance, also used in unifapprox 
	
	public int maxstep; // max number of iterations in perronvector calculations
	
	public int minClustSize; // minimal number of nodes in a cluster
	
	public double minClustScore; // minimal score of a cluster

	public PFClustering (){
		
	}
	
	public PFClustering (Network net){
		this.network = net;
		// set default parameters
		this.tolerance = 1E-5;
		this.maxstep = 1000;
		this.minClustSize = 1;
		this.minClustScore = 0.0;
	}
	
	public PFClustering (Network net, double tolerance, int maxstep, int minClustSize, double minClustScore){
		this.network = net;
		this.tolerance = tolerance;
		this.maxstep = maxstep;
		this.minClustSize = minClustSize;
		this.minClustScore = minClustScore;
	}
	
	/**
	 * Perform nodeClustering algorithm
	 * @param p
	 * @return
	 */
	public void nodeClustering (){
		this.pfClustering(true);
	}
	
	/**
	 * Perform edgeClustering algorithm
	 * @param p
	 * @return
	 */
	public void edgeClustering (){
		this.pfClustering(false);
	}
	
	/**
	 * Perform Perron-Frobenius based node or edge clustering
	 * @param p
	 * @param node
	 */
	public void pfClustering(boolean node) {
		if (!this.network.issym)
			System.out.println("Error: pfClustering needs symmetric network.");
		this.clustersX = new ArrayList<HashSet<String>>();
		double[] v;
		double score;
		// create deep copy of network
		ArrayList<Edge> edges = new ArrayList<Edge>(this.network.edges.size());
		for (Edge e : this.network.edges)
			edges.add(e);
		Network remainingNet = new Network(edges,this.network.issym); 
		// clustering loop
		while (!remainingNet.edges.isEmpty()){
			// reset v,w
			v = remainingNet.perronvector(this.tolerance, this.maxstep);
			// HashSet<String> clusterX = remainingNet.idxToNames(remainingNet.unifapprox(v, this.tolerance));
			HashSet<String> clusterX = remainingNet.idxToNames(remainingNet.bestCluster(v, this.tolerance));
			clustersX.add(clusterX);
			score = this.network.clusterScore(clusterX);
			if (node)
				remainingNet.removeNodes(clusterX);
			else
				remainingNet.removeEdges(clusterX);
			System.out.println("Found cluster with score " + score + " and " + clusterX.size() + " nodes. Remaining network has " + remainingNet.nodes.size()
					+ " nodes and "+ remainingNet.edges.size() + " edges.");
		}
		
		// post-processing
		this.selectClustersSize();
		this.selectClustersScore();
		System.out.println("Found " + this.clustersX.size() + " clusters with at least " + this.minClustSize + " nodes and score higher than " + this.minClustScore + ".");
	}
	

	
	/**
	 * 
	 */
	public void selectClustersSize(){
		HashSet<HashSet<String>> remove = new HashSet<HashSet<String>>();
		for (HashSet<String> clust : this.clustersX)
			if (clust.size()<this.minClustSize)
				remove.add(clust);
		for (HashSet<String> clust : remove)
			this.clustersX.remove(clust);
	}
	
	/**
	 * 
	 */
	public void selectClustersScore(){
		HashSet<HashSet<String>> remove = new HashSet<HashSet<String>>();
		for (HashSet<String> clust : this.clustersX)
			if (this.network.clusterScore(clust) < this.minClustScore)
				remove.add(clust);
		for (HashSet<String> clust : remove)
			this.clustersX.remove(clust);
	}
	
	/**
	 * Write clusters to file
	 * @param dir
	 * @param file
	 */
	public void writeClustersX(String dir, String file){
		try {
			File f = new File(dir, file);
			PrintWriter pw = new PrintWriter(f);
			int ctr = 0;
			for (HashSet<String> cluster : this.clustersX) {
				for (String node : cluster){
					String s = String.format("%s\t%d", node, ctr); 
					pw.println(s);
				}
				ctr++;
			}
			pw.close();
		} catch (IOException e) {
			System.out.println(e);
		}
	}
	
}
