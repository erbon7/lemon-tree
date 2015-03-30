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


import java.io.*;
import java.util.*;

import cern.colt.list.DoubleArrayList;

/**
 * A simple class for weighted networks. By convention, in an undirected network each edge appears
 * only once, a boolean issym keeps track whether the edges are symmetric or not
 * 
 * @author tom
 *
 */

public class Network {


	public ArrayList<String> nodes; // map from integers to nodes
	
	public HashMap<String,Integer> nodeidx; // map from nodes to integers
	
	public ArrayList<Edge> edges;

	public boolean issym; // TRUE for symmetric (undirected) networks
	
	public JavaSparseArray A; // adjacency matrix

	
	public Network(){
		this.edges = new ArrayList<Edge>();
		this.setProperties();
	}
	
	public Network(boolean issym){
		this.edges = new ArrayList<Edge>();
		this.setProperties();
		this.issym = issym;
	}
	
	public Network(ArrayList<Edge> edges, boolean issym){
		this.edges = edges;
		this.issym = issym;
		this.setProperties();	
	}
	
	
	/**
	 * Read network from file
	 * @param dir
	 * @param file
	 */
	public Network (String dir, String file, boolean issym) {
		this.edges = new ArrayList<Edge>();
		try {
			Scanner fileScanner = new Scanner(new File(dir, file));
			fileScanner.useDelimiter("\\n");
			// walk through file, skip comment lines (starting with #)
			while (fileScanner.hasNext()){
				Scanner line = new Scanner(fileScanner.next().toString());
				line.useDelimiter("\\t");
				// read from
				String from = line.next();
				from = from.trim();
				// skip if comment line
				if (!from.startsWith("#")){
					//read tp
					String to = line.next();
					to = to.trim();
					if (line.hasNext()){
						Double w = Double.parseDouble(line.next());
						this.edges.add(new Edge(from,to,w));
					} else {
						this.edges.add(new Edge(from,to,1.0));
					}
				}
				line.close();
			}
			this.issym = issym;
			for (Edge e:this.edges){
				e.issym = issym;
			}
			this.setProperties();
			System.out.println("Read network with " + this.nodes.size() + " nodes and " + this.edges.size() + " edges.");
			fileScanner.close();
		} catch (FileNotFoundException e) {
			System.out.println(e);
		}
	}
	
	
	/**
	 * Create sparse adjacency matrix
	 */
	public void setMatrix(){
		// create list of neighbors for all nodes
		HashMap<String,HashSet<Edge>> neighborsOf = new HashMap<String,HashSet<Edge>>();
		for (Edge e: this.edges){
			if (!neighborsOf.containsKey(e.from))
				neighborsOf.put(e.from, new HashSet<Edge>());
			if (!neighborsOf.containsKey(e.to))
				neighborsOf.put(e.to, new HashSet<Edge>());
			neighborsOf.get(e.from).add(e);
			neighborsOf.get(e.to).add(e);
		}
		
		// fill index and value arrays
		int[][] indexes = new int[this.nodes.size()][];
		double[][] values = new double[this.nodes.size()][];
		for (int k=0; k<this.nodes.size(); k++){
			String nd = this.nodes.get(k);
			HashSet<Edge> nb = neighborsOf.get(nd);
			indexes[k] = new int[nb.size()+1];
			values[k] = new double[nb.size()+1];
			// set 1 on diagonal
			indexes[k][0] = k;
			values[k][0] = 1.0;
			int ctr=1;
			for (Edge e : nb){
				String other;
				if (e.from.equals(nd))
					other = e.to;
				else
					other = e.from;
				indexes[k][ctr] = this.nodeidx.get(other);
				values[k][ctr] = e.weight;
				ctr++;
			}
		}		
		this.A = new JavaSparseArray (values, indexes);
	}
	
	/**
	 * Set network properties
	 */
	public void setProperties (){
		this.setNodes();
		this.setMatrix();
	}
	
	/**
	 * Collects all nodes of the network in one set.
	 * We should actually sort them alphabetically
	 */
	public void setNodes () {
		HashSet<String> nds = new HashSet<String>(); // with hashset we don't have to worry about duplicates
		for (Edge e : this.edges){
			nds.add(e.from);
			nds.add(e.to);
		}
		this.nodes = new ArrayList<String>();
		for (String node : nds)
			this.nodes.add(node);
		// set inverse map
		this.nodeidx = new HashMap<String,Integer>();
		for (int k=0; k<this.nodes.size(); k++)
			this.nodeidx.put(this.nodes.get(k), k);
	}
	
	/**
	 * Remove edges between a subset of nodes from the network
	 * @param cluster
	 */
	public void removeEdges(HashSet<String> clusterX){
		ArrayList<Edge> tobeRemoved = new ArrayList<Edge>();
		for (Edge e: this.edges)
			if (clusterX.contains(e.from) && clusterX.contains(e.to))
				tobeRemoved.add(e);
		for (Edge e : tobeRemoved)
			this.edges.remove(e);
		this.setProperties();

	}
	
	/**
	 * Remove edges between a subset of nodes from the network
	 * @param cluster
	 */
	public void removeEdges(HashSet<String> clusterX, HashSet<String> clusterY){
		ArrayList<Edge> tobeRemoved = new ArrayList<Edge>();
		for (Edge e: this.edges)
			if (clusterX.contains(e.from) && clusterY.contains(e.to))
				tobeRemoved.add(e);
		for (Edge e : tobeRemoved)
			this.edges.remove(e);
		this.setProperties();

	}
	
	/**
	 * Remove edges below certain weight
	 * @param w
	 */
	public void removeEdges(double w){
		ArrayList<Edge> tobeRemoved = new ArrayList<Edge>();
		for (Edge e: this.edges)
			if (e.weight < w)
				tobeRemoved.add(e);
		for (Edge e : tobeRemoved)
			this.edges.remove(e);
		this.setProperties();
	}
	
	/**
	 * Remove a subset of nodes from the network
	 * @param cluster
	 */
	public void removeNodes(HashSet<String> clusterX){
		ArrayList<Edge> tobeRemoved = new ArrayList<Edge>();
		for (Edge e: this.edges)
			if (clusterX.contains(e.from) || clusterX.contains(e.to))
				tobeRemoved.add(e);
		for (Edge e : tobeRemoved)
			this.edges.remove(e);
		this.setProperties();

	}	
	
	/**
	 * Remove a subset of nodes from the network
	 * @param cluster
	 */
	public void removeNodes(HashSet<String> clusterX, HashSet<String> clusterY){
		ArrayList<Edge> tobeRemoved = new ArrayList<Edge>();
		for (Edge e: this.edges)
			if (clusterX.contains(e.from) || clusterY.contains(e.to))
				tobeRemoved.add(e);
		for (Edge e : tobeRemoved)
			this.edges.remove(e);
		this.setProperties();

	}	

	/**
	 * Add an edge or update weight of an existing edge with or without setting other properties
	 * @param e
	 */
	public void addEdge (Edge e, boolean setall){
		if (setall){
			boolean exists = false;
			for (Edge f : this.edges){
				if (f.equals(e)){
					f.weight += e.weight;
					exists = true;
					break;
				}
			}
			if (!exists)
				this.edges.add(e);
			this.setProperties();
		} else { // only use if we are sure the edge does not exist yet
			this.edges.add(e);
		}
	}
	
	/**
	 * Sum of edge weights for a certain node
	 * @param node
	 * @return
	 */
	public double totalWeight (String node){
		int ix = this.nodeidx.get(node);
		double w = 0.0;
		for (double w2 : this.A.getValueArray()[ix])
			w += w2;
		return w;
	}
	
	/**
	 * Score of a cluster = number of edges / number of nodes
	 * @param clust
	 * @return
	 */
	public double clusterScore(HashSet<String> clust){
		double[] v = new double[this.nodes.size()];
		for (int k=0; k<v.length; k++)
			if (clust.contains(this.nodes.get(k)))
				v[k]=1.0;
			else
				v[k]=0.0;
		double[] w = this.A.matrixvector(v);
		double score = 0.0;
		for (int k=0; k<v.length; k++)
			score += v[k]*w[k];
//		double score = 0.0;
//		for (Edge e : this.edges)
//			if (clust.contains(e.from) && clust.contains(e.to))
//				score += e.weight;
		score = score/(double)clust.size();
		return score;
	}
	
	/**
	 * Score of a cluster = number of edges / number of nodes
	 * @param v vector of 0 and 1's
	 * @return
	 */
	public double clusterScore(double[] v){
		double[] w = this.A.matrixvector(v);
		double score = 0.0;
		for (int k=0; k<v.length; k++)
			score += v[k]*w[k];
		score = score/Math.pow(norm(v),2);
		return score;
	}
	
	/**
	 * Compute Perron vector for an adjacency matrix defined by a list of weighted edges
	 * Perron vector is the dominant eigenvector for a symmetric nonnegative matrix. The current code does
	 * not support asymmetric matrices.
	 * @param edgelist
	 * @param p
	 * @return
	 */
	public double[] perronvector (double tolerance, int maxstep){
		// Initialize v to node with highest total degree
		double dmax = 0;
		String nodemax = new String();
		for (String node : this.nodes){
			double nb = this.totalWeight(node);
			if (nb>dmax){
				dmax = nb;
				nodemax = node;
			}
		}
		double[] v = new double[this.nodes.size()];
		v[this.nodeidx.get(nodemax)] = 1.0;
		double mu = 1.0, munew;

		// Power algorithm until convergence
		int step=0;
		double diff=1;
		while (diff>tolerance && step<maxstep){
			step++;
			// Matrix - vector multiplications
			v = this.A.matrixvector(v);
			munew = norm(v);
			// normalize
			normalize(v, munew);
			// Convergence parameter
			diff = Math.abs(1.0 - munew/mu);
			// Update mu
			mu = munew;
		}
		if (step>=maxstep)
			System.out.println("Max. number of steps reached, diff = " + diff);
		return v;
	}

	/**
	 * Create best cluster from vector weights
	 * @param v
	 * @param tol
	 * @return
	 */
	public HashSet<Integer> bestCluster (double[] v, double tolerance){
		// get list of values above tolerance level
		DoubleArrayList vals = new DoubleArrayList();
		for (int k=0; k<v.length; k++)
			if (v[k]>tolerance){
				vals.add(v[k]);
			}
		// sort in descending order
		vals.quickSort();
		vals.reverse();

		// compute scores
		double[] elements = new double[this.nodes.size()];
		double maxscore=0.0, maxval=0.0, newscore;
		for (int k=0; k<vals.size(); k++){
			for (int m=0; m<v.length; m++){
				if (v[m]>=vals.get(k))
					elements[m] = 1.0;
			}
			newscore = this.clusterScore(elements);
			if (newscore >= maxscore){
				maxscore = newscore;
				maxval = vals.get(k);
			}
		}
		// set final cluster
		HashSet<Integer> cluster = new HashSet<Integer>();
		for (int m=0; m<v.length; m++){
			if (v[m]>=maxval)
				cluster.add(m);
		}
		
		return cluster;
	}
	
	/**
	 * Convert weight vector to node set
	 * @param v
	 * @param p
	 * @return
	 */
	public static HashSet<Integer> unifapprox (double[] v, double tolerance){
		// get list of values above tolerance level
		DoubleArrayList vals = new DoubleArrayList();
		for (int k=0; k<v.length; k++)
			if (v[k]>tolerance){
				vals.add(v[k]);
			}
		// sort in descending order
		vals.quickSort();
		vals.reverse();

		// initialize uniform approx vector
		double[] v1 = new double[vals.size()];
		v1[0] = 1.0;
		for (int k=1; k<v1.length; k++)
			v1[k] = 0.0;
		double[] diff = new double[vals.size()];
		for (int k=0; k<diff.length; k++)
			diff[k] = vals.get(k)-v1[k];
		double minValue=norm(diff);
		int minCtr=1;

		for(int ctr=2;ctr<=vals.size();ctr++){
			// update v en diff
			for (int k=0; k<ctr; k++)
				v1[k] = 1.0;
			normalize(v1);
			for (int k=0; k<diff.length; k++)
				diff[k] = vals.get(k)-v1[k];
			double nValue=norm(diff);
			if (nValue<minValue){
				minValue=nValue;
				minCtr=ctr;
			}
		}
		double cut = vals.get(minCtr-1);
		HashSet<Integer> clusterI = new HashSet<Integer>();
		for (int k=0; k<v.length; k++)
			if (v[k]>=cut){
				clusterI.add(k);
			}
		
		return clusterI;

	}
	
	/**
	 * Normalize vector to have norm 1
	 * @param v
	 * @param p
	 */
	public static void normalize (double[] v, double nrm){
		// normalize
		for (int k=0; k<v.length; k++)
			v[k] = v[k]/nrm;
	}
	
	/**
	 * Normalize vector to have norm 1
	 * @param v
	 * @param p
	 */
	public static void normalize (double[] v){
		// compute norm
		double nrm = norm(v);
		// normalize
		for (int k=0; k<v.length; k++)
			v[k] = v[k]/nrm;
	}
	
	/**
	 * Norm of a vector
	 * @param v
	 * @param p
	 * @return
	 */
	public static double norm (double[] v){
		double nrm = 0.0;
		for (double x : v)
			nrm += Math.pow(x,2);
		nrm = Math.sqrt(nrm);
		return nrm;
	}

	
	/**
	 * Convert set of node indices to set of node names
	 * @param clI
	 * @return
	 */
	public HashSet<String> idxToNames(HashSet<Integer> clI){
		HashSet<String> cl = new HashSet<String>();
		for (Integer k : clI)
			cl.add(this.nodes.get(k));
		return cl;
	}

	
}
