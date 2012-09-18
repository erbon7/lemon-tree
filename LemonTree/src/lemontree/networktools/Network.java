/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.networktools;

import java.io.*;
import java.util.*;

import lemontree.utils.SetOperation;
import cern.colt.matrix.impl.*;
import lemontree.modulenetwork.*;

/**
 * A class for networks between nodes of type String.
 * 
 * @author tomic
 *
 */

public class Network {

	/** 
	 * For each regulator the set of target genes and their weight.
	 */
	public HashMap<String, HashMap<String, Double>> edges; 
	
	/**
	 * For each gene the set of regulators and their weight
	 */
	public HashMap<String, HashMap<String, Double>> inverseEdges;
	
	/**
	 * For each gene the set of signed regulators. A signed regulator
	 * is "regname+", "regname-" or "regname0" for activation, repression 
	 * or unknown
	 */
	public HashMap<String, HashSet<String>> signedInverseEdges;

	/** 
	 * Number of nodes in network.
	 */
	public int numNodes;
	
	/**
	 * Number of edges in network.
	 */
	public int numEdges;
	
	/**
	 * List of all nodes in network.
	 */
	public HashSet<String> nodeList;
	
	/**
	 * List of all edges, sorted by weight
	 */
	public List<Edge> edgeList;
	
	/**
	 * List of lists of genes who share the same sources.
	 */
	public ArrayList<ArrayList<String>> commonSources;
	
	/**
	 * List of modules in the network. Modules are defined as sets of nodes
	 * sharing the same signed regulators.
	 */
	public ArrayList<ArrayList<String>> modules;
	

	/**
	 * Empty constructor.
	 * 
	 * @author tomic
	 */
	public Network () {
		this.edges = new HashMap<String, HashMap<String, Double>>();
		this.inverseEdges = new HashMap<String, HashMap<String, Double>>();
		this.setProperties();
	}
	
	/**
	 * Construct network from a set of edges.
	 * 
	 * @param edges
	 * 
	 * @author tomic
	 */
	public Network (HashMap<String, HashMap<String, Double>> edges) {
		this.edges = edges;
		this.setProperties();
	}
	
	/**
	 * Construct network from a file.
	 * 
	 * @param dir directory where file is located.
	 * @param file name of file.
	 * 
	 * @author tomic
	 */
	public Network (String dir, String file) {
		this.edges = new HashMap<String, HashMap<String, Double>>();
		try {
			Scanner fileScanner = new Scanner(new File(dir, file)).useDelimiter("\\n");
			// walk through file, skip comment lines (starting with #)
			while (fileScanner.hasNext()){
				Scanner line = new Scanner(fileScanner.next().toString()).useDelimiter("\\t");
				// read regulator
				String regname = line.next();
				regname = regname.trim();
				// skip if comment line
				if (!regname.startsWith("#")){
					//read gene name
					String genename = line.next();
					genename = genename.trim();
					if (!this.edges.keySet().contains(regname))
						this.edges.put(regname, new HashMap<String, Double>());
					if (line.hasNext()){
						Double w = Double.parseDouble(line.next());
						this.edges.get(regname).put(genename, w);
					} else {
						this.edges.get(regname).put(genename, 1.0);
					}
				}
			}
			this.setProperties();
		} catch (FileNotFoundException e) {
			System.out.println(e);
		}
	}
	
	/**
	 * Constructs a network from different network files stored in the same directory. Edge weights
	 * are summed over the different networks.
	 * 
	 * @param netDir
	 */
	public Network(String netDir){
		// get all files in netDir
		File dir = new File(netDir);
		File[] files = dir.listFiles();
		// initialize with 1st file
		System.out.println("Start with " + files[0].getName());
		Network net = new Network(netDir, files[0].getName());
		// add other files
		for (int k=1; k<files.length; k++){
			System.out.println("... adding " + files[k].getName());
			Network net2 = new Network(netDir, files[k].getName());
			net.mergeNet(net2);
		}
		this.edges = net.edges;
		this.setProperties();
	}
	
	/**
	 * Merges this network with another one, summing over the edge weights
	 * @param net
	 */
	public void mergeNet(Network net){
		for (String reg : net.edges.keySet()){
			if (this.edges.containsKey(reg)){
				for (String tgt : net.edges.get(reg).keySet()){
					double w1 = net.edges.get(reg).get(tgt);
					if (this.edges.get(reg).containsKey(tgt)){
						double w0 = this.edges.get(reg).get(tgt);
						this.edges.get(reg).put(tgt, w0 + w1);
					} else {
						this.edges.get(reg).put(tgt, w1);
					}
				}
			} else {
				this.edges.put(reg, net.edges.get(reg));
			}
		}
		this.setProperties();
	}
	
	/**
	 * Sets network properties.
	 * 
	 * @author tomic
	 */
	public void setProperties () {
		this.computeNodeList();
		this.computeNumEdges();
		this.setInverseEdges();
		this.setEdgeList();
	}

	/**
	 * Collects all nodes of the network in one set.
	 * 
	 * @author tomic
	 */
	public void computeNodeList () {
		this.nodeList = new HashSet<String>();
		for (String reg : this.edges.keySet()){
			this.nodeList.add(reg);
			for (String gene : this.edges.get(reg).keySet())
				this.nodeList.add(gene);
		}
		this.numNodes = this.nodeList.size();
	}
	
	/**
	 * Counts the number of edges in the network.
	 * 
	 * @author tomic
	 */
	public void computeNumEdges(){
		int num = 0;
		for (String reg : this.edges.keySet())
			num += this.edges.get(reg).keySet().size();
		this.numEdges = num;
	}
	
	/**
	 * Computes the highest weight in the network.
	 * @return
	 * 
	 * @author tomic
	 */
	public double maxWeight(){
		double max=0;
		for (String src : this.edges.keySet())
			for (String tgt: this.edges.get(src).keySet())
				if (this.edges.get(src).get(tgt) > max)
					max = this.edges.get(src).get(tgt);
    	return max;
	}
	
	/**
	 * Computes the lowest weight in the network.
	 * @return
	 * 
	 * @author tomic
	 */
	public double minWeight(){
		double min=this.maxWeight();
		for (String src : this.edges.keySet())
			for (String tgt: this.edges.get(src).keySet())
				if (this.edges.get(src).get(tgt) < min)
					min = this.edges.get(src).get(tgt);
    	return min;
	}
	
	/**
	 * Scales the weight of the edges with a given factor.
	 * 
	 * @author tomic
	 */
	public void scaleEdgeWeights(double scale){
		for (String src : this.edges.keySet())
			for (String tgt : this.edges.get(src).keySet()){
				double w = this.edges.get(src).get(tgt);
				this.edges.get(src).put(tgt, w/scale);
			}
		this.setProperties();
	}
	
	
	/**
	 * Sets the edge list.
	 * 
	 * @author tomic
	 */
	public void setEdgeList(){
		this.edgeList = new ArrayList<Edge>();
		for (String src : this.edges.keySet())
			for (String tgt : this.edges.get(src).keySet()){
				double w = this.edges.get(src).get(tgt);
				this.edgeList.add(new Edge(src, tgt, w));
			}
	}
	
	/**
	 * Sorts the edges by their weight.
	 *
	 * @author tomic
	 */
	public void sortEdges(){
		Collections.sort(this.edgeList);
	}
	
	/**
	 * Checks if the edge list contains a certain edge. Only the source and
	 * target node are checked, not the weight.
	 * 
	 * @param e
	 * @return
	 * 
	 * @author tomic
	 */
	public boolean containsEdge(Edge e){
		boolean b = false;
		for (int k=0; k<this.edgeList.size(); k++){
			if (this.edgeList.get(k).equals(e)){
				b = true;
				break;
			}
		}
		return b;
	}
	
	
	/**
	 * Returns network with given number of highest weight edges.
	 *
	 * @param rank
	 *
	 * @author tomic
	 */
	public Network getHighestRankEdges(int rank){
		Network net = new Network();
		this.sortEdges();
		for (int k=0; k<rank; k++){
			String src = this.edgeList.get(k).source;
			String tgt = this.edgeList.get(k).target;
			double w = this.edgeList.get(k).weight;
			if (!net.edges.containsKey(src))
				net.edges.put(src, new HashMap<String, Double>());
			net.edges.get(src).put(tgt, w);
		}
		net.setProperties();
		return net;
	}
	
	/**
	 * Returns network with given number of highest weight edges, and at
	 * most one regulator per gene
	 *
	 * @param rank
	 *
	 * @author tomic
	 */
	public Network getHighestRankEdgesBestPerGene(int rank){
		Network net = new Network();
		this.sortEdges();
		ArrayList<String> assignedTargets = new ArrayList<String>();
		int k = 0;
		int numTgt = 0;
		while (numTgt < rank){
			String src = this.edgeList.get(k).source;
			String tgt = this.edgeList.get(k).target;
			if (!assignedTargets.contains(tgt)){
				assignedTargets.add(tgt);
				double w = this.edgeList.get(k).weight;
				if (!net.edges.containsKey(src))
					net.edges.put(src, new HashMap<String, Double>());
				net.edges.get(src).put(tgt, w);
				numTgt += 1;
			}
			k += 1;
		}
		net.setProperties();
		return net;
	}
	
	/**
	 * Sets the inverse edges.
	 * 
	 * @author tomic
	 */
	public void setInverseEdges(){
		this.inverseEdges = new HashMap<String, HashMap<String, Double>>();
		for (String node : this.nodeList){
			this.inverseEdges.put(node, new HashMap<String, Double>());
			for (String source : this.edges.keySet())
				if (this.edges.get(source).containsKey(node)){
					double w = this.edges.get(source).get(node);
					this.inverseEdges.get(node).put(source, w);
				}
		}
		// some nodes may not have any regulators
		HashSet<String> removeNodes = new HashSet<String>();
		for (String node : this.inverseEdges.keySet())
			if (this.inverseEdges.get(node).isEmpty())
				removeNodes.add(node);
		for (String node : removeNodes)
			this.inverseEdges.remove(node);
	}
	
	
	/**
	 * Sets the signed inverse edges.
	 * 
	 * @param M
	 * 
	 * @author tomic
	 */
	public void setSignedInverseEdges(ModuleNetwork M){
		this.signedInverseEdges = new HashMap<String, HashSet<String>>();
		this.setInModuleEdges(M);
		for (Edge e : this.edgeList){
			if (!this.signedInverseEdges.containsKey(e.target))
				this.signedInverseEdges.put(e.target, new HashSet<String>());
			String regname;
			if (e.sign==1)
				regname = e.source + "+";
			else if (e.sign==-1)
				regname = e.source + "-";
			else
				regname = e.source + "0";
			this.signedInverseEdges.get(e.target).add(regname);
		}
	}
	
	/**
	 * Sets the edges if the inverse edges are known.
	 * 
	 * @author tomic
	 */
	public void setEdgesFromInverse(){
		this.edges = new HashMap<String, HashMap<String, Double>>();
		for (String tgt : this.inverseEdges.keySet()){
			for (String src : this.inverseEdges.get(tgt).keySet()){
				// check if it's already there
				if (!this.edges.keySet().contains(src))
					this.edges.put(src, new HashMap<String,Double>());
				double w = this.inverseEdges.get(tgt).get(src);
				this.edges.get(src).put(tgt, w);
			}
		}
		this.setProperties();
	}
	
	/**
	 * Finds the common edges between the first rank edges of this
	 * network and a second network.
	 * 
	 * @param net
	 * @param rank
	 * @return
	 * 
	 * @author tomic
	 */
	public ArrayList<Edge> commonEdges(Network net, int rank){
		Network net1 = this.getHighestRankEdges(1);
		Network net2 = net.getHighestRankEdges(1);

		ArrayList<Edge> commonEdges = new ArrayList<Edge>();

		if (net2.containsEdge(net1.edgeList.get(0))){
			commonEdges.add(net1.edgeList.get(0));
		}
		for (int k=1; k<rank; k++){
			Edge e1 = this.edgeList.get(k);
			net1.edgeList.add(e1);
			Edge e2 = net.edgeList.get(k);
			net2.edgeList.add(e2);
			if (net2.containsEdge(e1)){
				commonEdges.add(e1);
			}
			if (net1.containsEdge(e2)){
				for (int l=0; l<net1.edgeList.size(); l++){
					if (net1.edgeList.get(l).equals(e2)){
						commonEdges.add(net1.edgeList.get(l));
						break;
					}
				}
			}
		}
		return commonEdges;
	}
	
	/**
	 * Removes edges with weight below certain value.
	 * 
	 * @param min
	 * 
	 * @author tomic
	 */
	public void removeEdges(double min){
		// collect source nodes that have no targets left
		HashSet<String> removeSources = new HashSet<String>();
		for (String source : this.edges.keySet()){
			// collect targets of source that need to be removed
			HashSet<String> removeTargets = new HashSet<String>();
			for (String target : this.edges.get(source).keySet())
				if (this.edges.get(source).get(target) < min)
					removeTargets.add(target);
			// remove them
			for (String target : removeTargets)
				this.edges.get(source).remove(target);
			// check if there are any left for this source
			if (this.edges.get(source).isEmpty())
				removeSources.add(source);
		}
		// remove sources
		for (String source : removeSources)
			this.edges.remove(source);
		// reset network properties
		this.setProperties();
	}
	
	/**
	 * Removes self-interactions from the network. Useful if the reference
	 * net does not contain those interactions.
	 *
	 * @author tomic
	 */
	public void removeSelfInteractions(){
		for (String src : this.edges.keySet())
			this.edges.get(src).remove(src);
		this.setProperties();
	}
	
	/**
	 * Removes regulators with only min or less targets.
	 * @param min
	 * 
	 * @author tomic
	 */
	public void removeSingletons(int min){
		HashSet<String> remove = new HashSet<String>();
		for (String src : this.edges.keySet())
			if (this.edges.get(src).size() <= min)
				remove.add(src);
		for (String src : remove)
			this.edges.remove(src);
		this.setProperties();
	}
	
	/**
	 * Removes edges if the source is not in the regulator list of M and/or
	 * the target is not in the gene list of M.
	 * 
	 * @param M
	 */
	public void removeNotInModuleNetwork(ModuleNetwork M){
		System.out.println(this.numEdges + "\t" + this.numNodes);
		HashSet<String> remove = new HashSet<String>();
		for (String src : this.edges.keySet()){
			if (!M.regulatorSet.contains(M.geneMap.get(src)))
				remove.add(src);
			else {
				HashSet<String> non = new HashSet<String>();
				for (String tgt : this.edges.get(src).keySet()){
					if (!M.geneSet.contains(M.geneMap.get(tgt)))
						non.add(tgt);
				}
				for (String tgt : non)
					this.edges.get(src).remove(tgt);
			}
		}
		for (String src : remove)
			this.edges.remove(src);
		this.setProperties();
		System.out.println(this.numEdges  + "\t" + this.numNodes);
	}
	
	/**
	 * Keeps edges such that the resulting network consists only of
	 * SIM, multi-output FFL and DOR. This is achieved by a triple loop
	 * over the edges:
	 * <ol>
	 * 	<li>Keep for each gene only its best regulator (SIM)</li>
	 * 	<li>Keep for each gene its second best regulator if it is itself
	 * 	regulated by the first regulator (FFL)</li>
	 * 	<li>Keep additional regulators if they already regulate other genes together with 
	 * 	the existing regulators (DOR)</li>
	 * </ol>
	 * 
	 * @param numFFL number of times to repeat step 2
	 * @param numDOR number of times to repeat step 3
	 * 
	 * @author tomic
	 */
	public void keepOnlyNetworkMotifs (int numFFL, int numDOR){
		System.out.println(this.numEdges);
		// sort edges
		this.sortEdges();
		// rebuild network using inverse edges
		this.inverseEdges = new HashMap<String, HashMap<String, Double>>();
		// SIM loop: assign best regulator to each node
		System.out.println("Creating SIM's ...");
		for (Edge e : this.edgeList){
			if (!this.inverseEdges.containsKey(e.target)){
				this.inverseEdges.put(e.target, new HashMap<String, Double>());
				this.inverseEdges.get(e.target).put(e.source, e.weight);
			}
		}
		System.out.println(this.inverseEdges.keySet().size());
		// FFL loop
		for (int k=0; k<numFFL; k++){
			System.out.println("Creating FFL's ...");
			for (Edge e : this.edgeList){
				// check if regulator already has a parent
				if (this.inverseEdges.containsKey(e.source)){
					// get current regulator of source and target
					Set<String> regS = this.inverseEdges.get(e.source).keySet();
					Set<String> regT = this.inverseEdges.get(e.target).keySet();
					// compute intersection
					ArrayList<String> regST = SetOperation.intersection(SetOperation.set2list(regS), SetOperation.set2list(regT));
					// add edge if ffl
					if (!regST.isEmpty()){
						this.inverseEdges.get(e.target).put(e.source, e.weight);
					}
				}
			}
		}
		// DOR loop
		for (int k=0; k<numDOR; k++){
			System.out.println("Creating DOR's ...");
			for (Edge e : this.edgeList){
				// get current regulators of target
				Set<String> regT = this.inverseEdges.get(e.target).keySet();
				// loop over all current parent sets
				for (String tgt : this.inverseEdges.keySet()){
					// regulators of tgt
					Set<String> regTt = this.inverseEdges.get(tgt).keySet();
					// check if regTt contains source
					if (regTt.contains(e.source)){
						// check if regTt and regT intersect
						if (!SetOperation.intersection(SetOperation.set2list(regT), SetOperation.set2list(regTt)).isEmpty()){
							this.inverseEdges.get(e.target).put(e.source, e.weight);
							break;
						}
					}
				}
			}
		}
		this.setEdgesFromInverse();
		this.setProperties();
		System.out.println(this.numEdges);
	}
	
	
	/**
	 * Finds sets of nodes who share the same source nodes.
	 * 
	 * @author tomic
	 */
	public void setCommonSources(){
		this.commonSources = new ArrayList<ArrayList<String>> ();
		// loop over nodes with sources
    	for (String node: this.inverseEdges.keySet()){
    		int ind = 0;
    		boolean foundOne = false;
    		// source nodes of current node
    		Set<String> set1 = this.inverseEdges.get(node).keySet();
    		// loop over current sets of genes with common sources
    		for (int i=0; i<commonSources.size(); i++){
    			// source nodes of current set
    			Set<String> set2 = this.inverseEdges.get(commonSources.get(i).get(0)).keySet();
    			if (set1.equals(set2)){ // check for equality
    				foundOne = true;
    				ind = i;
    				break;
    			}
    		}
    		if (foundOne){
    			// add node to a list in common sources, if one is found
    			this.commonSources.get(ind).add(node);
    		} else { // otherwise create new singleton list in common sources
    			ArrayList<String> single = new ArrayList<String>();
    			single.add(node);
    			this.commonSources.add(single);
    		}
    	}
	}
	
	/**
	 * Returns a list of edges between genes which share at least 1.
	 * regulator. Source and target of edges are in alphabetical order
	 * 
	 * @param rank
	 * @return
	 * 
	 * @author tomic
	 */
	public List<Edge> coRegEdges(){
		List<Edge> coregList = new ArrayList<Edge>();
		for (String src : this.edges.keySet()){
			ArrayList<String> tgts = new ArrayList<String>();
			for (String tgt : this.edges.get(src).keySet())
				tgts.add(tgt);
			Collections.sort(tgts);
			for (int k=0; k<tgts.size()-1; k++){
				for (int l=k+1; l<tgts.size(); l++){
					Edge e = new Edge(tgts.get(k), tgts.get(l));
					coregList.add(e);
				}
			}
		}
		return coregList;
	}
	
	/**
	 * Sets the modules in the network
	 *
	 * @author tomic
	 */
	public void setModules(){
		this.modules = new ArrayList<ArrayList<String>>();
		// loop over signed inverse edges
		for (String node : this.signedInverseEdges.keySet()){
	   		int ind = 0;
    		boolean foundOne = false;
    		// source nodes of current node
    		Set<String> set1 = this.signedInverseEdges.get(node);
    		// loop over current sets of genes with common sources
    		for (int i=0; i<this.modules.size(); i++){
    			// source nodes of current set
    			Set<String> set2 = this.signedInverseEdges.get(this.modules.get(i).get(0));
    			if (set1.equals(set2)){ // check for equality
    				foundOne = true;
    				ind = i;
    				break;
    			}
    		}
    		if (foundOne){
    			// add node to a list in common sources, if one is found
    			this.modules.get(ind).add(node);
    		} else { // otherwise create new singleton list in common sources
    			ArrayList<String> single = new ArrayList<String>();
    			single.add(node);
    			this.modules.add(single);
    		}
		}
	}
	
	/**
	 * Validates the edges by comparison with a given network
	 * 
	 * @param net
	 * 
	 * @author tomic
	 */
	public void validateEdges(Network net){
		for (Edge e : this.edgeList){
			e.validate(net);
		}
		Collections.sort(this.edgeList);
	}

	/**
	 * Sets the modules of the targets of each edge. 
	 * 
	 * @param M
	 * 
	 * @author tomic
	 */
	public void setInModuleEdges(ModuleNetwork M){
		for (Edge e : this.edgeList){
			e.setInModule(M);
		}
		Collections.sort(this.edgeList);
	}
	
	/**
	 * Sets the rank of each edge in another network.
	 * 
	 * @param net
	 * 
	 * @author tomic
	 */
	public void setRankEdges(Network net){
		for (Edge e : this.edgeList){
			e.setRank(net);
		}
		Collections.sort(this.edgeList);
	}
	
	/**
	 * Computes recall - precision for first rank edges
	 * @param rank
	 * @return row index: edge rank, column index: 0=recall, 1=precision, 2=total tp,
	 * 3=total fp according to reference net
	 * 
	 * @author tomic
	 */
	public DenseDoubleMatrix2D computeRecallPrecision(Network refNet, int rank){
		DenseDoubleMatrix2D out = new DenseDoubleMatrix2D(rank, 4);
		Network net = this.getHighestRankEdges(rank);
		net.validateEdges(refNet);
		Collections.sort(net.edgeList);
		double tp=0, fp=0, ne=0, precision, recall;
		for (int k=0; k<rank; k++){
			Edge e = net.edgeList.get(k);
			if (e.valid == 1){
				tp += 1.0;
				ne += 1.0;
			} else if (e.valid == -1){
				fp += 1.0;
				ne += 1.0;
			} 
			if (ne == 0)
				continue;
			recall = tp/refNet.numEdges;
			precision = tp/ne;
			out.set(k, 0, recall);
			out.set(k, 1, precision);
			out.set(k, 2, tp);
			out.set(k, 3, fp);
		}
		return out;
	}

	
	/**
	 * Prints edges with weight higher than some minimum to a file.
	 * 
	 * @param file
	 * @param min
	 * 
	 * @author tomic
	 */
	public void printEdges (String file, double min){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			for (String src : this.edges.keySet())
				for (String tgt : this.edges.get(src).keySet()){
					double w = this.edges.get(src).get(tgt);
					if (w >= min){
						pw.println(src + "\t" + tgt + "\t" + w);
					}
				}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints edges with highest rank to a file
	 * 
	 * @param file
	 * @param min
	 * 
	 * @author tomic
	 */
	public void printEdges (String file, int rank){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			this.sortEdges();
			for (int k=0; k<rank; k++)
				pw.println(this.edgeList.get(k).toString());
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints edges with highest rank to a file
	 * 
	 * @param file
	 * @param rank
	 * @param sign +1 or -1
	 * 
	 * @author tomic
	 */
	public void printEdges (String file, int rank, int sign){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			this.sortEdges();
			for (int k=0; k<rank; k++){
				Edge e = this.edgeList.get(k);
				if (e.sign == sign)
					pw.println(e.toString());
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints edges with highest rank to a file, sorted by regulator and indicates if they belong
	 * to a reference net or not
	 * 
	 * @param file
	 * @param min
	 * 
	 * @author tomic
	 */
	public void printEdges (String file, int rank, Network refNet){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			Network net = this.getHighestRankEdges(rank);
			net.validateEdges(net);
			for (String src : net.edges.keySet()){
				ArrayList<Edge> elst = new ArrayList<Edge>();
				for (String tgt : net.edges.get(src).keySet()){
					Edge e = new Edge(src, tgt, net.edges.get(src).get(tgt));
					e.validate(refNet);
					elst.add(e);
				}
				Collections.sort(elst);
				for (Edge e : elst)
					pw.println(e.toString());
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints edges with highest rank to a file, sorted by regulator, indicates if they belong
	 * to a reference net or not, and prints to which modules in a module network the targets 
	 * belong.
	 * 
	 * @param file
	 * @param min
	 * 
	 * @author tomic
	 */
	public void printEdges (String file, int rank, Network refNet, ModuleNetwork M){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			Collections.sort(this.edgeList);
			for (int k=0; k<rank; k++){
				Edge e = this.edgeList.get(k);
				e.validate(refNet);
				e.setInModule(M);
				pw.println(e.toString());
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints edges with highest rank to a file, sorted by regulator, 
	 * and prints to which modules in a module network the targets 
	 * belong.
	 * 
	 * @param file
	 * @param min
	 * 
	 * @author tomic
	 */
	public void printEdges (String file, int rank, ModuleNetwork M){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			Collections.sort(this.edgeList);
			for (int k=0; k<rank; k++){
				Edge e = this.edgeList.get(k);
				e.setInModule(M);
				pw.println(e.toString());
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints the edges that are common between this network and another one
	 * at each rank level.
	 * 
	 * @param file
	 * @param net
	 * @param rank
	 */
	public void printCommonEdgesRank(String file, Network net, int rank){
		this.sortEdges();
		net.sortEdges();
		Network net1 = this.getHighestRankEdges(1);
		Network net2 = net.getHighestRankEdges(1);
		int com = 0;
		if (net.containsEdge(this.edgeList.get(0)))
			com = 1;
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			for (int k=1; k<rank; k++){
				Edge e1 = this.edgeList.get(k);
				Edge e2 = net.edgeList.get(k);
				if (e1.equals(e2))
					com += 1;
				else if (net1.containsEdge(e2))
					com += 1;
				net1.edgeList.add(e1);
				net2.edgeList.add(e2);
				pw.println(k+1 + "\t" + com);
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Prints the modules in the network with their regulators.
	 * 
	 * @param file
	 * 
	 * @author tomic
	 */
	public void printModules (String file){
		this.setModules();
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			for (ArrayList<String> set : this.modules){
				String s = set.get(0);
				for (String reg : this.signedInverseEdges.get(set.get(0)))
					s += "\t" + reg;
				pw.println(s);
				for (int i=1; i<set.size(); i++){
					pw.println(set.get(i));
				}
				pw.println();
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
	/**
	 * Creates a new network where each target of a source is an indirect target
	 * in the original network and the weight is the path length. 
	 * 
	 * @param maxLength 
	 * @return
	 * 
	 * @author tomic
	 */
	public Network indirectEdges (int maxLength){
		Network net = new Network();
		
		for (String src : this.edges.keySet()){
			net.edges.put(src, new HashMap<String, Double>());
			for (String tgt : this.edges.get(src).keySet())
				net.edges.get(src).put(tgt, 1.0);
			Set<String> reachableNodes = this.edges.get(src).keySet();
			// loop until we find target or can't continue
			int length = 1;
			while (!reachableNodes.isEmpty() && length<maxLength) {
				Set<String> nextReachable = new HashSet<String>();
				length += 1;
				for (String rgene: reachableNodes) // current genes
					if (this.edges.keySet().contains(rgene))
						for (String tgene: this.edges.get(rgene).keySet()) { // children of current genes
							nextReachable.add(tgene);
						}
				reachableNodes = nextReachable;
				for (String tgt : reachableNodes){
					net.edges.get(src).put(tgt, (double)length);
				}
			}
		}
		
		net.setProperties();
		
		return net;
	}
	
	/**
	 * Tries to find a directed path between 2 nodes.
	 * 
	 * @param source source node.
	 * @param target target node
	 * @return length of path of one is found, -1 otherwise, 0 if one of the nodes does not belong to the network
	 * 
	 * @author tomic
	 */
	public int existsPath(String source, String target) {
		if (this.nodeList.contains(source) && this.nodeList.contains(target)){
			int length = 1;
			// initially start with children of source
			Set<String> reachableNodes = this.edges.get(source).keySet();
			// loop until we find target or can't continue
			while (!reachableNodes.isEmpty() && !reachableNodes.contains(target)) {
				Set<String> nextReachable = new HashSet<String>();
				length += 1;
				for (String rgene: reachableNodes) // current genes
					if (this.edges.keySet().contains(rgene))
						for (String tgene: this.edges.get(rgene).keySet()) { // children of current genes
							nextReachable.add(tgene);
						}
				reachableNodes = nextReachable;
			}
			if (reachableNodes.isEmpty()){
				return -1;
			}
			else if (reachableNodes.contains(target)){
				return length;
			}
			else {// java not smart enough to see that this wouldn't happen
				return -1;
			}
		} else 
			return 0;
	
	}


}
