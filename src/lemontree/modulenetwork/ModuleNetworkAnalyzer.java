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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import cern.colt.matrix.impl.*;
import lemontree.networktools.*;

/**
 * This class collects several methods for analyzing module networks.
 * 
 * @author tomic
 *
 */

public class ModuleNetworkAnalyzer {

	/**
	 * Module network to be analyzed.
	 */
	public ModuleNetwork moduleNetwork;
	
	/**
	 * Gene network obtained from module network. Connects each regulator
	 * of a module to all the genes in that module. For multiple regulation trees 
	 * per module and/or multiple module sets, the gene network edges are weighted
	 * by {@link Module#regulatorWeightsMax}.
	 */
	public Network geneNetworkMax;
	
	/**
	 * Gene network obtained from module network. Connects each regulator
	 * of a module to all the genes in that module. For multiple regulation trees 
	 * per module and/or multiple module sets, the gene network edges are weighted
	 * by {@link Module#regulatorWeightsSum}.
	 */
	public Network geneNetworkSum;
	
	/**
	 * Random gene network obtained from module network. Connects the randomly
	 * assigned regulators to the genes in the module. For multiple regulation trees 
	 * per module and/or multiple module sets, the network edges are weighted
	 * by {@link Module#regulatorWeightsRandom}.
	 */
	public Network randomGeneNetwork;
	
	/**
	 * Regulator network obtained from {@link ModuleNetworkAnalyzer#geneNetworkMax}. 
	 * Subnetwork of gene network which only connects the regulators (= nodes with
	 * targets). 
	 */
	public Network regulatorNetworkMax;
	
	/**
	 * Gene network linking genes that belong to the same module. For multiple
	 * module sets the edges are weighted by the number of occurrences.
	 */
	public Network coexpressionNetwork;
	
	/**
	 * Reference network with known interactions.
	 */
	public Network referenceNetwork;
	
	/**
	 * Subnetwork of {@link ModuleNetworkAnalyzer#geneNetworkMax} containing only the 
	 * edges between genes of {@link ModuleNetworkAnalyzer#referenceNetwork}.
	 */
	public Network geneInRefNetworkMax;
	
	/**
	 * Empty constructor
	 */
	public ModuleNetworkAnalyzer() {
		
	}
	
	/**
	 * Construct analyzer from module network
	 * 
	 * @param M ModuleNetwork
	 * @param weighted use weighted edges or not
	 */
	public ModuleNetworkAnalyzer(ModuleNetwork M, boolean weighted) {
		this.moduleNetwork = M;
		if (weighted){
			this.setGeneNetworkMax();
			this.setGeneNetworkSum();
//			this.setRandomGeneNetwork();
			// scale edge weights
			double max = this.geneNetworkMax.maxWeight();
			this.geneNetworkMax.scaleEdgeWeights(max);
			this.geneNetworkMax.setProperties();
			double max2 = this.geneNetworkSum.maxWeight();
			this.geneNetworkSum.scaleEdgeWeights(max2);
			this.geneNetworkMax.setProperties();
//			this.randomGeneNetwork.scaleEdgeWeights(max2);
//			this.randomGeneNetwork.setProperties();
		} else {
			this.setGeneNetworkSimple();
		}
//		this.setRegulatorNetwork();
		this.referenceNetwork = new Network();
		this.geneInRefNetworkMax = new Network();
	}
	
	/**
	 * Construct analyzer from module network and reference network.
	 * 
	 * @param M ModuleNetwork
	 * @param weighted use weighted edges or not
	 * @param dir directory where reference network is located.
	 * @param refNetFile name of reference network file.
	 */
	public ModuleNetworkAnalyzer(ModuleNetwork M, boolean weighted, String dir, 
			String refNetFile){
		this.moduleNetwork = M;
		if (weighted){
			this.setGeneNetworkMax();
			this.setGeneNetworkSum();
//			this.setRandomGeneNetwork();
			// scale edge weights
			double max = this.geneNetworkMax.maxWeight();
			this.geneNetworkMax.scaleEdgeWeights(max);
			this.geneNetworkMax.setProperties();
			double max2 = this.geneNetworkSum.maxWeight();
			this.geneNetworkSum.scaleEdgeWeights(max2);
			this.geneNetworkMax.setProperties();
//			this.randomGeneNetwork.scaleEdgeWeights(max2);
//			this.randomGeneNetwork.setProperties();
		} else {
			this.setGeneNetworkSimple();
		}
//		this.setRegulatorNetwork();
		this.referenceNetwork = new Network(dir, refNetFile);
		this.setGeneInRefNetwork();
	}
	
	
	/**
	 * Construct analyzer from module network (for data, genes, etc.) and 
	 * gene network file (containing matrix of edge weights). This used to 
	 * analyze the output of other programs (e.g., CLR).
	 * 
	 * @param networkFile
	 */
	public ModuleNetworkAnalyzer(ModuleNetwork M, String dir, 
			String refNetFile, String networkFile){
		this.moduleNetwork = M;
		this.readGeneNetwork(networkFile);
		this.setRegulatorNetwork();
		this.referenceNetwork = new Network(dir, refNetFile);
		this.setGeneInRefNetwork();
	}
	
	
	/**
	 * Creates the gene network, no weighted edges. Assumes 1 set of modules and
	 * 1 tree per module.
	 * 
	 * @author tomic
	 */
	public void setGeneNetworkSimple(){
		this.geneNetworkMax = new Network();
		for (Module mod : this.moduleNetwork.moduleSet){
			for (Gene reg : mod.parents){
				if (!this.geneNetworkMax.edges.containsKey(reg.name))
					this.geneNetworkMax.edges.put(reg.name, new HashMap<String, Double>());
				for (Gene gene : mod.genes)
					this.geneNetworkMax.edges.get(reg.name).put(gene.name, 1.0);
			}
		}
		this.geneNetworkMax.setProperties();
	}

	/**
	 * Sets the gene network with weighted edges using {@link Module#regulatorWeightsMax}.
	 * 
	 * @author tomic
	 */
	public void setGeneNetworkMax() {
		this.geneNetworkMax = new Network();
		// check for multiple module sets
		if (this.moduleNetwork.moduleSets.isEmpty())
			this.moduleNetwork.moduleSets.add(this.moduleNetwork.moduleSet);
		for (ArrayList<Module> modset : this.moduleNetwork.moduleSets){
			for (Module mod : modset){
				for (Gene reg : mod.regulatorWeightsMax.keySet()){
					// add reg if it's not already there
					if (!this.geneNetworkMax.edges.containsKey(reg.name))
						this.geneNetworkMax.edges.put(reg.name, new HashMap<String,Double>());
					for (Gene gene : mod.genes){
						// if the edge is not yet in the network, add it with reg weight for this module
						if (!this.geneNetworkMax.edges.get(reg.name).containsKey(gene.name)){
							this.geneNetworkMax.edges.get(reg.name).put(gene.name, mod.regulatorWeightsMax.get(reg));
						} else {
							// if the edge is already there, add the reg weight for this module
							double w = this.geneNetworkMax.edges.get(reg.name).get(gene.name);
							this.geneNetworkMax.edges.get(reg.name).put(gene.name, w+mod.regulatorWeightsMax.get(reg));
						}
					}
				}
			}
		} 
   		this.geneNetworkMax.setProperties();
	}
	
	/**
	 * Sets the gene network with weighted edges using {@link Module#regulatorWeightsSum}.
	 * 
	 * @author tomic
	 */
	public void setGeneNetworkSum() {
		this.geneNetworkSum = new Network();
		// check for multiple module sets
		if (this.moduleNetwork.moduleSets.isEmpty())
			this.moduleNetwork.moduleSets.add(this.moduleNetwork.moduleSet);
		for (ArrayList<Module> modset : this.moduleNetwork.moduleSets){
			for (Module mod : modset){
				for (Gene reg : mod.regulatorWeightsSum.keySet()){
					// add reg if it's not already there
					if (!this.geneNetworkSum.edges.containsKey(reg.name))
						this.geneNetworkSum.edges.put(reg.name, new HashMap<String,Double>());
					for (Gene gene : mod.genes){
						// if the edge is not yet in the network, add it with reg weight for this module
						if (!this.geneNetworkSum.edges.get(reg.name).containsKey(gene.name)){
							this.geneNetworkSum.edges.get(reg.name).put(gene.name, mod.regulatorWeightsSum.get(reg));
						} else {
							// if the edge is already there, add the reg weight for this module
							double w = this.geneNetworkSum.edges.get(reg.name).get(gene.name);
							this.geneNetworkSum.edges.get(reg.name).put(gene.name, w+mod.regulatorWeightsSum.get(reg));
						}
					}
				}
			}
		} 
   		this.geneNetworkSum.setProperties();
	}
	

	/**
	 * Reads the gene network from a file.
	 * @param dir
	 * @param networkFile
	 * 
	 * @author tomic
	 */
	public void readGeneNetwork (String networkFile){
		this.geneNetworkMax = new Network();
		try {
			Scanner fileScanner = new Scanner(new File(networkFile)).useDelimiter("\\n");
			// walk through file
			int k=0;
			while (fileScanner.hasNext()){
				Scanner line = new Scanner(fileScanner.next().toString()).useDelimiter("\\s");
				// we can process all values by only processing regulator lines
				Gene reg = this.moduleNetwork.geneSet.get(k);
				if (this.moduleNetwork.regulatorSet.contains(reg)){
					this.geneNetworkMax.edges.put(reg.name, new HashMap<String, Double>());
					int l=0;
					while (line.hasNext()){
						Gene gene = this.moduleNetwork.geneSet.get(l);
						double w = line.nextDouble();
						// add edge between gene(k) and gene(l) with weight w
						if (w>0.0)
							this.geneNetworkMax.edges.get(reg.name).put(gene.name, w);
						l++;
					}
				}
				k++;
			}
		} catch (FileNotFoundException e) {
			System.out.println(e);
		}
		this.geneNetworkMax.setProperties();
	}

	
	/**
	 * Truncates edges in geneNetwork with weight below certain value.
	 * 
	 * @param min
	 * 
	 * @author tomic
	 */
	public void truncGeneNetwork (double min){
		this.geneNetworkMax.removeEdges(min);
		if (!this.geneInRefNetworkMax.edges.isEmpty())
			this.geneInRefNetworkMax.removeEdges(min);
		if (!this.regulatorNetworkMax.edges.isEmpty())
			this.regulatorNetworkMax.removeEdges(min);
	}

	/**
	 * Sets the random gene network with weighted edges using {@link Module#regulatorWeightsRandom},
	 * keeping only edges with minimal weight.
	 * 
	 * @author tomic
	 */
	public void setRandomGeneNetwork () {
		this.randomGeneNetwork = new Network();
		// check for multiple module sets
		if (this.moduleNetwork.moduleSets.isEmpty())
			this.moduleNetwork.moduleSets.add(this.moduleNetwork.moduleSet);
		for (int k=0; k<this.moduleNetwork.moduleSets.size(); k++){
			ArrayList<Module> modset = this.moduleNetwork.moduleSets.get(k);
			for (Module mod : modset){
				for (Gene reg : mod.regulatorWeightsRandom.keySet()){
					// add reg if it's not already there
					if (!this.randomGeneNetwork.edges.containsKey(reg.name))
						this.randomGeneNetwork.edges.put(reg.name, new HashMap<String,Double>());
					for (Gene gene : mod.genes){
						// if the edge is not yet in the network, add it with reg weight for this module
						if (!this.randomGeneNetwork.edges.get(reg.name).containsKey(gene.name)){
							this.randomGeneNetwork.edges.get(reg.name).put(gene.name, mod.regulatorWeightsRandom.get(reg));
						} else {
							// if the edge is already there, add the reg weight for this module
							double w = this.randomGeneNetwork.edges.get(reg.name).get(gene.name);
							this.randomGeneNetwork.edges.get(reg.name).put(gene.name, w+mod.regulatorWeightsRandom.get(reg));
						}
					}
				}
			}
		} 
		this.randomGeneNetwork.setProperties();
	}
	
	
	/**
	 * Creates the regulator network from the gene network.
	 * 
	 * @author tomic
	 */
	public void setRegulatorNetwork(){
		this.regulatorNetworkMax = new Network();
		for (String reg : this.geneNetworkMax.edges.keySet()){
			this.regulatorNetworkMax.edges.put(reg, new HashMap<String, Double>());
			for (String gene : this.geneNetworkMax.edges.get(reg).keySet()){
				if (this.geneNetworkMax.edges.containsKey(gene)){
					double w = this.geneNetworkMax.edges.get(reg).get(gene);
					this.regulatorNetworkMax.edges.get(reg).put(gene, w);
				}
			}
		}
		// some regulators may have zero regulator targets
		HashSet<String> removeRegs = new HashSet<String>();
		for (String reg : this.regulatorNetworkMax.edges.keySet())
			if (this.regulatorNetworkMax.edges.get(reg).isEmpty())
				removeRegs.add(reg);
		for (String reg : removeRegs)
			this.regulatorNetworkMax.edges.remove(reg);
		this.regulatorNetworkMax.setProperties();
	}

	/**
	 * Sets the coexpression network. This is an undirected network, by convention
	 * we create the lower triangular half of the connectivity matrix
	 * 
	 * @author tomic
	 */
	public void setCoexpressionNetwork(){
		this.coexpressionNetwork = new Network();
		// check for multiple module sets
		if (this.moduleNetwork.moduleSets.isEmpty())
			this.moduleNetwork.moduleSets.add(this.moduleNetwork.moduleSet);
		for (ArrayList<Module> modset : this.moduleNetwork.moduleSets){
			for (Module mod : modset){
				for (Gene gene1 : mod.genes){
					for (Gene gene2 : mod.genes){
						String source, target;
						if (gene1.number >= gene2.number){
							source = gene1.name;
							target = gene2.name;
						}
						else {
							source = gene2.name;
							target = gene1.name;
						}
						if (!this.coexpressionNetwork.edges.containsKey(source)){
							this.coexpressionNetwork.edges.put(source, new HashMap<String, Double>());
						}
						if (!this.coexpressionNetwork.edges.get(source).containsKey(target)){
							this.coexpressionNetwork.edges.get(source).put(target, 1.0);
						} else {
							double num = this.coexpressionNetwork.edges.get(source).get(target);
							this.coexpressionNetwork.edges.get(source).put(target, num+1);
						}
					}
				}
			}
		} 
   		this.coexpressionNetwork.setProperties();
	}
	
	/**
	 * Creates the {@link ModuleNetworkAnalyzer#geneInRefNetworkMax} from the 
	 * {@link ModuleNetworkAnalyzer#geneNetworkMax} and the 
	 * {@link ModuleNetworkAnalyzer#referenceNetwork}
	 *
	 * @author tomic
	 */
	public void setGeneInRefNetwork(){
		this.geneInRefNetworkMax = new Network();
		for (String reg : this.geneNetworkMax.edges.keySet()){
			if (this.referenceNetwork.nodeList.contains(reg)){
				this.geneInRefNetworkMax.edges.put(reg, new HashMap<String, Double>());
				for (String gene : this.geneNetworkMax.edges.get(reg).keySet())
					if (this.referenceNetwork.nodeList.contains(gene)){
						double w = this.geneNetworkMax.edges.get(reg).get(gene);
						this.geneInRefNetworkMax.edges.get(reg).put(gene, w);
					}
			}
		}
		// some regulators may have zero targets
		HashSet<String> removeRegs = new HashSet<String>();
		for (String reg : this.geneInRefNetworkMax.edges.keySet())
			if (this.geneInRefNetworkMax.edges.get(reg).isEmpty())
				removeRegs.add(reg);
		for (String reg : removeRegs)
			this.geneInRefNetworkMax.edges.remove(reg);
		this.geneInRefNetworkMax.setProperties();
	}
	

	/**
	 * Computes recall - precision for first rank edges
	 * @param rank
	 * @return row index: edge rank, column index: 0=recall, 1=precision, 2=total tp,
	 * 3=total fp according to reference net
	 */
	public DenseDoubleMatrix2D computeRecallPrecision(int rank){
		return this.geneNetworkMax.computeRecallPrecision(this.referenceNetwork, rank);
	}
	
	/**
	 * Prints recall - precision data for different edge weight cutoffs
	 * 
	 * @param file
	 * @param min
	 * 
	 * @author tomic
	 */
	public void printRecallPrecision (String file, int rank){
		try{
			File output = new File(file);
			PrintWriter pw = new PrintWriter(output);
			Network net = this.geneNetworkSum.getHighestRankEdges(rank);
			net.validateEdges(this.referenceNetwork);
			Collections.sort(net.edgeList);
			double tp=0, fp=0, ne=0, precision, recall;
			for (int k=0; k<net.edgeList.size(); k++){
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
				recall = tp/(double)this.referenceNetwork.numEdges;
				precision = tp/ne;
				pw.println(k + "\t" + recall + "\t" + precision + "\t" + (int)tp  + "\t" + (int)fp);
			}
			pw.close();
		} catch (IOException e){
			System.out.println(e);
		}
	}
	
}
