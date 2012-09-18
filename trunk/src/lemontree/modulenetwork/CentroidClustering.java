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
import java.io.PrintWriter;
import java.util.*;

import de.unifreiburg.frias.omics.pfclustering.*;

public class CentroidClustering {

	public ModuleNetwork moduleNetwork;

	public ArrayList<HashSet<Gene>> centroidClusters;

	public Network coClustNetwork;
	
	public boolean nodeClustering; // node (TRUE) or edge (FALSE) clustering
	
	public int minClustSize; // minimal number of nodes in a cluster
	
	public double minClustScore; // minimal score of a cluster
	
	public double tolerance; // convergence tolerance in clustering algorithm
	
	public int maxstep; // max number of iterations in perronvector calculations
	

	public CentroidClustering(){

	}

	public CentroidClustering (ModuleNetwork moduleNetwork){
		this.moduleNetwork = moduleNetwork;
		this.setCoClusteringNetwork(0.25);
		this.nodeClustering = true;
		this.minClustSize = 4;
		this.minClustScore = 2.0;
		this.tolerance = 1E-5;
		this.maxstep = 1000;
	}

	public CentroidClustering (ModuleNetwork moduleNetwork, boolean nodeClustering, double minWeight){
		this.moduleNetwork = moduleNetwork;
		this.setCoClusteringNetwork(minWeight);
		this.nodeClustering = nodeClustering;
		this.minClustSize = 4;
		this.minClustScore = 2.0;
		this.tolerance = 1E-5;
		this.maxstep = 1000;
	}
	
	public CentroidClustering (ModuleNetwork moduleNetwork, boolean nodeClustering, double minWeight, int minClustSize, double minClustScore){
		this.moduleNetwork = moduleNetwork;
		this.setCoClusteringNetwork(minWeight);
		this.nodeClustering = nodeClustering;
		this.minClustSize = minClustSize;
		this.minClustScore = minClustScore;
		this.tolerance = 1E-5;
		this.maxstep = 1000;
	}

	public CentroidClustering (ModuleNetwork moduleNetwork, boolean nodeClustering, double minWeight, int minClustSize, double minClustScore, double tolerance, int maxstep){
		this.moduleNetwork = moduleNetwork;
		this.setCoClusteringNetwork(minWeight);
		this.nodeClustering = nodeClustering;
		this.minClustSize = minClustSize;
		this.minClustScore = minClustScore;
		this.tolerance = tolerance;
		this.maxstep = maxstep;
	}
	
	public void doCentroidClustering (){
		PFClustering clust = new PFClustering(this.coClustNetwork, this.tolerance, this.maxstep, this.minClustSize, this.minClustScore);
		// node- or edge-clustering
		if (this.nodeClustering)
			clust.nodeClustering();
		else
			clust.edgeClustering();
		// Convert string-based clusters to gene clusters
		this.centroidClusters = new ArrayList<HashSet<Gene>>();
		for (HashSet<String> cluster : clust.clustersX){
			HashSet<Gene> geneCluster = new HashSet<Gene>();
			for (String genename : cluster){
				Gene gene = this.moduleNetwork.geneMap.get(genename);
				geneCluster.add(gene);
			}
			this.centroidClusters.add(geneCluster);
		}
	}

	/**
	 * Set the coclustering network: an undirected network between genes with weight equal 
	 * to their coclustering frequency, keeping only edges with a minimal weight
	 */
	public void setCoClusteringNetwork(double minWeight){
		System.out.print("Init CoClustering...");
		this.coClustNetwork = new Network(true);
		if (this.moduleNetwork.moduleSets.isEmpty())
			this.moduleNetwork.moduleSets.add(this.moduleNetwork.moduleSet);
		double numSets = this.moduleNetwork.moduleSets.size();

		this.moduleNetwork.setInModules();
		ArrayList<Gene> genes = this.moduleNetwork.geneSet;
		// loop over gene pairs
		for (int k=0; k<genes.size()-1; k++)
			for (int l=k+1; l<genes.size(); l++){
				double nr = 0.0;
				for (int m=0; m<numSets; m++) {
					if (this.moduleNetwork.inModules.get(genes.get(k).name).get(m).equals(this.moduleNetwork.inModules.get(genes.get(l).name).get(m)))
						nr += 1.0;
				}
				double w = nr/numSets;
				if (w > minWeight){
					Edge e = new Edge(genes.get(k).name, genes.get(l).name, w , true);
					this.coClustNetwork.addEdge(e,false);
				}
			}

		// set network properties
		this.coClustNetwork.setProperties();
		System.out.println(" [ok]");
	}
	
	/**
	 * Print clusters to file.
	 * 
	 * @param filename
	 */
	public void printClusters (String filename) {
		System.out.print("Printing clusters...");
		
		try {
			PrintWriter out = new PrintWriter(new File(filename));
			for (int i=0;i<centroidClusters.size();i++) {
				for (Gene g : centroidClusters.get(i)) {
					out.println(g.name+"\t"+i);
				}
			}
			out.close();
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		System.out.println(" [ok]");
	}

}
