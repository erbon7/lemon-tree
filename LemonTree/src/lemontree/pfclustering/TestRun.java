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


public class TestRun {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		// Read network
		String dir = "/Users/tom/Eclipse/EdgeClustering/test";
		//String outFile = "STAGE_VF_centroid.txt";
		String file = "ppi_yu_trimax.txt";
		Network net = new Network(dir, file, true);
		PFClustering clust = new PFClustering(net,1e-5,500,4,2.0);
		clust.edgeClustering();
		//clust.writeClustersX(dir, outFile);
	}

}
