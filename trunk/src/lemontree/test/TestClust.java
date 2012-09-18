package lemontree.test;

import java.io.File;
import java.io.PrintWriter;

import lemontree.pfclustering.Edge;

import lemontree.modulenetwork.CentroidClustering;
import lemontree.modulenetwork.Gene;
import lemontree.modulenetwork.ModuleNetwork;

public class TestClust {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
//		String matrix_file = "/Users/eric/work/sysbiol/n6/agilent_mRNA/dataCS";
//		String cluster_file_list = "/Users/eric/files";
//		String output_file = "/Users/eric/out.txt";
		
		String matrix_file = args[0];
		String cluster_file_list = args[1];
		String output_file = args[2];
		String node_clustering = args[3];  
		
		ModuleNetwork m = new ModuleNetwork();
		m.readExpressionMatrix(matrix_file, null);
		m.readMultipleClusters(cluster_file_list);
		boolean b = true;
		if (node_clustering.equalsIgnoreCase("false"))
			b = false;
		
		double minWeight = 0.33;
		int minClustSize = 10;
		int minClustScore = 2;
		
		CentroidClustering cc = new CentroidClustering(m,b,minWeight,minClustSize,minClustScore);
		
//		try {
//			PrintWriter out = new PrintWriter(new File("pairwise_prob_025_java"));
//			for (int i=0;i< cc.coClustNetwork.edges.size(); i++) {
//				Edge e = cc.coClustNetwork.edges.get(i);
//				out.println(e.from+" "+e.to+" "+e.weight);
//			}
//			out.close();
//		}
//		catch (Exception e) {
//			e.printStackTrace();
//		}
//		System.exit(1);
		
		cc.doCentroidClustering();
		
		try {
			PrintWriter out = new PrintWriter(new File(output_file));
			for (int i=0;i<cc.centroidClusters.size();i++) {
				for (Gene g : cc.centroidClusters.get(i)) {
					out.println(g.name+"\t"+i);
				}
			}
			out.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}

	}

}
