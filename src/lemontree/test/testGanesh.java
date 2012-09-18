package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;
import lemontree.modulenetwork.Gene;
import java.util.*;

public class testGanesh {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String datasetFile ="/Users/eric/work/lemone_test/all_NF.txt";
		String geneFile ="/Users/eric/work/lemone_test/gene_list";
		ModuleNetwork m = new ModuleNetwork();
		//m.initGaneshTask(fileName, 0.1f, 0.0f, 0.1f, 0.1f);
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile, geneFile);
		
		int init_num_clust=0;
		int num_runs=1;
		int burn_in=5;
		int num_steps=10;
		int sample_steps=10;
		double score_gain=0.0;
		boolean use_bayesian_score=true;
		
		m.gibbsSamplerGenes(init_num_clust, num_runs, num_steps, burn_in, sample_steps, score_gain, use_bayesian_score);
		m.writeClusters("/Users/eric/cluster2.txt");
	}

}
