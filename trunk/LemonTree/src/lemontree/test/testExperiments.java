package lemontree.test;

import lemontree.modulenetwork.Module;
import lemontree.modulenetwork.ModuleNetwork;

public class testExperiments {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String prefix = "/bioinfo/users/ebonnet/TreeLemon/gbm2/";
		String datasetFile = prefix+"reg_cn/data";
		String regulatorsFile= prefix+"reg_cn/del_list";
		String clusterFile=prefix+"tight_clusters/tc_10_033.txt";

		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile, null);
		m.readClusters(clusterFile);
		m.readRegulators(regulatorsFile);
		m.initStatisticsAndScore();
		//m.setDataMeanAndSDFromModuleset();
		
		/*
		 * sample experiments and save xml file
		 */
//		int num_runs=1;
//		int num_steps=1100;//1100;
//		int burn_in=100;//100;
//		int sample_steps=100;//100;
//		double score_gain = 0.0;
//		boolean use_bayesian_score=true;
//		m.gibbsSamplerExpts(num_runs, num_steps, burn_in, sample_steps, score_gain, use_bayesian_score);
//		m.writeRegTreeXML("/Users/eric/reg.xml.gz");

		/*
		 * Read from xml and assign regulators
		 */
		
		m.readRegTreeXML(prefix+"experiments/exp.xml.gz");
		
		for (Module mod : m.moduleSet) {
			mod.hierarchicalTree.testLevel(4);
		}
		
		m.writeRegTreeXML(prefix+"experiments/exp_lev4.xml.gz");
		
//		double beta_reg=20;
//		int num_reg=10;
//		m.assignRegulatorsNoAcyclStoch(beta_reg, num_reg);
//		m.printRegulators("/Users/eric/reg.txt", true, false);
	
	}

}
