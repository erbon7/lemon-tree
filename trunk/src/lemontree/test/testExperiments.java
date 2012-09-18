package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;

public class testExperiments {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String datasetFile ="/Users/eric/work/lemone_test/all_NF.txt";
		String regulatorsFile="/Users/eric/work/lemone_test/reg_mir_list.txt";
		String clusterFile="/Users/eric/work/lemone_test/tc4.txt";

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
		
		m.readRegTreeXML("/Users/eric/reg.xml.gz");
		double beta_reg=20;
		int num_reg=10;
		m.assignRegulatorsNoAcyclStoch(beta_reg, num_reg);
		m.printRegulators("/Users/eric/reg.txt", true, false);
	}

}
