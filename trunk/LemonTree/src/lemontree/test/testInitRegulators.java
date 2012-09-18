package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;

public class testInitRegulators {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String datasetFile ="/Users/eric/work/lemone_test/data/all_prostate.txt";
		String regulatorsFile="/Users/eric/work/lemone_test/data/reg_mir_list.txt";
		String clusterFile="/Users/eric/work/lemone_test/data/tc4.txt";

		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile, null);
		m.readClusters(clusterFile);
		m.readRegulators(regulatorsFile);
		m.initStatisticsAndScore();
		m.setDataMeanAndSDFromModuleset();

		int num_runs=1;
		int num_steps=2;//1100;
		int burn_in=1;//100;
		int sample_steps=2;//100;
		double score_gain = 0.0;
		boolean use_bayesian_score=true;
		m.gibbsSamplerExpts(num_runs, num_steps, burn_in, sample_steps, score_gain, use_bayesian_score);
		double beta_reg=20;
		int num_reg=10;
		m.assignRegulatorsNoAcyclStoch(beta_reg, num_reg);
		m.printRegulators("/Users/eric/reg.txt", false, false);
		m.printRandomRegulators("/Users/eric/rand.txt", false);
		m.writeRegTreeXML("/Users/eric/reg.xml.gz");
	}

}
