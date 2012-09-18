package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;

public class TestRangeReg {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String datasetFile ="/Users/eric/work/lemone_test/all_NF.txt";
		String regulatorsFile="/Users/eric/work/lemone_test/reg_mir_list.txt";
		String clusterFile="/Users/eric/work/lemone_test/tc4.txt";

		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile,null);
		m.readClusters(clusterFile);
		m.readRegulators(regulatorsFile);
		m.initStatisticsAndScore();

		// read xml and assign range regulators
		/*m.readRegTreeXML("/Users/eric/reg.xml.gz");
		
		double beta_reg=20;
		int num_reg=10;
		m.assignRegulatorsNoAcyclStoch(beta_reg, num_reg, 0, 0);
		m.printRegulators("/Users/eric/reg.txt", true, false);
		m.writeRegTreeXML("/Users/eric/reg.xml.gz");
		*/
		
		m.readRegTreeXML("/Users/eric/reg.xml.gz");
		m.setTestSplits();
		m.printRegulators("/Users/eric/reg.txt", true, false);
		
	}

}
