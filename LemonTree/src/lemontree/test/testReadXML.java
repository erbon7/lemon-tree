package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;

public class testReadXML {

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
		m.setDataMeanAndSDFromModuleset();

		m.readRegTreeXML("/Users/eric/reg.xml.gz");
		m.setTestSplits();
		m.printRegulators("/Users/eric/test.txt", false, false);
		m.printRandomRegulators("/Users/eric/rand.txt", false);
		
		m.setTopRegulators("/Users/eric/test.txt");
		m.setModuleMeanSigma();
		m.checkExperiments();

	}
}
