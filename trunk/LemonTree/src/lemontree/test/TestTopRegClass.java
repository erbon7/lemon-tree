package lemontree.test;

import lemontree.utils.DrawModules;
import lemontree.modulenetwork.ModuleNetwork;

public class TestTopRegClass {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		
		String datasetFile ="/Users/eric/work/lemone_test/all_NF.txt";
		String regulatorsFile="/Users/eric/work/lemone_test/reg_mir_list.txt";
		String clusterFile="/Users/eric/work/lemone_test/tc4.txt";

		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile,null);
		m.readClusters(clusterFile);
		m.readRegulators(regulatorsFile);
		m.initStatisticsAndScore();

		m.readRegTreeXML("/Users/eric/reg.xml.gz");
		m.setTestSplits();
		m.setTopRegulators("/Users/eric/reg.txt");
		m.setModuleMeanSigma();
		m.checkExperiments();
		
		DrawModules dm = new DrawModules(m);
		dm.drawAllModules();

	}

}
