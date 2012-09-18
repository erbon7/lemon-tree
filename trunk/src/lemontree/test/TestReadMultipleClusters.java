package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;

public class TestReadMultipleClusters {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String datasetFile ="/Users/eric/work/lemone_test/all_NF.txt";
		//String geneFile ="/Users/eric/work/lemone_test/gene_list";
		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile, null);
		m.readMultipleClusters("/Users/eric/work/lemone_test/cluster_list");
		System.out.println(m.moduleSets.size());
		System.out.println(m.moduleSet.size());
	}

}
