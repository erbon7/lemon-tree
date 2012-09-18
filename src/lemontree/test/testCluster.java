package lemontree.test;

import lemontree.modulenetwork.ModuleNetwork;

public class testCluster {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String fileName ="/Users/eric/work/lemone_test/tc4.txt";
		ModuleNetwork m = new ModuleNetwork();
		m.readExpressionMatrix("/Users/eric/work/lemone_test/all_NF.txt",null);
		m.readClusters(fileName);
		
	}

}
