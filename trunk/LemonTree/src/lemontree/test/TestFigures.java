package lemontree.test;

import java.util.ArrayList;

import lemontree.modulenetwork.Module;
import lemontree.modulenetwork.ModuleNetwork;
import lemontree.modulenetwork.TreeNode;
import lemontree.utils.DrawModules;

public class TestFigures {

	public static void main(String[] args) {
		
		String data_prefix = "/Users/eric/gbm_test/data_dir/";
		String data_file = data_prefix + "data";
		String cluster_file = data_prefix + "tc_abs_10_033.txt";
		String tree_file = data_prefix + "exp_abs.xml.gz";
		
		//String reg_file = data_prefix + "amp_reg_list";
		//String top_regulators = data_prefix + "amp_top_reg_list.txt";
		
		String reg_file = data_prefix + "del_reg_list";
		String top_regulators = data_prefix + "del_top_reg_list.txt";
		
		
		boolean use_regulator_mean = true;
		int cut_level = 1;
		
		ModuleNetwork M = new ModuleNetwork();
		//read expression data, genes, clusters and regulators from files
		M.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		M.readExpressionMatrix(data_file,null);
		M.readClusters(cluster_file);
		M.readRegulators(reg_file);
		M.initStatisticsAndScore();
		M.setDataMeanAndSDFromModuleset();

		// read regulation trees from xml file
		M.readRegTreeXML(tree_file);
		M.setTestSplits();
		// set top regulators for each module
		M.setTopRegulatorClasses(top_regulators);
		// calculate mean and sigma for all modules
		M.setModuleMeanSigma();
		M.checkExperiments();
		// use module mean (default) or global mean for figures
		M.setGlobalMeanForFigures(false);
		// use individual regulators mean for figures (default false)
		M.setRegulatorlMeanForFigures(use_regulator_mean);
		if (use_regulator_mean==true)
			M.setRegulatorMeanSigma();
		// change gene names if a map file is given
//		if (map_file != null)
//			M.changeGeneNames(map_file);
		// cut trees to a certain level
		if (cut_level > 0) {
			for (Module mod : M.moduleSet) {
				for (TreeNode t : mod.hierarchicalTrees) {
					t.testLevel(cut_level);
				}
			}
		}

		DrawModules dm = new DrawModules(M);
		dm.setPrefix("/Users/eric/gbm_test/");
		ArrayList<Integer> mod_list = new ArrayList<Integer>();
		mod_list.add(18);
		dm.customDraw(mod_list);
	}
}
