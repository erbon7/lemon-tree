package lemontree.test;

import java.io.IOException;

import lemontree.modulenetwork.ModuleNetwork;
import java.io.*;

import lemontree.utils.DrawModules;

public class TestDrawMod {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		String datasetFile ="/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/data.all.merged";
//		String regulatorsFile="/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_data_merged";
//		String regulatorsType="/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_types";
//		String clusterFile="/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/n6_agilent_TC10.txt";
//		String regFile = "/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_st_tf_ki/st.tf.kin.list";
//		String exp = "/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_st_tf_ki/exp.xml.gz";
		
		String prefix = "/bioinfo/users/ebonnet/TreeLemon/tcga/figures/";
		String datasetFile = prefix + "data_matrix";
		String regulatorsFile= prefix + "reg_files_list";
		String regulatorsType=prefix + "reg_types";
		String clusterFile=prefix + "U133A_tc10.txt";
		String regFile = prefix + "st.tf.kin.list";
		String exp = prefix + "exp.xml.gz";
		
		String reg_param = prefix + "reg_param";
		String extra_info = prefix + "inter_U133A_AGI_fc_sign";
		
		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile,null);
		m.readClusters(clusterFile);
		
		m.readRegulators(regFile);
		
		m.readRegTreeXML(exp);
		m.initStatisticsAndScore();
//		m.setTopRegulators("/Users/eric/reg.txt");
		
		//m.setTopRegulatorClasses(regulatorsFile, regulatorsType);
		m.setTopRegulatorClasses(regulatorsFile);
		
		
		m.setModuleMeanSigma();
		m.setRegulatorMeanSigma();
		m.checkExperiments();
		m.setRegulatorlMeanForFigures(true);
		m.setGlobalMeanForFigures(true);
		
		DrawModules dm = new DrawModules(m);
		dm.setPrefix(prefix + "test/");
		
		
		// read extra parameters for each type of regulator
//		try {
//			BufferedReader br = new BufferedReader(new FileReader(reg_param));
//			String line;
//			System.out.println("Adding reg parameters...");
//			while((line = br.readLine()) != null) {
//				line.trim();
//				String[] tk = line.split("\\s+");
//				System.out.println(tk[0]+" "+tk[1]+" "+tk[2]);
//				dm.regClassMean.add(Double.parseDouble(tk[1]));
//				dm.regClassSd.add(Double.parseDouble(tk[2]));
//			}
//			br.close();
//			
//			
//			br = new BufferedReader(new FileReader(extra_info));
//			while((line = br.readLine()) != null) {
//				line.trim();
//				String[] tk = line.split("\\t");
//				//System.out.println(tk[0]+ " => "+tk[1]);
//				dm.extraInfo.put(tk[0], Integer.parseInt(tk[1]));
//			}
//			br.close();
//		}
//		catch (Exception e) {
//			e.printStackTrace();
//		}
		
		//dm.drawAllModules();
		dm.drawOneModule(0);
	}

}
