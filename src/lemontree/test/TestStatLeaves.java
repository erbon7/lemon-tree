package lemontree.test;


import java.io.*;
import java.text.*;
import java.util.*;
import cern.jet.stat.Probability;
import BiNGO.HypergeometricDistribution;
import lemontree.modulenetwork.*;


public class TestStatLeaves {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		
		String datasetFile ="/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/data.all.merged";
		String clusterFile="/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/n6_agilent_TC10.txt";
		String regFile = "/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_st_tf_ki/st.tf.kin.list";
		String exp = "/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_st_tf_ki/exp.xml.gz";
		//String reg = "/bioinfo/users/ebonnet/TreeLemon/n6/agilent_mRNA/reg_st_tf_ki/reg.xml.gz";
		
		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1f, 0.0, 0.1f, 0.1f);
		m.readExpressionMatrix(datasetFile,null);
		m.readClusters(clusterFile);
		
		m.readRegulators(regFile);
		m.initStatisticsAndScore();

		m.readRegTreeXML(exp);
		//m.readRegTreeXML(reg);
		m.initStatisticsAndScore();

		m.setConditionMap();
		annotateExperiments(m);
		
		try {
			statLeaves(m);
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}


	/*
	 * Set a tissue type for each condition (= each sample, = each column)
	 */
	public static void annotateExperiments (ModuleNetwork modnet) {

		// experiments type first col = sample name, second = tissue type 
		String expFile = "/bioinfo/users/ebonnet/TreeLemon/n6/experiments.txt";

		// annotate experiments
		try {
			BufferedReader input = new BufferedReader(new FileReader(new File(expFile)));
			String line;
			while((line = input.readLine()) != null) {
				line.trim();
				String[] tk = line.split("\t");
				String exp_name = tk[0];
				String tissue = tk[1];
				modnet.conditionMap.get(exp_name).tissue = tissue;
			}
			input.close();
		}
		catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}

	}

	static void statLeaves(ModuleNetwork m) throws IOException {

		// Number of samples, = conditions = columns in the matrix
		int bigN = 60;

		// tissue categories
		HashSet<String> tissueSet = new HashSet<String>();
		tissueSet.add("BR");
		tissueSet.add("CO");
		tissueSet.add("CNS");
		tissueSet.add("LC");
		tissueSet.add("LE");
		tissueSet.add("ME");
		tissueSet.add("OV");
		tissueSet.add("PR");
		tissueSet.add("RE");

		// count for each tissue type in complete conditions set
		HashMap<String,Integer> smallN = new HashMap<String,Integer>();
		smallN.put("LC", 9);
		smallN.put("CO", 7);
		smallN.put("OV", 7);
		smallN.put("CNS", 6);
		smallN.put("RE", 8);
		smallN.put("ME", 10);
		smallN.put("BR", 5);
		smallN.put("LE", 6);
		smallN.put("PR", 2);

		
		// output file
		BufferedWriter out = new BufferedWriter(new FileWriter("/bioinfo/users/ebonnet/TreeLemon/n6/output.txt"));
		
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMaximumFractionDigits(2); 
		nf.setMinimumFractionDigits(2);
		
		for (Module mod : m.moduleSet) {

			// for each tree
			int numTree = 0;
			for (TreeNode root : mod.hierarchicalTrees) {

				// list of leaves (internal nodes + leaves)
				ArrayList<TreeNode> nodeList = new ArrayList<TreeNode>();
				root.gatherInternals(nodeList);
				//root.gatherLeafList(leafList);

				// for each leave 
				for (int i=0; i<nodeList.size();i++) {
					TreeNode t = nodeList.get(i);
					if (t.nodeStatus.equals("internal")) {

						TreeNode tr = t.rightChild;
						TreeNode tl = t.leftChild;
						
						// gather information about the right and left conditions
						double[] mean = new double[2];
						double[] sd = new double[2];
						double[] N = new double[2];
						
						mean[0] = t.rightChild.leafDistribution.gaussianParam[0];
						sd[0] = t.rightChild.leafDistribution.gaussianParam[1];
						N[0] = t.rightChild.leafDistribution.statistics[0];
						
						mean[1] = t.leftChild.leafDistribution.gaussianParam[0];
						sd[1] = t.leftChild.leafDistribution.gaussianParam[1];
						N[1] = t.leftChild.leafDistribution.statistics[0];
						
						StringBuffer[] str = new StringBuffer[2];
						str[0] = new StringBuffer();
						str[1] = new StringBuffer();
						
						for (int j=0;j < t.rightChild.leafDistribution.condSet.size(); j++) {
							Experiment e = m.conditionSet.get(t.rightChild.leafDistribution.condSet.get(j));
							str[0].append(e.tissue+":");
						}
						
						for (int j=0;j < t.leftChild.leafDistribution.condSet.size(); j++) {
							Experiment e = m.conditionSet.get(t.leftChild.leafDistribution.condSet.get(j));
							str[1].append(e.tissue+":");
						}
						
						ArrayList<TreeNode> the_childs = new ArrayList<TreeNode>();
						the_childs.add(tr);
						the_childs.add(tl);

						int ct=0;
						for (TreeNode child : the_childs) {
							
							if (child.leafDistribution.condSet.size()>5) {

								HashMap<String,Integer> tissueCount = new HashMap<String, Integer>();
								
								// initiate count
								for (String type : tissueSet)
									tissueCount.put(type, 0);

								// Nb of conditions in the leave
								int leafSize = child.leafDistribution.condSet.size();

								// Count tissue categories in the leave
								for (int j=0;j<leafSize;j++) {
									// Experiment ID
									int expNum = child.leafDistribution.condSet.get(j);
									// tissue code
									String tissue = m.conditionSet.get(expNum).tissue;
									// increment tissue count
									int c = tissueCount.get(tissue);
									tissueCount.put(tissue, c+1);
								}
								
								int second_child = ct==0 ? 1 : 0;
								
								int leafSizeSecondChild = the_childs.get(second_child).leafDistribution.condSet.size();
								
								// calculate enrichment for each tissue category
								for (String ts : tissueSet) {
									int smallX = tissueCount.get(ts);
									if (leafSize >= 5 && leafSizeSecondChild >= 5 && smallX > 2) {
										HypergeometricDistribution hg = new HypergeometricDistribution(smallX,leafSize,smallN.get(ts),bigN);
										//double prob = hg.calculate();
										double prob = Double.parseDouble(hg.calculateHypergDistr());
										
										double oddsRatio = (smallX * 1.0 / leafSize * 1.0) / (smallN.get(ts) * 1.0 / bigN * 1.0);
										
										// write results to output file
										/* Fields:
										 * 0: module number
										 * 1: tree number
										 * 2: tree node number
										 * 3: tissue type
										 * 4: proportion of tissue in child leaf : leaf size
										 * 5: hypergeometric test p-value
										 * 6: odds ratio
										 * 7: mean first child
										 * 8: standard dev. first child
										 * 9: mean second child
										 * 10: standard dev second child
										 * 11: t-test p-value
										 * 12: string of tissue types for first child
										 * 13: string of tissue type for second child
										 * 
										 */
										out.write(mod.number+"\t"+numTree+"\t"+i+"\t"+ts+"\t"+smallX+":"+leafSize+"\t"+prob+"\t"+oddsRatio+"\t");
										out.write(nf.format(mean[ct])+"\t"+nf.format(sd[ct])+"\t"+nf.format(mean[second_child])+"\t"+nf.format(sd[second_child])+"\t");
										
										double[] res = new double[3];
										welchTest(N[ct],mean[ct],sd[ct],N[second_child],mean[second_child],sd[second_child],res);
										out.write(res[0]+"\t");
										
										// strings of tissues labels
										out.write(str[ct]+"\t"+str[second_child]);
										out.newLine();
										
									}
								}
							}
							ct++;
						}
					}
				}
				++numTree;
			}
		}

		out.close();
	}
	
	static void testWelchTest () {
		
		double N1 = 10.0;
		double N2 = 14.0;
		double mean1 = 5.5;
		double mean2 = 13.5;
		double sd1 = 3.02765;
		double sd2 = 4.1833;
		double[] res = new double[3];
		
		welchTest(N1,mean1,sd1,N2,mean2,sd2,res);
		
		System.out.println("p-value: "+res[0]);  // 1.855e-5
		System.out.println("t-statistic: "+res[1]); // -5.4349
		System.out.println("degrees of freedom: "+res[2]); // 21.982
	}
	
	/**
	 *  Welch t-test for 2 samples of unequal size and variances
	 *  
	 * @param N1
	 * @param mean1
	 * @param sd1
	 * @param N2
	 * @param mean2
	 * @param sd2
	 * @param res results array 0:p-value, 1:t-statistic, 2:degrees of freedom
	 */
	static void welchTest (double N1, double mean1, double sd1, double N2, double mean2, double sd2, double[] res) {
		
		// variance
		double var1 = Math.pow(sd1, 2);
		double var2 = Math.pow(sd2, 2);
		
		// t-statistic 
		res[1] = (mean1 - mean2) / Math.sqrt(var1/N1 + var2/N2);
		
		// degrees of freedom
		res[2] = Math.pow(var1/N1+var2/N2, 2);
		res[2] = res[2] / (Math.pow(var1,2)/(Math.pow(N1, 2)*(N1-1)) + Math.pow(var2,2)/(Math.pow(N2, 2)*(N2-1))); 
		
		// p-value, two-tailed
		if (res[1] > 0)
			res[0] = (1 - Probability.studentT(res[2],res[1])) * 2;
		else 
			res[0] = Probability.studentT(res[2],res[1]) * 2;
	}

}
