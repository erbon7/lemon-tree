/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.modulenetwork;

import java.io.IOException;

import lemontree.utils.BiNGO;
import lemontree.utils.DrawModules;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import lemontree.modulenetwork.Globals; //revamp


/**
 * 
 * Run LeMoNe with command line options
 * 
 * @author erbon
 *
 */
public class RunCli {
	

	/**
	 * Parse command line options and run LeMoNe 
	 * 
	 * @param args command-line arguments string
	 *            
	 */
	public static void main(String[] args) {
		
		// set dummy values for those parameters, they'll be filled later
		String task = null;
		String gene_file = null;
		String data_file = null;
		String reg_file = null;
		String tree_file = null;
		String output_file = null;
		String range = null;
		int num_steps = 0;
		int burn_in = 0;
		int sample_steps = 0;
		String cluster_file = null;
		String column_file = null; //revamp
		String go_annot_file = null;
		String go_ref_file = null;
		String top_regulators = null;
		String map_file = null;
		String go_ontology_file = null;
		String draw_experiment_color = null;
		
		
		// set default values for those parameters, users can override them
		double alpha = 0.1;
		double beta = 0.1;
		double mu = 0.0;
		double lambda = 0.1;
		double score_gain = 0.0;
		int init_num_clust = 0;
		int num_runs = 1;
		boolean use_bayesian_score = true;
		int num_reg = 10;
		double beta_reg = 20;
		String go_p_value = "0.05";
		String go_namespace = "biological_process";
		boolean use_global_mean = false;
		boolean use_regulator_mean = false;
		int cut_level = 0;
		double min_weight = 0.25;
		int min_clust_size = 10;
		int min_clust_score = 2;
		Boolean node_clustering = true;
		boolean draw_experiment_names = true;
		double reassign_thr = 0.0;

		// create the different options
		Options opts = new Options();
		opts.addOption("task", true, "task to perform");
		opts.addOption("gene_file", true, "gene file");
		opts.addOption("data_file", true, "data file (genes)");
		opts.addOption("reg_file", true, "regulators file");
		opts.addOption("tree_file", true, "tree file");
		opts.addOption("output_file", true, "output file name");
		opts.addOption("num_steps", true, "number of steps (Gibbs sampler)");
		opts.addOption("burn_in", true, "number of burn-in steps (Gibbs sampler)");
		opts.addOption("sample_steps", true, "sample steps interval (Gibbs sampler)");
		opts.addOption("cluster_file", true, "cluster file name");
		opts.addOption("num_clust", true, "number of clusters");
		opts.addOption("column_file", true, "column weight file name"); //revamp
		opts.addOption("alpha", true, "alpha0 parameter value");
		opts.addOption("beta", true, "beta0 parameter value");
		opts.addOption("mu", true, "mu0 parameter value");
		opts.addOption("lambda", true, "lambda0 parameter value");
		opts.addOption("score_gain", true, "score gain cutoff value");
		opts.addOption("init_num_clust", true, "initial number of clusters (Gibbs sampler)");
		opts.addOption("num_runs", true, "number of runs (Gibbs sampler)");
		opts.addOption("num_reg", true, "maximum number of regulators assigned for each node");
		opts.addOption("beta_reg", true, "beta parameter value for regulators assignment");
		opts.addOption("num_split",true, "number of splits for the module set");
		opts.addOption("prefix",true, "java command prefix for the split option command line");
		opts.addOption("range",true, "module set range for the assignment of the regulators");
		opts.addOption("go_annot_file", true, "GO custom annotation file");
		opts.addOption("go_ontology_file", true, "GO ontology file name");
		opts.addOption("go_ref_file", true, "GO refence gene list file name");
		opts.addOption("go_p_value", true, "GO p-value cutoff");
		opts.addOption("go_namespace", true, "GO namespace");
		opts.addOption("go_annot_def", false, "GO annotation file flag");
		opts.addOption("matlab", false, "Matlab format for output files");
		opts.addOption("help", false, "help");
		opts.addOption("h", false, "help");
		opts.addOption("top_regulators", true, "Top regulators file name");
		opts.addOption("use_global_mean", false, "Use global mean for the figures");
		opts.addOption("use_regulator_mean", false, "Use regulator mean for the figures");
		opts.addOption("all_regulators", false, "Print all regulators");
		opts.addOption("map_file", true, "Gene names map file");
		opts.addOption("cut_level", true, "Regulation tree cut level");
		opts.addOption("min_weight", true, "Tight clusters minimum weight");
		opts.addOption("min_clust_size", true, "Tight clusters minimum cluster size");
		opts.addOption("min_clust_score", true, "Tight clusters minimum cluster score");
		opts.addOption("node_clustering", true, "Perform node clustering (true) or edge clustering (false)");
		opts.addOption("draw_experiment_names", true, "Draw experiment names in the figures");
		opts.addOption("draw_experiment_color", true, "Draw experiment color codes in the figures");
		opts.addOption("reassign_thr", true, "Node re-assignment threshold");

		// build a parser object and parse the command line (!)
		CommandLineParser parser = new PosixParser();
		try {
			
			CommandLine cmd = parser.parse(opts, args);
			if (cmd.hasOption("min_weight"))
				min_weight = Double.parseDouble(cmd.getOptionValue("min_weight"));
			
			if (cmd.hasOption("min_clust_size"))
				min_clust_size = Integer.parseInt(cmd.getOptionValue("min_clust_size"));
			
			if (cmd.hasOption("min_clust_score"))
				min_clust_score = Integer.parseInt(cmd.getOptionValue("min_clust_score"));
			
			if (cmd.hasOption("task"))
				task = cmd.getOptionValue("task");
			
			if (cmd.hasOption("data_file"))
				data_file = cmd.getOptionValue("data_file");
			
			if (cmd.hasOption("tree_file"))
				tree_file = cmd.getOptionValue("tree_file");
			
			if (cmd.hasOption("gene_file"))
				gene_file = cmd.getOptionValue("gene_file");
			
			if (cmd.hasOption("reg_file"))
				reg_file = cmd.getOptionValue("reg_file");
			
			if (cmd.hasOption("output_file"))
				output_file = cmd.getOptionValue("output_file");
			
			if (cmd.hasOption("cluster_file")) //revamp
				cluster_file = cmd.getOptionValue("cluster_file");

			if (cmd.hasOption("column_file"))
				column_file = cmd.getOptionValue("column_file");

			if (cmd.hasOption("alpha"))
				alpha = Double.parseDouble(cmd.getOptionValue("alpha"));
			
			if (cmd.hasOption("beta"))
				beta = Double.parseDouble(cmd.getOptionValue("beta"));
			
			if (cmd.hasOption("lambda"))
				lambda = Double.parseDouble(cmd.getOptionValue("lambda"));
			
			if (cmd.hasOption("mu"))
				mu = Double.parseDouble(cmd.getOptionValue("mu"));
			
			if (cmd.hasOption("score_gain"))
				score_gain = Double.parseDouble(cmd.getOptionValue("score_gain"));
			
			if (cmd.hasOption("num_steps"))
				num_steps = Integer.parseInt(cmd.getOptionValue("num_steps"));
			
			if (cmd.hasOption("burn_in"))
				burn_in = Integer.parseInt(cmd.getOptionValue("burn_in"));
			
			if (cmd.hasOption("sample_steps"))
				sample_steps = Integer.parseInt(cmd.getOptionValue("sample_steps"));
			
			if (cmd.hasOption("init_num_clust"))
				init_num_clust = Integer.parseInt(cmd.getOptionValue("init_num_clust"));
			
			if (cmd.hasOption("num_reg"))
				num_reg = Integer.parseInt(cmd.getOptionValue("num_reg"));
			
			if (cmd.hasOption("beta_reg"))
				beta_reg = Double.parseDouble(cmd.getOptionValue("beta_reg"));
			
			if (cmd.hasOption("range"))
				range = cmd.getOptionValue("range");
			
			if (cmd.hasOption("go_annot_file"))
				go_annot_file = cmd.getOptionValue("go_annot_file");
			
			if (cmd.hasOption("go_ontology_file"))
				go_ontology_file = cmd.getOptionValue("go_ontology_file");
			
			if (cmd.hasOption("go_ref_file"))
				go_ref_file = cmd.getOptionValue("go_ref_file");
			
			if (cmd.hasOption("go_p_value"))
				go_p_value = cmd.getOptionValue("go_p_value");

			if (cmd.hasOption("go_namespace"))
				go_namespace = cmd.getOptionValue("go_namespace");
			
			if (cmd.hasOption("top_regulators"))
				top_regulators = cmd.getOptionValue("top_regulators");
			
			if (cmd.hasOption("use_global_mean"))
				use_global_mean = true;
			
			if (cmd.hasOption("use_regulator_mean"))
				use_regulator_mean = true;
			
			if (cmd.hasOption("map_file"))
				map_file = cmd.getOptionValue("map_file");
			
			if (cmd.hasOption("cut_level"))
				cut_level = Integer.parseInt(cmd.getOptionValue("cut_level"));
			
			if (cmd.hasOption("node_clustering"))
				if (cmd.getOptionValue("node_clustering").equalsIgnoreCase("false"))
					node_clustering = false;
			
			if (cmd.hasOption("draw_experiment_names"))
				if (cmd.getOptionValue("draw_experiment_names").equalsIgnoreCase("false"))
					draw_experiment_names = false;
			
			if (cmd.hasOption("draw_experiment_color"))
				draw_experiment_color = cmd.getOptionValue("draw_experiment_color");
			
			if (cmd.hasOption("reassign_thr"))
				reassign_thr = Double.parseDouble(cmd.getOptionValue("reassign_thr"));

		}
		catch (ParseException exp) {
			System.out.println("Error while parsing command line:");
			System.out.println();
			exp.printStackTrace();
			System.exit(1);
		}
		
		
		// print header
		printBanner();
		
		// something has to be done, we need a task to be set
		if (task == null)
			Die("Error: task option must be set.");
		
		// ---------------------------------------------------------------
		// ganesh task: 2-way clustering of genes using the gibbs sampler
		// ---------------------------------------------------------------
		if (task.equalsIgnoreCase("ganesh")) {

			// those parameters must be set
			if (data_file == null)
				Die("Error: data_file option must be set.");
			if (output_file == null) 
				Die("Error: output_file option must be set.");
			
			// set default values if the user did not change them
			if (num_runs == 0)
				num_runs = 1;
			if (num_steps == 0)
				num_steps = 100;
			if (burn_in == 0)
				burn_in = 50;
			if (sample_steps == 0)
				sample_steps = 100;

			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:               " + task);
			System.out.println("data_file:          " + data_file);
			System.out.println("gene_file:          " + gene_file);
			System.out.println("output_file:        " + output_file);
			System.out.println("lambda:             " + lambda);
			System.out.println("mu:                 " + mu);
			System.out.println("alpha:              " + alpha);
			System.out.println("beta:               " + beta);
			System.out.println("num_steps:          " + num_steps);
			System.out.println("burn_in:            " + burn_in);
			System.out.println("sample_steps:       " + sample_steps);
			System.out.println("score_gain:         " + score_gain);
			
			// Create ModuleNetwork object
			ModuleNetwork M = new ModuleNetwork();
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			M.readExpressionMatrix(data_file, gene_file);
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			// Gibbs sample different module sets with one tree per module
			M.gibbsSamplerGenes(init_num_clust, num_runs, num_steps, burn_in, sample_steps, score_gain, use_bayesian_score);
			// write results to text file
			M.writeClusters(output_file);
		 }
		 //-------------------------------------------------------------------------------------
		 // tight_clusters task: node clustering to produce tight clusters
		 //-------------------------------------------------------------------------------------
		 else if (task.equalsIgnoreCase("tight_clusters")) {
			 
			 if (data_file == null)
					Die("Error: data_file option must be set.");
			 if (cluster_file == null)
				 Die("Error: cluster_file option must be set.");
			 if (output_file == null)
				 Die("Error: output_file option must be set.");
			 
			 System.out.println("Parameters");
			 System.out.println("----------");
			 System.out.println("task:               " + task);
			 System.out.println("data_file:          " + data_file);
			 System.out.println("cluster_file:       " + cluster_file);
			 System.out.println("output_file:        " + output_file);
			 System.out.println("node_clustering:    " + node_clustering);
			 System.out.println("min_weight:         " + min_weight);
			 System.out.println("min_clust_size:     " + min_clust_size);
			 System.out.println("min_clust_score:    " + min_clust_score);
			 System.out.println();
			 
			 ModuleNetwork M = new ModuleNetwork();
			 M.readExpressionMatrix(data_file, null);
			 M.readMultipleClusters(cluster_file);

			 // find tight clusters with node clustering algorithm
			 CentroidClustering cc = new CentroidClustering(M, node_clustering, min_weight, min_clust_size, min_clust_score);
			 cc.doCentroidClustering();
			 cc.printClusters(output_file);
			 
		 }
		 //-------------------------------------------------------------------------------------
		 // regulators task: learn regulation programs (gibbs sampling exp. + assign regulators)
		 //-------------------------------------------------------------------------------------
		 else if (task.equalsIgnoreCase("regulators")) {
			
			// those parameters must be set
			if (data_file == null)
				Die("Error: data_file option must be set.");
			if (reg_file == null)
				Die("Error: reg_file option must be set.");
			if (cluster_file == null)
				Die("Error: cluster_file option must be set.");
			if (output_file == null)
				Die("Error: output_file option must be set.");
			
			// set default values if the user did not change them
			if (num_steps == 0)
				num_steps = 1100;
			if (burn_in == 0)
				burn_in = 100;
			if (sample_steps == 0)
				sample_steps = 100;

			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:               " + task);
			System.out.println("data_file:          " + data_file);
			System.out.println("reg_file:           " + reg_file);
			System.out.println("cluster_file:       " + cluster_file);
			System.out.println("output_file:        " + output_file);
			System.out.println("lambda:             " + lambda);
			System.out.println("mu:                 " + mu);
			System.out.println("alpha:              " + alpha);
			System.out.println("beta:               " + beta);
			System.out.println("num_runs:           " + num_runs);
			System.out.println("num_steps:          " + num_steps);
			System.out.println("burn_in:            " + burn_in);
			System.out.println("sample_steps:       " + sample_steps);
			System.out.println("score_gain:         " + score_gain);
			System.out.println("num_reg:            " + num_reg);
			
			// create module network object
			ModuleNetwork M = new ModuleNetwork();
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			M.readExpressionMatrix(data_file, null);
			M.readClusters(cluster_file);
			M.readRegulators(reg_file);
			M.initStatisticsAndScore();
			M.setDataMeanAndSDFromModuleset();

			// cluster experiments using the gibbs sampler
			M.gibbsSamplerExpts(num_runs, num_steps, burn_in, sample_steps, score_gain, use_bayesian_score);
			// assign regulators
			M.assignRegulatorsNoAcyclStoch(beta_reg, num_reg);
			// write results as text file with all regulators, top 1%, random regulators and regulations trees as xml
			M.printRegulators(output_file+".allreg.txt", true, false);
			M.printRegulators(output_file+".topreg.txt", false, false);
			M.printRandomRegulators(output_file+".randomreg.txt", false);
			M.writeRegTreeXML(output_file+".xml.gz");
		} 
		//----------------------------------------------------------
		// experiments task: cluster conditions using gibbs sampling 
		//----------------------------------------------------------
		 else if (task.equalsIgnoreCase("experiments")) {
			
			// those parameters must be set
			if (data_file == null)
				Die("Error: data_file option must be set.");
			if (cluster_file == null)
				Die("Error: cluster_file option must be set.");
			if (output_file == null)
				Die("Error: output_file option must be set.");
			
			// set default values if the user did not change them
			if (num_steps == 0)
				num_steps = 1100;
			if (burn_in == 0)
				burn_in = 100;
			if (sample_steps == 0)
				sample_steps = 100;

			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:               " + task);
			System.out.println("data_file:          " + data_file);
			System.out.println("cluster_file:       " + cluster_file);
			System.out.println("output_file:        " + output_file);
			System.out.println("lambda:             " + lambda);
			System.out.println("mu:                 " + mu);
			System.out.println("alpha:              " + alpha);
			System.out.println("beta:               " + beta);
			System.out.println("num_steps:          " + num_steps);
			System.out.println("burn_in:            " + burn_in);
			System.out.println("sample_steps:       " + sample_steps);
			System.out.println("score_gain:         " + score_gain);
			
			
			// read data and clusters
			ModuleNetwork M = new ModuleNetwork();
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			M.readExpressionMatrix(data_file, null);
			M.readClusters(cluster_file);
			M.initStatisticsAndScore();

			// cluster experiments using the gibbs sampler
			M.gibbsSamplerExpts(num_runs, num_steps, burn_in, sample_steps, score_gain, use_bayesian_score);

			// write results as xml file
			M.writeRegTreeXML(output_file);
		}
		//---------------------------------------------------------------
		// split_reg task: assign regulators for a given range of modules
		//---------------------------------------------------------------
		else if (task.equalsIgnoreCase("split_reg")) {
			
			// those parameters must be set
			if (data_file == null)
				Die("Error: data_file option must be set.");
			if (reg_file == null)
				Die("Error: reg_file option must be set.");
			if (cluster_file == null)
				Die("Error: cluster_file option must be set.");
			if (tree_file == null)
				Die("Error: tree_file option must be set.");
			if (output_file == null)
				Die("Error: output_file option must be set.");
			if (range == null)
				Die("Error: range option must be set.");
			
			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:               " + task);
			System.out.println("data_file:          " + data_file);
			System.out.println("reg_file:           " + reg_file);
			System.out.println("cluster_file:       " + cluster_file);
			System.out.println("output_file:        " + output_file);
			System.out.println("num_reg:            " + num_reg);
			System.out.println("beta_reg:           " + beta_reg);
			System.out.println("range:              " + range);
			
			ModuleNetwork M = new ModuleNetwork();
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			M.readExpressionMatrix(data_file, null);
			M.readClusters(cluster_file);
			M.readRegulators(reg_file);
			M.initStatisticsAndScore();
			M.readRegTreeXML(tree_file);
			
			String[] val = range.split(":");
			int start_module = Integer.parseInt(val[0]);
			int stop_module = Integer.parseInt(val[1]);
			
			// assign regulators
			M.assignRegulatorsNoAcyclStoch(beta_reg, num_reg, start_module, stop_module);
	        
			// write results
			M.printRegulators(output_file+".allreg.txt", true, false);
			M.printRandomRegulators(output_file+".randomreg.txt", false);
			M.writeRegTreeXML(output_file+".xml.gz");
		}
		//----------------------------------------------------------------------------
		// go_annotation task: GO annotation of a cluster file
		//----------------------------------------------------------------------------
		else if (task.equalsIgnoreCase("go_annotation")) {
			
			// those parameters must be set
			if (cluster_file == null)
				Die("Error: cluster_file parameter must be set.");
			if (output_file == null)
				Die("Error: output_file option must be set.");
			if (go_annot_file == null)
				Die("Error: go_annot_file option must be set.");
			if (go_ontology_file == null)
				Die("Error: go_ontology_file option must be set.");
			if (go_ref_file == null)
				Die("Error: go_ref_file option must be set.");
			
			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:               " + task);
			System.out.println("cluster_file:       " + cluster_file);
			System.out.println("output_file:        " + output_file);
			System.out.println("go_annot_file:      " + go_annot_file);
			System.out.println("go_ontology_file:   " + go_ontology_file);
			System.out.println("go_ref_file:        " + go_ref_file);
			System.out.println("go_p_value:         " + go_p_value);
			System.out.println("go_namespace:       " + go_namespace);
			System.out.println("map_file            " + map_file);
			
			BiNGO b = new BiNGO(go_annot_file, go_ontology_file, go_p_value, go_namespace);
			
			try {
				b.GOstats(cluster_file, go_ref_file, output_file, map_file);
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
		//----------------------------------------------------------------------------
		// figures task: create eps figures for each module
		//----------------------------------------------------------------------------
		else if (task.equalsIgnoreCase("figures")) {
			
			// those parameters must be set
			if (top_regulators == null)
				Die("Error: top_regulators option must be set.");
			if (data_file == null)
				Die("Error: data_file option must be set.");
			if (reg_file == null)
				Die("Error: reg_file option must be set.");
			if (cluster_file == null)
				Die("Error: cluster_file option must be set.");
			if (tree_file == null)
				Die("Error: tree_file option must be set.");

			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:                  " + task);
			System.out.println("data_file:             " + data_file);
			System.out.println("reg_file:              " + reg_file);
			System.out.println("cluster_file:          " + cluster_file);
			System.out.println("tree_file:             " + tree_file);
			System.out.println("top_regulators:        " + top_regulators);
			System.out.println("use_regulator_mean:    " + use_regulator_mean);
			System.out.println("use_global_mean:       " + use_global_mean);
			System.out.println("map_file:              " + map_file);
			System.out.println("cut_level:             " + cut_level);
			System.out.println("draw_experiment_names: " + draw_experiment_names);
			System.out.println("draw_experiment_color: " + draw_experiment_color);
			
			ModuleNetwork M = new ModuleNetwork();
			//read expression data, genes, clusters and regulators from files
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
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
			M.setGlobalMeanForFigures(use_global_mean);
			// use individual regulators mean for figures (default false)
			M.setRegulatorlMeanForFigures(use_regulator_mean);
			if (use_regulator_mean==true)
				M.setRegulatorMeanSigma();
			// change gene names if a map file is given
			if (map_file != null)
		          M.changeGeneNames(map_file);
			 // cut trees to a certain level
	        if (cut_level > 0) {
	          for (Module mod : M.moduleSet) {
	            for (TreeNode t : mod.hierarchicalTrees) {
	              t.testLevel(cut_level);
	            }
	          }
	        }

	       
	        DrawModules dm = new DrawModules(M);
	        
	        if (draw_experiment_color != null) {
				M.setExperimentColor(draw_experiment_color);
				dm.enableExperimentColor();
	        }
			
			if (draw_experiment_names == false) {
				dm.unsetDrawExperimentNames();
			}
			
			dm.drawAllModules();
		}
		//----------------------------------------------------------------------------
		// topdown task: run "old" heuristic algo
		//----------------------------------------------------------------------------
		else if (task.equalsIgnoreCase("topdown")) {
			int maxParents = 3;
			double epsConvergence = 1E-3;
			// Create ModuleNetwork object
			ModuleNetwork M = new ModuleNetwork();
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			M.readExpressionMatrix(data_file, gene_file);
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			M.readRegulators(reg_file);
			// Top-down search
			M.heuristicSearchMaxTopDown(maxParents,epsConvergence);
			// write results as xml file
			M.writeRegTreeXML(output_file);
		}

		//----------------------------------------------------------------------------
		// revamp task: reassign nodes if score is above a certain threshold
		//----------------------------------------------------------------------------
		else if (task.equalsIgnoreCase("revamp")) {

			double scoregain = 0.0;
			boolean useBHCscore = true;
			double epsConvergence = 1E-3;
			

			// those parameters must be set
			if (data_file == null)
				Die("Error: data_file parameter must be set.");
			if (cluster_file == null)
				Die("Error: cluster_file parameter must be set.");
			if (output_file == null)
				Die("Error: output_file option must be set.");
			

			System.out.println("Parameters");
			System.out.println("----------");
			System.out.println("task:               " + task);
			System.out.println("data_file:          " + data_file);
			System.out.println("cluster_file:       " + cluster_file);
			System.out.println("column_file:        " + column_file); //revamp
			System.out.println("output_file:        " + output_file);
			System.out.println("node_clustering:    " + node_clustering);
			System.out.println("reassign_thr:       " + reassign_thr);


			// Create ModuleNetwork object and initialise to into one cluster
			ModuleNetwork M = new ModuleNetwork();
			M.setNormalGammaPriors(lambda, mu, alpha, beta);
			// if a "gene_file" list is given, add only the genes that are in the list
			M.readExpressionMatrixRevamp(data_file, cluster_file, column_file); //revamp
			M.readClusters(cluster_file);
			//M.readRegulators(reg_file);

			// run node reassignment
			M.heuristicSearchMaxRevamp(reassign_thr, scoregain, useBHCscore, epsConvergence);

			// write results as xml file
			M.writeRegTreeXML(output_file+".xml.gz");
			//M.writeXML(output_file);

			// write results to text file
			M.writeClusters(output_file+".txt");
			

		}

		else {
			System.out.println("task option '"+task+"' unknown.");
			System.out.println();
		}
	}
	
	public static void Die (String msg) {
		System.out.println(msg);
		System.exit(1);
	}
	
	public static void printBanner () {
		System.out.println("");
		System.out.println("LemonTree - Version 3.1.0");
		System.out.println("-------------------------");
		System.out.println("");
	}
}
