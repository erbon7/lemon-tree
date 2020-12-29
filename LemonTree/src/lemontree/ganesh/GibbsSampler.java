/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.ganesh;

import java.util.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;

public class GibbsSampler {

	public double FACTOR_0 = 1;
	Boolean flag = false;
	public double ROUND = 5;
	public double epsConv = 0.001;
	public double[][] data;
	public int row;
	public int column;
	public List<Cluster> ClusterSet;
	public double dataMean;
	public double dataSigma;
	public double lambda0 = 0.1, mu0 = 0.0, alpha0 = 0.1, beta0 = 0.1;
	public int num_cluster;
	List<HashSet<Integer>> rows;
	List<List<HashSet<Integer>>> partitions;
	public List<Row> AllRows;

	/**
	 * Method to read data for clustering and to calculate data mean and std.
	 * 
	 * @param data data matrix
	 */

	public GibbsSampler(double[][] data) {
		this.data = data;
		this.SetVariables();
	}

	/**
	 * Method to read data and initial clusters for clustering and to calculate
	 * data mean and std.
	 * 
	 * @param data data matrix
	 * @param rows initial clusters.
	 */
	public GibbsSampler(double[][] data, List<HashSet<Integer>> rows) {
		this.rows = rows;
		this.data = data;
		this.SetVariables();
	}

	/**
	 * Method to read data,initial clusters,partitions for clustering and to
	 * calculate data mean and std.
	 * 
	 * @param data data matrix
	 * @param rows initial clusters
	 * @param partitions clusters of columns for each row cluster
	 */
	public GibbsSampler(double[][] data, List<HashSet<Integer>> rows,
			List<List<HashSet<Integer>>> partitions) {
		this.rows = rows;
		this.partitions = partitions;
		this.data = data;
		this.SetVariables();
	}

	/**
	 * Performs a Gibbs sample run and returns cluster membership numbers for
	 * each row of the input matrix and for each sample. Since the output is
	 * used in Matlab, clusters are numbered starting from 1.
	 * 
	 * @param numSteps total number of steps for the procedure
	 * @param burnIn burn in steps
	 * @param sampleStep sample steps
	 * @return for each row, for each sample, the number of the cluster it belongs to
	 * @author tomic
	 */
	public double[][] Run(int numSteps, int burnIn, int sampleStep) {
		// compute number of samples
		int numSamples = 0;
		for (int k = 1; k <= numSteps; k++)
			if (k > burnIn && (numSteps - k) % sampleStep == 0)
				numSamples += 1;
		// initialize output matrix
		double[][] inCluster = new double[this.row][numSamples];
		// initialize sampler
		System.out.println("3");
		this.AssignRandomrowSinglecolumn();
		// start sampling
		int sampleCount = 0;
		System.out.println("4");
		for (int k = 1; k <= numSteps; k++) {
			this.Cluster();
			if (k > burnIn && (numSteps - k) % sampleStep == 0) {
				// store sample
				for (int j = 0; j < this.ClusterSet.size(); j++) {
					for (int r : this.ClusterSet.get(j).RowSet) {
						inCluster[r][sampleCount] = j + 1;
					}
				}
				// increase sample counter
				sampleCount += 1;
			}
		}

		return inCluster;
	}

	/**
	 * calculates data mean and std.
	 * 
	 * 2007-12-21 tomic : changed to loop only once over all data
	 */
	public void SetVariables() {
		double nrDataPoints = 0;
		double sumSquare = 0.0;
		for (int i = 0; i < this.data.length; i++) {
			for (int j = 0; j < this.data[i].length; j++) {
				if (!Double.isNaN(this.data[i][j])) {
					this.dataMean += this.data[i][j];
					sumSquare += Math.pow(this.data[i][j], 2);
					nrDataPoints += 1;
				}
			}
		}
		this.dataMean /= (double) nrDataPoints;
		this.dataSigma = Math.sqrt(sumSquare - nrDataPoints
				* Math.pow(this.dataMean, 2))
				/ Math.sqrt((double) nrDataPoints);

		this.row = this.data.length;
		this.column = this.data[0].length;
	}

	/**
	 * clusters data in both directions
	 */
	public void Cluster2way() {
		this.Cluster();
		for (Cluster clust : ClusterSet)
			clust.GibbsSampler(column, data);
	}

	/**
	 * moves each row once using the Gibbs sampler and merges clusters also
	 * using the Gibbs sampler
	 */
	public double Cluster() {
		AllRows = new ArrayList<Row>(row);
		for (int i = 0; i < row; i++) {
			Row newRow = new Row(i, num_cluster);
			AllRows.add(newRow);
		}
		double loglike = 0;
		loglike = this.loglikelihoodsum();
		for (int k = 0; k < row; k++) {
			int k1 = (int) (row * Math.random());
			this.operation(k1);
		}
		loglike = loglikelihoodsum();
		// merge clusters
    int j = 0;
    int total = ClusterSet.size();
		for (int m = 0; m < total; m++) {
			for (int k = 0; k < ClusterSet.size(); k++) {
				ClusterSet.get(k).number = k;
			}
			if (!this.operationCluster(ClusterSet.get(j))) {
        // the cluster was not removed
        // move the pointer forward
        j++;
      }
		}
		loglike = loglikelihoodsum();
		return (loglike);
	}

	/**
	 * checks whether each cluster would stay or merge into any other
	 */
	public double clusterClusters() {
		double large = 0, loglike = 0;
		loglike = this.loglikelihoodsum();
		//int count = 1;
		do {
			large = loglike;
		//	count += 1;
      int j = 0;
      int total = ClusterSet.size();
			for (int m = 0; m < total; m++) {
				for (int k = 0; k < ClusterSet.size(); k++) {
					ClusterSet.get(k).number = k;
				}
				if (!this.operationCluster(ClusterSet.get(j))) {
          // the cluster was not removed
          // move the pointer forward
          j++;
        }
			}
			loglike = loglikelihoodsum();
		} while (Math.abs(large - loglike) > epsConv * row * column);

		return (loglike);
	}

	/**
	 * calculates score difference when the two clusters clust1 and clust2 are
	 * merged compard to when seperate
	 * 
	 * @param clust1 first cluster
	 * @param clust2 second cluster
	 */
	public double scoreDiff(Cluster clust1, Cluster clust2) {
		double scoreDiff = 0;
		double totalscore = 0;
		Cluster clustTmp = new Cluster(lambda0, mu0, alpha0, beta0);
		clustTmp.RowSet = new HashSet<Integer>();
		for (Integer row : clust1.RowSet)
			clustTmp.RowSet.add(row);
		for (Integer row : clust2.RowSet)
			clustTmp.RowSet.add(row);
		clustTmp.SingleColumnAssign(column);
		totalscore = clustTmp.logLikelihood(data);
		scoreDiff = totalscore
				- (clust1.logLikelihood(data) + clust2.logLikelihood(data));
		return (scoreDiff);
	}

	/**
	 * moves each row once using the Gibbs sampler keeping the number of
	 * clusters fixed
	 */
	public double ClusterFixed() {
		AllRows = new ArrayList<Row>(row);
		for (int i = 0; i < row; i++) {
			Row newRow = new Row(i, num_cluster);
			AllRows.add(newRow);
		}
		double loglike;
		loglike = this.loglikelihoodsum();
		for (int k = 0; k < row; k++) {
			int k1 = (int) (row * Math.random());
			this.operationFixedclusters(k1);
		}
		loglike = this.loglikelihoodsum();
		return (loglike);
	}

	/**
	 * this function can be used when each cluster does not have different
	 * partitions like in case of clustering experiments where all genes go in
	 * one partition
	 * 
	 * @param row index of row to be moved
	 *            
	 */
	public void operationPartition(int row) {
		int currentCluster = 0;

		// find cluster containing given row
		for (int i = 0; i < num_cluster; i++) {
			if (ClusterSet.get(i).RowSet.contains(row)) {
				currentCluster = i;
				break;
			}
		}

		// remove cluster if the given row is only element
		if (ClusterSet.get(currentCluster).RowSet.size() == 1) {
			ClusterSet.remove(currentCluster);
			num_cluster = num_cluster - 1;
		}
		// else remove row from current cluster and update the statistics of
		// cluster in which the row is
		else {
			ClusterSet.get(currentCluster).RowSet.remove(row);
			ClusterSet.get(currentCluster).loglikeRowminus(row, data);
		}

		DoubleMatrix1D ratio = new DenseDoubleMatrix1D(num_cluster + 1);
		// create a vector with probability of gene to be in each cluster
		ratio.set(0, FACTOR_0);
		for (int i = 0; i < num_cluster; i++) {
			//int nsize = ClusterSet.get(i).RowSet.size();
			ratio.set(i + 1, bayesRatio(row, i));
		}
		double max_ratio = ratio.viewSorted().get(ratio.size() - 1);
		for (int j = 0; j < ratio.size(); j++) {
			ratio.set(j, Math.exp(ratio.get(j) - max_ratio));
		}
		double sum = ratio.zSum();
		for (int j = 0; j < ratio.size(); j++) {
			ratio.set(j, ratio.get(j) / sum);
		}
		AllRows.get(row).ratio = new DenseDoubleMatrix1D(num_cluster + 1);
		AllRows.get(row).ratio = ratio;
		// outcome is the no. of cluster to which the gene belongs
		// if 0 then the gene forms its own cluster and update statistics
		// otherwise add it to appro. cluster and update statistics
		int outcome = compare(ratio);
		if (outcome == 0) {
			num_cluster = num_cluster + 1;
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(ClusterSet.size() - 1).RowSet = new HashSet<Integer>();
			ClusterSet.get(ClusterSet.size() - 1).RowSet.add(row);
			if (flag)
				ClusterSet.get(ClusterSet.size() - 1)
						.RandomColumnAssign(column);
			else {
				ClusterSet.get(ClusterSet.size() - 1)
						.SingleColumnAssign(column);
			}
			ClusterSet.get(ClusterSet.size() - 1).logLikelihood(data);
		} else {
			int new_cluster = outcome - 1;
			ClusterSet.get(new_cluster).RowSet.add(row);
			ClusterSet.get(new_cluster).loglikeRowplus(row, data);
		}
	}

	/**
	 * function to find the suitable cluster for a given row
	 * 
	 * @param row index of the row to be moved.
	 *            
	 */
	public void operation(int row) {
		int currentCluster = 0;

		// find cluster containing given row
		for (int i = 0; i < num_cluster; i++) {
			if (ClusterSet.get(i).RowSet.contains(row)) {
				currentCluster = i;
				break;
			}
		}

		Cluster cluster1 = new Cluster(lambda0, mu0, alpha0, beta0);
		cluster1.RowSet = new HashSet<Integer>();
		cluster1.RowSet.add(row);
		cluster1.Columns = new ArrayList<Column>();
		for (Column col : ClusterSet.get(currentCluster).Columns) {
			Column col1 = new Column(lambda0, mu0, alpha0, beta0);
			col1.ColumnSet = new HashSet<Integer>();
			for (Integer i : col.ColumnSet)
				col1.ColumnSet.add(i);
			cluster1.Columns.add(col1);
		}

		// remove cluster if the given row is only element
		if (ClusterSet.get(currentCluster).RowSet.size() == 1) {
			// System.out.println(" cluster removed "+currentCluster);
			ClusterSet.remove(currentCluster);
			num_cluster = num_cluster - 1;
		}
		// else remove row from current cluster and update the statistics of
		// cluster in which the row is
		else {
			ClusterSet.get(currentCluster).RowSet.remove(row);
			ClusterSet.get(currentCluster).loglikeRowminus(row, data);
			// ClusterSet.get(currentCluster).logLikelihood(data);
		}

		DoubleMatrix1D ratio = new DenseDoubleMatrix1D(num_cluster + 1);
		// create a vector with probability of gene to be in each cluster
		ratio.set(0, FACTOR_0);
		double Singleloglike = cluster1.logLikelihood(data);
		for (int i = 0; i < num_cluster; i++) {
			//int nsize = ClusterSet.get(i).RowSet.size();
			// System.out.println(i+"\t"+bayesRatioFixedCluster(row,i));
			ratio.set(i + 1, bayesRatioFixedCluster(row, i) - Singleloglike);
		}
		double max_ratio = ratio.viewSorted().get(ratio.size() - 1);
		for (int j = 0; j < ratio.size(); j++) {
			ratio.set(j, Math.exp(ratio.get(j) - max_ratio));
		}
		double sum = ratio.zSum();
		for (int j = 0; j < ratio.size(); j++) {
			ratio.set(j, ratio.get(j) / sum);
		}
		AllRows.get(row).ratio = new DenseDoubleMatrix1D(num_cluster + 1);
		AllRows.get(row).ratio = ratio;
		// outcome is the no. of cluster to which the gene belongs
		// if 0 then the gene forms its own cluster and update statistics
		// otherwise add it to appro. cluster and update statistics

		int outcome = compare(ratio);

		if (outcome == 0) {
			num_cluster = num_cluster + 1;
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(ClusterSet.size() - 1).RowSet = new HashSet<Integer>();
			ClusterSet.get(ClusterSet.size() - 1).RowSet.add(row);
			if (flag)
				ClusterSet.get(ClusterSet.size() - 1)
						.RandomColumnAssign(column);
			else {
				ClusterSet.get(ClusterSet.size() - 1)
						.SingleColumnAssign(column);
			}
			ClusterSet.get(ClusterSet.size() - 1).logLikelihood(data);
		} else {
			int new_cluster = outcome - 1;
			// System.out.println(new_cluster);
			ClusterSet.get(new_cluster).RowSet.add(row);
			ClusterSet.get(new_cluster).loglikeRowplus(row, data);
			// ClusterSet.get(new_cluster).logLikelihood(data);
		}
	}

	/**
	 * operation without changing the number of clusters
	 * 
	 * @param row index of the row to be moved.
	 */
	public void operationFixedclusters(int row) {
		int currentCluster = 0;

		// find cluster containing given row
		for (int i = 0; i < num_cluster; i++) {
			if (ClusterSet.get(i).RowSet.contains(row)) {
				currentCluster = i;
				break;
			}
		}

		// remove cluster if the given row is only element
		if (ClusterSet.get(currentCluster).RowSet.size() == 1) {
			ClusterSet.remove(currentCluster);
			num_cluster = num_cluster - 1;
		}
		// else remove row from current cluster and update the statistics of
		// cluster in which the row is
		else {
			ClusterSet.get(currentCluster).RowSet.remove(row);
			ClusterSet.get(currentCluster).loglikeRowminus(row, data);
		}

		DoubleMatrix1D ratio = new DenseDoubleMatrix1D(num_cluster);
		// create a vector with probability of gene to be in each cluster
		if (currentCluster != num_cluster)
			ratio.set(currentCluster, Math.exp(ClusterSet.get(currentCluster)
					.logLikelihoodDiffFixed(row, data)));
		for (int i = 0; i < num_cluster; i++) {
			//int nsize = ClusterSet.get(i).RowSet.size();
			ratio.set(i, Math.exp(bayesRatioFixedCluster(row, i)));
		}
		double sum = ratio.zSum();
		for (int j = 0; j < ratio.size(); j++) {
			ratio.set(j, ratio.get(j) / sum);
		}
		int outcome = compare(ratio);
		int new_cluster = outcome;
		ClusterSet.get(new_cluster).RowSet.add(row);
		ClusterSet.get(new_cluster).loglikeRowplus(row, data);

	}

	/**
	 * operation to merge clusters
	 * 
	 * @param clust cluster to be moved.
	 */
	public boolean operationCluster(Cluster clust) {
		DoubleMatrix1D ratio = new DenseDoubleMatrix1D(ClusterSet.size());
		// create a vector of probability of merging given cluster with evry
		// other
		ratio.set(clust.number, 1);
		int i = 0;
		for (Cluster clust1 : ClusterSet) {
			if (clust1.number != clust.number)
				ratio.set(i, Math.exp(scoreDiff(clust, clust1)));
			i++;
		}
		double sum = ratio.zSum();
		for (int j = 0; j < ratio.size(); j++) {
			ratio.set(j, ratio.get(j) / sum);
		}
		// outcome is the no. of cluster to which give cluster merges
		int outcome = compare(ratio);
		if (outcome != clust.number) {
			for (Integer row : clust.RowSet)
				ClusterSet.get(outcome).RowSet.add(row);
			ClusterSet.remove(clust.number);
			num_cluster = num_cluster - 1;
      return true;
		}
    return false;
	}

	/**
	 * used with GibbsSampler(data, rows, partitions) stores initial clusters in
	 * both directions into proper variables
	 */
	public void InputAssign() {
		flag = true;
		num_cluster = rows.size();
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < rows.size(); i++) {
			for (Integer j : rows.get(i))
				ClusterSet.get(i).RowSet.add(j);
		}
		for (int i = 0; i < partitions.size(); i++) {
			int numColumn = partitions.get(i).size();
			ClusterSet.get(i).Columns = new ArrayList<Column>(numColumn);
			for (int j = 0; j < numColumn; j++) {
				Column col = new Column(lambda0, mu0, alpha0, beta0);
				col.ColumnSet = new HashSet<Integer>();
				for (Integer k : partitions.get(i).get(j))
					col.ColumnSet.add(k);
				ClusterSet.get(i).Columns.add(col);
			}
		}
	}

	/**
	 * used with GibbsSampler(data, rows) assigns rows according to the input
	 * for each row cluster, puts all the columns in one cluster
	 */
	public void AssignRowsSinglecolumn() {
		flag = false;
		num_cluster = rows.size();
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < rows.size(); i++) {
			for (Integer j : rows.get(i))
				ClusterSet.get(i).RowSet.add(j);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).SingleColumnAssign(column);
		}
	}

	/**
	 * used with GibbsSampler(data, rows) assigns rows according to the input
	 * for each row cluster, puts all the columns in one cluster
	 */
	public void AssignRowsEachcolumn() {
		flag = false;
		num_cluster = rows.size();
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < rows.size(); i++) {
			for (Integer j : rows.get(i))
				ClusterSet.get(i).RowSet.add(j);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).ColumnAssign(column);
		}
	}

	/**
	 * forms a single cluster of all the data
	 */
	public void AssignSinglerowSingleColumn() {
		flag = false;
		num_cluster = 1;
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
		ClusterSet.add(cluster);
		ClusterSet.get(0).RowSet = new HashSet<Integer>();
		for (int i = 0; i < row; i++) {
			ClusterSet.get(0).RowSet.add(i);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).SingleColumnAssign(column);
		}
	}

	/**
	 * used with GibbsSampler(data, rows) assigns rows according to the input
	 * for each row cluster, assigns columns randomly to sqrt(num. columns)
	 * clusters
	 */
	public void AssignRowsRandomcolumn() {
		flag = true;
		num_cluster = rows.size();
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < rows.size(); i++) {
			for (Integer j : rows.get(i))
				ClusterSet.get(i).RowSet.add(j);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).RandomColumnAssign(column);
		}
	}

	/**
	 * assigns rows randomly to sqrt(num. rows) clusters for each row cluster,
	 * assigns columns randomly to sqrt(num. columns) clusters
	 */
	public void RandomAssignTwoway(int num) {
		flag = true;
		// num_cluster = (int)Math.sqrt(this.row);
		if (num == 0)
			this.num_cluster = (int) Math.sqrt(this.row);
		else
			this.num_cluster = num;
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < row; i++) {
			int j = (int) Math.floor(Math.random() * num_cluster);
			ClusterSet.get(j).RowSet.add(i);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).RandomColumnAssign(column);
		}
	}

	/**
	 * assigns rows randomly to sqrt(num. rows) clusters for each row cluster,
	 * assigns columns to a single cluster
	 */
	public void AssignRandomrowSinglecolumn() {
		flag = false;
		num_cluster = (int) Math.sqrt(this.row);
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < row; i++) {
			int j = (int) Math.floor(Math.random() * num_cluster);
			ClusterSet.get(j).RowSet.add(i);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).SingleColumnAssign(column);
		}
	}

	/**
	 * assigns rows randomly to sqrt(num. rows) clusters for each row cluster,
	 * assigns each column to its own cluster
	 */
	public void AssignRandomrowEachcolumn() {
		flag = false;
		num_cluster = this.row;
		this.ClusterSet = new ArrayList<Cluster>(num_cluster);
		for (int i = 0; i < num_cluster; i++) {
			Cluster cluster = new Cluster(lambda0, mu0, alpha0, beta0);
			ClusterSet.add(cluster);
			ClusterSet.get(i).RowSet = new HashSet<Integer>();
		}
		for (int i = 0; i < row; i++) {
			ClusterSet.get(i).RowSet.add(i);
		}
		for (int i = 0; i < num_cluster; i++) {
			ClusterSet.get(i).ColumnAssign(column);
		}
	}

	/**
	 * total likelihood of finding data for given partition
	 */
	public double loglikelihoodsum() {
		double likelihood = 0;
		for (int i = 0; i < num_cluster; i++) {
			likelihood += ClusterSet.get(i).logLikelihood(data);
		}
		return (likelihood);
	}

	/**
	 * Gibbs sampler to determine to which cluster the row joins given an array
	 * of probabilities to join each cluster
	 * 
	 * @param ratio array of probability to join the cluster
	 * @return the cluster chosen to join
	 */
	public int compare(DoubleMatrix1D ratio) {
		int result = 0;
		double compare = Math.random();
		double partialsum = 0;
		for (int j = 0; j < ratio.size(); j++) {
			partialsum = partialsum + ratio.get(j);
			if (partialsum > compare) {
				result = j;
				break;
			}

		}
		return (result);
	}

	/**
	 * ratio of log of probability that row joins num_cluster cluster to row
	 * forming its own cluster
	 * 
	 * @param row index of the row to be moved
	 * @param num_cluster index of cluster to which row joins
	 * @return ratio
	 */
	public double bayesRatio(int row, int num_cluster) {
		//double n = ClusterSet.get(num_cluster).RowSet.size();
		double ratio = ClusterSet.get(num_cluster).logLikelihoodDiff(row, data);
		return (ratio);
	}

	public double bayesRatioFixedCluster(int row, int i) {
		double ratio2 = ClusterSet.get(i).logLikelihoodDiffFixed(row, data);
		return (ratio2);
	}

	public double bayesRatioFixedClusterOnegene(int row, int currentCluster, int i) {
		Cluster clust = new Cluster(lambda0, mu0, alpha0, beta0);
		clust.RowSet = new HashSet<Integer>();
		clust.RowSet.add(row);
		if (flag)
			clust.SingleColumnAssign(column);
		else
			clust.ColumnAssign(column);
		double ratio1 = clust.logLikelihood(data);
		double ratio2 = ClusterSet.get(i).logLikelihoodDiffFixed(row, data);
		return (ratio1 - ratio2);
	}
}
