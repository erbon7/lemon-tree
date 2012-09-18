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

public class Cluster {
	
	public int number; 
	public HashSet<Integer> RowSet;
	public List<Column> Columns;
	public int countTotal;
	public double sumTotal;
	public double sum2Total;
	public int countTotalTmp;
	public double sumTotalTmp;
	public double sum2TotalTmp;
	public double loglikelihood;
	double lambda0;
	double alpha0;
	double beta0;
	double mu0;
	
	public Cluster(double lambda0, double mu0, double alpha0, double beta0){
		this.lambda0=lambda0;
		this.mu0=mu0;
		this.alpha0=alpha0;
		this.beta0=beta0;
	}
	
	/**
     * randomly assigns columns to sqrt(num cluster) number of clusters
     * @param column
     * 					total number of columns 
     */
	public void RandomColumnAssign(int column){
		int num_cluster = (int)Math.sqrt(column);
		this.Columns = new ArrayList<Column>(num_cluster);
		for(int i=0;i<num_cluster;i++)
		{
			Column col = new Column(lambda0,mu0,alpha0,beta0);
			col.ColumnSet = new HashSet<Integer>();
			Columns.add(col);
		}
		for(int i=0;i<column;i++){
			int j = (int) Math.floor(Math.random()*num_cluster);
			Columns.get(j).ColumnSet.add(i);
			}	
	}
	
	/**
     * use Gibbs Sampler algorithm to cluster columns
     * @param column
     * 					total number of columns 
     * @param data
     * 					whole datamatrix
     */
	public void GibbsSampler(int column,double data[][]){
		double data1[][] = new double[column][RowSet.size()];
		int j=0;
		for(Integer row : RowSet)
		{
			for(int i=0;i<column;i++){
			data1[i][j]=data[row][i];
		}
			j++;
		}
		List<HashSet<Integer>> rows=new ArrayList<HashSet<Integer>>();
		HashSet<Integer> row;
		for(Column col : Columns){
			row=col.ColumnSet;
			rows.add(row);
		}
		GibbsSampler colCluster = new GibbsSampler(data1,rows);
		colCluster.AssignRowsSinglecolumn();
		for(int c=0;c<50;c++)
			colCluster.Cluster();
		this.Columns = new ArrayList<Column>(colCluster.ClusterSet.size());
		for(int i=0;i<colCluster.ClusterSet.size();i++)
		{
			Column col = new Column(lambda0,mu0,alpha0,beta0);
			col.ColumnSet = new HashSet<Integer>();
			col.ColumnSet = colCluster.ClusterSet.get(i).RowSet;
			Columns.add(col);
		}
	}
	
	/**
     * puts each columns in its own cluster
     * @param column
     * 					total number of columns 
     */	
	public void ColumnAssign(int column){
		this.Columns = new ArrayList<Column>(column);
		for(int j=0;j<column;j++)
			{
			Column col = new Column(lambda0,mu0,alpha0,beta0);
			col.ColumnSet = new HashSet<Integer>();
			col.ColumnSet.add(j);
			Columns.add(col);
			}	
	}
	
	/**
     * puts all columns in a single cluster
     * @param column
     * 					total number of columns 
     */	
	public void SingleColumnAssign(int column){
		this.Columns = new ArrayList<Column>();
		Column col = new Column(lambda0,mu0,alpha0,beta0);
		col.ColumnSet = new HashSet<Integer>();
		for(int j=0;j<column;j++)
			col.ColumnSet.add(j);
		Columns.add(col);
	}
	
	/**
     * loglikelihood sum for the cluster 
     * @param data
     * 				whole datamatrix 
     */		
	public double logLikelihood(double[][] data){
		loglikelihood = 0;
		countTotal=0;
		sumTotal=0;
		sum2Total=0;
		for(Column column : Columns ){
			double loglike = column.logLikelihood(RowSet,data);
			countTotal+=column.count;
			sumTotal+=column.sum;
			sum2Total+=column.sum2;
			loglikelihood+=loglike;
		}
		return(loglikelihood);
}
	
	/**
     * loglikelihood sum when a row is removed from the cluster
     * @param data
     * 				whole datamatrix
     */		
	public double loglikeRowminus(int row,double[][] data){
		double loglikelihood = 0;
		for(Column column : Columns ){
			double loglike = column.loglikeRowminus(row,data);
			loglikelihood+=loglike;
		}
		return(loglikelihood);
}

	/**
     * loglikelihood sum when a row is added to the cluster
     * @param data
     * 				 whole datamatrix
     */	
	public double loglikeRowplus(int row, double[][] data){
		double loglikelihood = 0;
		for(Column column : Columns ){
			double loglike = column.loglikeRowplus(row,data);
			loglikelihood+=loglike;
		}
		return(loglikelihood);
}
	
	/**
     * loglikelihood difference when a row is added to cluster
     * to when row is in its own cluster
     * @param data
     * 					whole datamatrix
     * @param row
     * 					index of row to be moved 
     */	
	public double logLikelihoodDiff(int row, double[][] data){
		double loglikelihood = 0;
		for(Column column : Columns ){
			double loglike = column.logLikelihoodDiff(row,RowSet,data);
			loglikelihood+=loglike;
		}
		return(loglikelihood);
}
	
	public double logLikelihoodDiffFixed(int row, double[][] data){
		double loglikelihood = 0;
		for(Column column : Columns ){
			double loglike = column.logLikelihoodDiffFixed(row,RowSet,data);
			loglikelihood+=loglike;
		}
		return(loglikelihood);
}

	
}
