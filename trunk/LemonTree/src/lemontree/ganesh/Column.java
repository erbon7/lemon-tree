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

import cern.jet.stat.Gamma;

public class Column {
	public HashSet<Integer> ColumnSet;
	HashSet<Integer> RowSetTmp;
	public double count;
	public double sum;
	public double sum2;
	public double lambda0, alpha0, beta0, mu0;
	public double loglikelihood;
	
	public Column(double lambda0, double mu0, double alpha0, double beta0){
		this.lambda0 = lambda0;
		this.mu0 = mu0;
		this.alpha0 = alpha0;
		this.beta0 = beta0;
		RowSetTmp = new HashSet<Integer>();
	}

	/**
	 * loglikelihood of the column
	 * @param RowSet
	 * 					list of rows in the cluster
	 * @param data
	 * 					whole datamatrix.
     */
public double logLikelihood(HashSet<Integer> RowSet, double[][] data){
	count=0;
	sum=0;
	sum2=0;
	for(Integer row : RowSet)
	{
		for(Integer column : ColumnSet)
		{
			if(!Double.isNaN(data[row][column]))
			{
			count++;
			sum+=data[row][column];
			sum2+=data[row][column]*data[row][column];
			}
		}
	}
	this.loglikelihood=logLikelihoodFunction(count,sum,sum2);
	if(Double.isNaN(this.loglikelihood)){
		this.loglikelihood=0;
	}
	return(loglikelihood);
}

/**
 * loglikelihood sum of column when a row is removed from the cluster
 * @param data
 * 				whole datamatrix
 */		
public double loglikeRowminus(int row, double[][] data){
	for(Integer column : ColumnSet)
	{
		if(!Double.isNaN(data[row][column]))
		{
		count--;
		sum-=data[row][column];
		sum2-=data[row][column]*data[row][column];
		}
	}
	this.loglikelihood=logLikelihoodFunction(count,sum,sum2);
	if(Double.isNaN(this.loglikelihood))
		this.loglikelihood=0;
	return(loglikelihood);
}

/**
 * loglikelihood sum of column when a row is added to the cluster
 * @param data
 * 				whole datamatrix
 */		
public double loglikeRowplus(int row, double[][] data){
	for(Integer column : ColumnSet)
	{
		if(!Double.isNaN(data[row][column]))
		{
		count++;
		sum+=data[row][column];
		sum2+=data[row][column]*data[row][column];
		}
	}
	this.loglikelihood=logLikelihoodFunction(count,sum,sum2);
	if(Double.isNaN(this.loglikelihood))
		this.loglikelihood=0;
	return(loglikelihood);
}

/**
 * difference between loglikelihood(row in the column of cluster) 
 * and (loglikelihood(column of cluster) + loglikelihood(row) )
 * @param row
 * 				index of row to be moved
 * @param RowSet
 * 				list of rows in a cluster
 * @param data
 * 				whole datamatrix
 */		
public double logLikelihoodDiff(int row,HashSet<Integer> RowSet, double[][] data){
	double countsingle =0;
	double counttotal=0;
	double sumsingle=0;
	double sumtotal=0;
	double sum2single=0;
	double sum2total=0;
	double loglike = 0;
	for(Integer column : ColumnSet)
	{
		if(!Double.isNaN(data[row][column]))
		{
		countsingle++;
		sumsingle+=data[row][column];
		sum2single+=data[row][column]*data[row][column];
		}
	}
	counttotal=count+countsingle;
	sumtotal=sum+sumsingle;
	sum2total=sum2+sum2single;
	double like2=logLikelihoodFunction(counttotal,sumtotal,sum2total); // loglikelihood row in column of cluster
	double likesingle=logLikelihoodFunction(countsingle,sumsingle,sum2single); // loglikelihood of only row
	loglike=like2-(loglikelihood+likesingle);
//	loglike=like2/counttotal-(loglikelihood/count+likesingle/countsingle);
	if(Double.isNaN(loglike))
		loglike=0;
	return(loglike);
}

/**
 * difference between loglikelihood(row in the cluster) and loglikelihood(cluster)
 * @param row
 * 				index of row to be moved
 * @param RowSet
 * 				list of rows in a cluster
 * @param data
 * 				whole datamatrix
 */	
public double logLikelihoodDiffFixed(int row,HashSet<Integer> RowSet, double[][] data){
	double countsingle =0;
	double counttotal=0;
	double sumsingle=0;
	double sumtotal=0;
	double sum2single=0;
	double sum2total=0;
	double loglike = 0;
	for(Integer column : ColumnSet)
	{
		if(!Double.isNaN(data[row][column]))
		{
		countsingle++;
		sumsingle+=data[row][column];
		sum2single+=data[row][column]*data[row][column];
		}
	}
	counttotal=count+countsingle;
	sumtotal=sum+sumsingle;
	sum2total=sum2+sum2single;
	double like2=logLikelihoodFunction(counttotal,sumtotal,sum2total);
	loglike=like2-loglikelihood;
	if(Double.isNaN(loglike))
		loglike=0;
	return(loglike);
}

/**
 * likelihood of each column in a given cluster
 * @param N
 * 			number of datapoints in a column
 * @param sum
 * 				sum of values of N datapoints
 * @param sum2
 * 				sum of square of value of N datapoints
 */
public double logLikelihoodFunction(double N, double sum, double sum2){
	double loglikelihood = 0 ;
	double lambda1 = lambda0 + N;
	double alpha1 = alpha0 + 0.5 * N;
	double beta1 = beta0 + 0.5
			* (sum2 - Math.pow(sum, 2) / N)
			+ lambda0 
			         * Math.pow(sum - mu0 * N, 2)
			/ (2 * lambda1 * N);

	loglikelihood = -0.5 * N * Math.log(2 * Math.PI) + 0.5
			* Math.log(lambda0) + alpha0 * Math.log(beta0)
			- Gamma.logGamma(alpha0) + Gamma.logGamma(alpha1) - alpha1
			* Math.log(beta1) - 0.5 * Math.log(lambda1);
	return(loglikelihood);
}

}
