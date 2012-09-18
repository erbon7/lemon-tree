/*
 * PFClustering - Java package for clustering weighted networks based on the Perron-Frobenius Theorem
 * 
 * Copyright (c) 2011 Tom Michoel
 *
 * This is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 * 
 */


package lemontree.pfclustering;
public class TestJSA{
	private double[][] values = {
					  {10.0,-2.0},
					  {3.0,9.0,3.0},
					  {7.0,8.0,7.0},
					  {3.0,8.0,7.0,5.0},
					  {8.0,9.0,9.0,13.0},
					  {4.0,2.0,-1.0}};

	private int[][] indexes = {
					   {0,4},
					   {0,1,5},
					   {1,2,3},
					   {0,2,3,4},
					   {1,3,4,5},
					   {1,4,5}};

	private double[] vector = {1.0,1.0,1.0,1.0,1.0,1.0};
	
	public static void main(String[] args){
		TestJSA testjsa = new TestJSA();
		testjsa.matrixvector();
	}

	public void matrixvector(){
		JavaSparseArray jsa = new JavaSparseArray(values,indexes);
		double[] vec = jsa.matrixvector(vector);
		write(vec);
	}

	public void write(double[] vec){
		for(int i = 0;i<vec.length;i++)System.out.println("value: "+vec[i]);
	}	
}

