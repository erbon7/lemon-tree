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

/**
  * This class implements a sparse matrix as 
  * Java Sparse Arrays.
  * @author Geir Gundersen
  * @version 1.0
  */
  public class JavaSparseArray{
	/**The 2D array that stores the values*/
  	private double[][] Avalue;
	/**The 2D array that stores the indexes*/
  	private int[][] Aindex;
	/**The dimension of the matrix*/
  	private int dimension;
  	
	private int m;
  	private int n;
   	
  	/**Constructor for the class
  	  */
  	  public JavaSparseArray(double[][] Avalue,int[][] Aindex){
  	  	this.Avalue = Avalue;
  	  	this.Aindex = Aindex;
  	  	dimension = Avalue.length;
  	  }
  	  
  	  /**A method for setting the value array
	  @param double[][] value
  	  */
  	  public void setValueArray(double[][] value){Avalue = value;}
  	  
  	  /**A method for setting the index array
	  @param int[][] index
  	  */
  	  public void setIndexArray(int[][] index){Aindex = index;} 
  	
  	
	 /**A method for getting the value array
	 @return double[][]
  	  */
  	  public double[][] getValueArray(){return Avalue;}
  	  
  	  /**A method for getting the index array
	  @return double[][]
  	  */
  	  public int[][] getIndexArray(){return Aindex;} 
	  
	  /**Returns an element from the matrix at position (i,j)
	  @param int i, int j
   	  @return double element
	  */
	  public double elementAt(int i, int j){
	  	double element = 0.0;
		boolean test = true;
		double[] value = Avalue[i];
		int[] index = Aindex[i];
		for(int k = 0;k<value.length&&test;k++){
			if(j==index[k]){
				test = false;
				element = value[k];
			}
		}
		return element;
	  }
	  
	  /**Sets an element from the matrix at position (i,j)
	  @param int i, int j, double v
	  */
	  public void setElementAt(int i, int j, double v){
	  	double element = 0.0;
		boolean test = true;
		double[] value = Avalue[i];
		int[] index = Aindex[i];
		for(int k = 0;k<value.length&&test;k++){
			if(j==index[k]){
				test = false;
				element = value[k];
			}
		}
		if(test)
			addElementAt(i, j, v);
		else
			value[j] = v;
	  }
	  
	  /**Removes an element from the matrix at position (i,j)
	  @param int i, int j, double v
	  */
	  public void removeElementAt(int i, int j){
	  	int stop = 0;
	  	double element = 0.0;
		boolean test = true;
		double[] value = Avalue[i];
		int[] index = Aindex[i];
		for(int k = 0;k<value.length&&test;k++){
			if(j==index[k]){
				test = false;
				stop = k;
				element = value[k];
			}
		}
		double[] value1 = new double[value.length-1];
		int[] index1 = new int[value.length-1];
		System.arraycopy(value, 0, value1, 0, stop-1);
		System.arraycopy(index, 0, index1, 0, stop-1);
		System.arraycopy(value, stop+1, value1, stop, value1.length);
		System.arraycopy(index, stop+1, index1, stop, value1.length);
		Avalue[i] = value1;
		Aindex[i] = index1;
	  }
	  
	  /**Adds an element at the matrix at position (i,j)
	  @param int i, int j, double v
	  */
	  public void addElementAt(int i, int j, double v){
	  	int stop = 0;
	  	double element = 0.0;
		boolean test = true;
		double[] value = Avalue[i];
		int[] index = Aindex[i];
		for(int k = 0;k<value.length&&test;k++){
			if(j==index[k]){
				test = false;
				element = value[k];
			}
		}
		double[] value1 = new double[value.length+1];
		int[] index1 = new int[value.length+1];
		System.arraycopy(value, 0, value1, 0, stop-1);
		System.arraycopy(index, 0, index1, 0, stop-1);
		value1[stop] = v;
		index1[stop] = j;
		System.arraycopy(value, stop, value1, stop+1, value1.length);
		System.arraycopy(index, stop, index1, stop+1, value1.length);
		Avalue[i] = value1;
		Aindex[i] = index1;
	  }

  	  /**A method for multiplying a matrix with a vector returning a vector
	  @param double[] b
   	  @return double[] c
  	  */
  	  public double[] matrixvector(double[] b){
  	  	double x = 0.00;
  	  	double[] value;
  	  	int[] index;
  	  	int k = 0;
  	  	int vlength = 0;
  	  	double[] c = new double[b.length];
  	  	int alength = Avalue.length;
  	  	for(int i=0;i<alength;i++){
  	  		value = Avalue[i];
  	  		index = Aindex[i];
  	  		for(int j=0;j<value.length;j++){
  	  			x += value[j]*b[index[j]];
  	  		}
  	  		c[i]=x;
  	  		x = 0;
  	  	}
    	  	return c;
  	   }
  	   
	/**A method for multiplying a vector with a matrix returning a vector
	@param double[] bvalue
   	@return double[] cvalue
  	*/   
	public double[] vectormatatrix(double[] bvalue){
		double[] valuerow = null;
		int[] indexrow = null;
		int alength = Avalue.length;
		double[] cvalue = new double[bvalue.length];
		double value = 0;
		for(int i = 0;i<alength;i++){
			valuerow = Avalue[i];
			indexrow = Aindex[i];
			double val = bvalue[i];
			int vrow = valuerow.length;
			for(int j = vrow-1;j>0;--j){
				cvalue[indexrow[j]] += value*valuerow[j];
			}
		}
		return cvalue;
	}

	/**A method for multiplying two matrices returning the result matrix
	@param JavaSparseArray B
   	@return JavaSparseArray C
  	*/   
	public JavaSparseArray times(JavaSparseArray B){
		double[][] Cvalue = new double[dimension][1];
		int[][] Cindex = new int[dimension][1];
		int[] temp = new int[dimension];
		double[] tempValue = new double[dimension];
		int[] tempIndex = new int[dimension];
		double[][] Bvalue = B.getValueArray();
		int[][] Bindex = B.getIndexArray();
		double scalar = 0;
		int len = -1;
		int index = 0;
		int jcol = 0;
		int jpos = 0;
		int nonzero = 0;
		double[] avalue;
		int[] aindex;
		double[] bvalue;
		int[] bindex;
		double[] cvalue;
		int[] cindex;
		for(int i = 0;i<temp.length;i++){temp[i]=-1;}
		for(int i = 0;i<Avalue.length;i++){
			avalue = Avalue[i];
			aindex = Aindex[i];
			for(int j = 0;j<avalue.length;j++){
				scalar = avalue[j];
				index = aindex[j];
				bvalue = Bvalue[index];
				bindex = Bindex[index];
				for(int k = 0;k<bvalue.length;k++){
					jcol = bindex[k];
					jpos = temp[jcol];
					if(jpos == -1){
						len++;
						nonzero++;
						tempIndex[len] = jcol;
						temp[jcol] = len;
						tempValue[len] = scalar*bvalue[k];

					}else{
						tempValue[jpos]+=scalar*bvalue[k];
					}
				}
			}
			cvalue = new double[len+1];
			cindex = new int[len+1];
			System.arraycopy(tempValue, 0, cvalue, 0, len+1);
			System.arraycopy(tempIndex, 0, cindex, 0, len+1);
			Cvalue[i] = cvalue;
			Cindex[i] = cindex;
			for(int ii = 0;ii<len+1;ii++){temp[tempIndex[ii]]=-1;}
			len = -1;
		}
		return new JavaSparseArray(Cvalue, Cindex);
	}
	
	/**A method for adding two matrices returning the result matrix
	@param JavaSparseArray B
   	@return JavaSparseArray C
  	*/      
	public JavaSparseArray addition(JavaSparseArray B){
		int[][] Bindex = B.getIndexArray();
		double[][] Bvalue = B.getValueArray();
		int[][] Cindex = new int[dimension][1];
		double[][] Cvalue = new double[dimension][1];
		int nonzero = 0;
		int num = 0;
		for(int i = 0;i<dimension;i++){
			int[] bindex = Bindex[i];
			double[] bvalue = Bvalue[i];
			int[] aindex = Aindex[i];
			double[] avalue = Avalue[i];
			boolean[] switchArray = new boolean[dimension];
			double[] tempValue = new double[dimension];
			for(int ii = 0;ii<aindex.length;ii++){
				switchArray[aindex[ii]] = true;
				tempValue[aindex[ii]] = avalue[ii];
			}
			num = aindex.length;
			for(int ii = 0;ii<bindex.length;ii++){
				if(!switchArray[bindex[ii]]){
					switchArray[bindex[ii]] = true;
					tempValue[bindex[ii]] = bvalue[ii];
					num++; 
				}
				else{
					tempValue[bindex[ii]] += bvalue[ii];
				}
			}
			int[] cindex = new int[num];
			double[] cvalue = new double[num];
			for(int ii = 0;ii<aindex.length;ii++){
				cindex[ii] = aindex[ii];
				switchArray[aindex[ii]] = false;
			}
			for(int ii = aindex.length,jj=0;ii<bindex.length;ii++,jj++){
				if(switchArray[bindex[jj]]){
					cindex[ii] = bindex[jj];
				}
			}
			for(int ii = 0;ii<cvalue.length;ii++){
				cvalue[ii] = tempValue[cindex[ii]];
			}
			Cindex[i] = cindex;
			Cvalue[i] = cvalue;
			nonzero += num;
		}
		return new JavaSparseArray(Cvalue, Cindex);	
	}

	/**A method for updating a matrix returning
	@param JavaSparseArray A
   	@return JavaSparseArray A (updated)
  	*/   
	public JavaSparseArray update(JavaSparseArray B){
		int[][] Bindex = B.getIndexArray();
		double[][] Bvalue = B.getValueArray();
		double[] tempValue = new double[dimension];
		int[] tempIndex = new int[dimension];
		int[] temp = new int[dimension];
		int len = -1;
		int index = 0;
		int jcol = 0;
		int jpos = 0;
		int nonzero = 0;
		int num = 0;
		for(int i = 0;i<temp.length;i++){temp[i] = -1;}
		for(int i = 0;i<dimension;i++){
			int[] aindex = Aindex[i];
			double[] avalue = Avalue[i];
			if(Bindex[i]!=null){
			for(int ii = 0;ii<aindex.length;ii++){
				len++;
				jcol = aindex[ii];
				tempIndex[len] = aindex[ii];
				temp[jcol] = len;
				tempValue[len] = avalue[ii];
			}
			int[] bindex = Bindex[i];
			double[] bvalue = Bvalue[i];
			for(int ii = 0;ii<bindex.length;ii++){
				jcol = bindex[ii];
				jpos = temp[jcol];
				if(jpos == -1){
					len++;
					tempIndex[len] = jcol;
					temp[jcol] = len;
					tempValue[len] = bvalue[ii]; 
				}
				else{
					tempValue[jpos] += bvalue[ii];
				}
			}
			int[] cindex = new int[len+1];
			double[] cvalue = new double[len+1];
			System.arraycopy(tempValue, 0, cvalue, 0,len+1);
			System.arraycopy(tempIndex, 0, cindex, 0,len+1);
			Aindex[i] = cindex;
			Avalue[i] = cvalue;
			for(int ii = 0;ii<len+1;ii++){temp[tempIndex[ii]] = -1;}
			len = -1;
			}
		}
		return new JavaSparseArray(Avalue, Aindex);	
	}
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
