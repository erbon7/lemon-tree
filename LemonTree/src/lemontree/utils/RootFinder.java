/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.utils;

import lemontree.modulenetwork.Globals; //revamp


/**
 * Finds a root of a real function of one variable in a given interval 
 * using the bisection method.
 * 
 * @author tomic
 *
 */
public class RootFinder {

	/**
	 * Function for which root has to be found.
	 */
	public RealFunction func;
	
	/**
	 * Lower bound of an interval containing the root.
	 */
	public double xMin;
	
	/**
	 * Upper bound of an interval containing the root.
	 */
	public double xMax;
	
	/**
	 * Accuracy of root.
	 */
	public double xAcc;
	
	/**
	 * The root.
	 */
	public double xRoot;
	
	/**
	 * Maximum number of iterations, regardless of accuracy
	 */
	private final int JMAX = 40;
	
	/**
	 * Maximum number of trials to find an xMax, given xMin
	 */
	private final int NTRY = 50;
	
	/**
	 * Geometric expansion factor to find an xMax, given xMin
	 */
	private final double FACTOR = 2;
	

	public RootFinder(){
		
	}
	
	public RootFinder(RealFunction func, double xMin, double xMax, double xAcc){
		this.func = func;
		this.xMin = xMin;
		this.xMax = xMax;
		this.xAcc = xAcc;
	}
	

	/**
	 * Bisection method to find the root.
	 * 
	 * @author tomic
	 */
	public void bisect() throws Exception {
		double dx, f, fmid, xmid, rt;
		f = func.evaluate(xMin);
		fmid = func.evaluate(xMax);
		//if (f*fmid >= 0.0)
		//	throw new Exception("Root must be bracketed for bisection.");
		// orient search so that f>0 lies at x+dx
		if (f < 0.0){
			rt = xMin;
			dx = xMax-xMin;
		} else {
			rt = xMax;
			dx = xMin-xMax;
		}
		// bisection loop
		for (int j=0; j<JMAX; j++){
			dx = 0.5*dx;
			xmid = rt+dx;
			fmid = func.evaluate(xmid);
			if (fmid <= 0.0)
				rt = xmid;
			if (Math.abs(dx) < xAcc || fmid == 0.0){
				xRoot = rt;
				return;
			}
		}
		// never get here
		throw new Exception("Too many bisections.");
	}
	
	/**
	 * Finds an xMax with opposite function sign as xMin, expanding geometrically.
	 * 
	 * @author tomic
	 */
	public void expandXMax() throws Exception {
		double f1, f2;
		f1 = func.evaluate(xMin);
		f2 = func.evaluate(xMax);
		for (int j=0; j<NTRY; j++){
			if (Globals.unweighted) {
				if ((f1*f2 < 0) || (Math.abs(f2) < xAcc)) 
					return;
			} else {
				if (f1*f2 < 0) // revamp: rm "|| (Math.abs(f2) < xAcc))"
					return;
			}
			xMin = xMax;
			xMax = FACTOR*xMax;
			f2 = func.evaluate(xMax);
		}
		// never get here
		throw new Exception("No xMax found.");
	}
}
