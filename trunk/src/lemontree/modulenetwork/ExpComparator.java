/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.modulenetwork;

import java.util.Comparator;

/**
 * Class to compare to experiments (integer id) using the expression
 * values of a certain gene
 * 
 * @author tomic
 *
 */

public class ExpComparator implements Comparator {
	
	/**
	 * Module network where data can be found
	 */
	public ModuleNetwork moduleNetwork;
	
	/**
	 * Gene used for comparison
	 */
	public Gene gene;
	
	public ExpComparator(ModuleNetwork mn, Gene g){
		this.moduleNetwork = mn;
		this.gene = g;
	}
	
	
	public int compare(Object o1, Object o2) {
		Integer m1 = (Integer)o1;
		Integer m2 = (Integer)o2;
		Double x1 = this.moduleNetwork.data[this.gene.number][m1];
		Double x2 = this.moduleNetwork.data[this.gene.number][m2];
		return x1.compareTo(x2);
	}

}
