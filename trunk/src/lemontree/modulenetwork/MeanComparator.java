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
 * Class to compare the means of two leaf nodes
 * 
 * @author tomic
 *
 */
public class MeanComparator implements Comparator {
	
	public int compare(Object o1, Object o2) {
		Double mu1 = ((TreeNode)o1).leafDistribution.gaussianParam[0];
		Double mu2 = ((TreeNode)o2).leafDistribution.gaussianParam[0];
		return mu1.compareTo(mu2);
	}

}
