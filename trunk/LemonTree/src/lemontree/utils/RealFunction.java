/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.utils;

/**
 * An interface for real functions of one variable.
 * 
 * @author tomic
 *
 */

public interface RealFunction {
	/**
	 * Returns function value at argument x.
	 */
	double evaluate (double x);
}
