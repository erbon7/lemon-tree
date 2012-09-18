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
 * Desribes the methods a class must have to run a simulated annealing procedure.
 * 
 * @author tomic
 * 
 */
public interface SimAnnealConfiguration {
	
	/**
	 * Replaces the current configuration by a random configuration.
	 * 
	 * @return new configuration.
	 */
	SimAnnealConfiguration randomInit();
	
	/**
	 * Makes a random change to the current configuration.
	 * 
	 * @return new configuration.
	 */
	SimAnnealConfiguration randomChange();
	
	/**
	 * Computes the function to be minimized.
	 * 
	 * @return value.
	 */
	double energyFunction();
	
}
