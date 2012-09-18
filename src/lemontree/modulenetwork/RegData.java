/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.modulenetwork;

/**
 * Simple class to store regulator informations and sort them
 * 
 * @author erbon
 *
 */
public class RegData implements Comparable<Object> {
	/**
	 * Regulator weight
	 */
	Double weight;
	
	/**
	 * Regulator name
	 */
	String name;
	
	/**
	 * Type: continuous or discrete
	 */
	String type;
	
	/**
	 * Constructor
	 * 
	 * @param w weight (double)
	 * @param g regulator (gene object)
	 * 
	 */
	public RegData(Double w, String n) {
		this.weight = w;
		this.name = n;
	}
	
	/**
	 * Constructor
	 * 
	 * @param w weight (double)
	 * @param g regulator (gene object)
	 * 
	 */
	public RegData(Double w, String n, boolean isDiscrete) {
		this.weight = w;
		this.name = n;
		if (isDiscrete)
			this.type = "d";
		else
			this.type = "c";
	}
	
	/**
	 * Comparator: compare RegData objects on weight values
	 */
	public int compareTo (Object o2) {
		int ret = 0;
		RegData r = (RegData) o2;
		if (weight < r.weight)
			ret = -1;
		else if (weight > r.weight)
			ret = 1;
		else if (weight.equals(r.weight))
			ret = 0;
		return ret;
	}
}
