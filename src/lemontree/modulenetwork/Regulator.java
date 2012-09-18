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
 * Simple class to handle the regulators of a given module
 * 
 * @author eric
 *
 */
public class Regulator {
	private Gene gene;
	private double score;
	
	public Regulator(Gene g, double sc) {
		this.gene = g;
		this.score = sc;
	}
	
	public Gene getGene() {
		return(this.gene);
	}
	
	public Double getScore() {
		return(this.score);
	}
}
