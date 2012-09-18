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
 * A simple class for weighted edges
 */

public class Edge {

	public String from;
	public String to;
	public double weight;
	public boolean issym; // TRUE for symmetric (undirected) edges
	
	public Edge (){
		
	}
	
	public Edge (String from, String to){
		this.from = from;
		this.to = to;
		this.weight = 1.0;
		this.issym = false;
	}
	
	public Edge (String from, String to, double weight){
		this.from = from;
		this.to = to;
		this.weight = weight;
		this.issym = false;
	}
	
	public Edge (String from, String to, double weight, boolean issym){
		this.from = from;
		this.to = to;
		this.weight = weight;
		this.issym = issym;
	}
	
	public boolean equals (Edge e){
		if (this.issym && e.issym){
			if ( (this.from.equals(e.from) && this.to.equals(e.to)) || (this.from.equals(e.to) && this.to.equals(e.from)))
				return true;
			else
				return false;
		} else {
			if (this.from.equals(e.from) && this.to.equals(e.to) && this.issym==e.issym)
				return true;
			else
				return false;
		}
	}
	
	/**
	 * Return inverse edge
	 */
	public Edge inverse (){
		return new Edge (this.to, this.from, this.weight, this.issym);
	}
	
	public String toString(){
		return this.from + "\t" + this.to + "\t" + this.weight;
	}
	
}
