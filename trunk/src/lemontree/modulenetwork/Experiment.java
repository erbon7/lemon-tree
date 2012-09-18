/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.modulenetwork;

import java.awt.Color;
import java.util.*;

/**
 * Class to store description and ID number for an experiment. The ID number
 * refers to the column index in the data matrix.
 * 
 * @author tomic
 *
 */
public class Experiment implements Comparable {

	/**
	 * Description of experiment
	 */
	public String name;
	
	/**
	 * ID number of experiment. This refers to the column index in the data matrix.
	 */
	public int number;
	
	/**
	 * Tissue 
	 */
	public String tissue;
	
	/**
	 * State (ex. normal / tumor)
	 */
	public String state;
	
	/** 
	 * A color object associated to each condition
	 */
	public Color col;
	
	/**
	 * a simple integer store an experiment ranking value.
	 */
	public int rank;
	
	public Experiment(){
		
	}
	
	public Experiment(String name, int number){
		this.name = name;
		this.number = number;
	}
	
	public int compareTo(Object o) {
		Experiment e = (Experiment) o;
		if (this.rank < e.rank)
			return -1;
		if (this.rank > e.rank)
			return +1;
		if (this.rank == e.rank)
			return 0;
		return 0;
	}
	
}
