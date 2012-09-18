/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.utils;

import java.util.*;
import cern.colt.matrix.impl.*;

/**
 * A class to perform simple operations on set (union, intersection, etc.)
 * 
 * @author tomic
 *
 */

public class SetOperation  {

	/**
	 * Converts an ArrayList to a HashSet
	 * 
	 * @author tomic
	 */
    public static <E> HashSet<E> list2set (ArrayList<E> list){
    	HashSet<E> set = new HashSet<E>();
    	for (E i : list)
    		set.add(i);
    	return set;
    }
    
    /**
     * Converts a HashSet to an ArrayList
     * 
     * @author tomic
     */
    public static <E> ArrayList<E> set2list (Set<E> set){
    	ArrayList<E> list = new ArrayList<E>();
    	for (E i : set)
    		list.add(i);
    	return list;
    }
    
    /**
     * Returns the intersection of two sets.
     * 
     * @param <E>
     * @param list1
     * @param list2
     * @return
     * 
     * @author tomic
     */
	public static <E> ArrayList<E> intersection(ArrayList<E> list1, ArrayList<E> list2){
		ArrayList<E> intsect = new ArrayList<E>();
		for (E i : list1)
			if (list2.contains(i))
				intsect.add(i);
		return intsect;
	}
	
	/**
	 * Returns the union of two sets
	 * 
	 * @param <E>
	 * @param list1
	 * @param list2
	 * @return
	 * 
	 * @author tomic
	 */
	public static <E> ArrayList<E> union(ArrayList<E> list1, ArrayList<E> list2){
		ArrayList<E> union= new ArrayList<E>();
		for (E i : list1)
			union.add(i);
		for (E j : list2)
			if (!union.contains(j))
				union.add(j);
		return union;
	}
	
	/**
	 * Returns the product partition of two partitions, i.e., the partition consisting of all 
	 * pairwise intersections.
	 * 
	 * @author tomic
	 */
	public static <E> ArrayList<ArrayList<E>> prodPart (ArrayList<ArrayList<E>> part1, 
				ArrayList<ArrayList<E>> part2){
		ArrayList<ArrayList<E>> intpart = new ArrayList<ArrayList<E>>();
		ArrayList<E> intsect = new ArrayList<E>();
		for (ArrayList<E> list1 : part1)
			for (ArrayList<E> list2 : part2){
				intsect = intersection(list1, list2);
				if (!intsect.isEmpty())
					intpart.add(intsect);
			}
		return intpart;
	}

	/**
	 * The relative sizes of each set in a partition.
	 * 
	 * @author tomic
	 */
	public static <E> DenseDoubleMatrix1D probPart (ArrayList<ArrayList<E>> part) {
		DenseDoubleMatrix1D prob = new DenseDoubleMatrix1D(part.size());
		int total = 0;
		for (ArrayList<E> list : part)
			total += list.size();
		for (int k=0; k<part.size(); k++){
			prob.set(k, (double)part.get(k).size()/(double)total);
		}
		return prob;
	}
	
	
	/**
	 * Conditional probabilities of partition partM given partR.
	 * 
	 * @author tomic
	 */
	public static <E> DenseDoubleMatrix2D condProbPart (ArrayList<ArrayList<E>> partR, 
			ArrayList<ArrayList<E>> partM) {
		DenseDoubleMatrix2D condprob = new DenseDoubleMatrix2D(partR.size(), partM.size());
		
		for (int k=0; k<partR.size(); k++)
			for (int l=0; l<partM.size(); l++)
				condprob.set(k, l, (double)intersection(partR.get(k),partM.get(l)).size()
						/(double)partR.get(k).size());
		
		return condprob;
	}
	

	/**
	 * Entropy of probability distribution.
	 * 
	 * @author tomic
	 */ 
	public static double entropy (DenseDoubleMatrix1D prob) {
		double ent = 0;
		for (int k=0; k<prob.size(); k++)
			ent += hEnt(prob.get(k));
		return ent;
	}
	
	/** 
	 * conditional entropy given probabilities
	 * 
	 * @param prob
	 * @param condprob
	 * @return
	 * 
	 * @author tomic
	 */
	public static double condEntropy (DenseDoubleMatrix1D prob, DenseDoubleMatrix2D condprob) {
		double ent = 0;
		for (int k=0; k<prob.size(); k++){
			double sum = 0;
			for (int l=0; l<condprob.columns(); l++)
				sum += hEnt(condprob.get(k,l));
			ent += prob.get(k) * sum;
		}
		return ent;
	}
	
	/** 
	 * entropy function h(p)=-p*log(p)
	 * 
	 * @param p
	 * @return
	 * 
	 * @author tomic
	 */
    public static double hEnt(double p){
        double plogp = 0.0;
        if (p > 0.0 && p < 1.0)
            plogp = -p*Math.log(p);
        return plogp;
    }
    

}
