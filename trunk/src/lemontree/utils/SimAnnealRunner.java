/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.utils;

import cern.colt.list.*;
import cern.jet.stat.*;
import java.util.ArrayList;
import java.io.*;
import java.util.Date;

/**
 * Runs a simulated annealing for combinatorial optimization.
 * 
 * To use: define a class implementing the interface {@link SimAnnealConfiguration} 
 * and pass an instance of that class as an argument to the method 
 * {@link SimAnnealRunner#run}. 
 * 
 * @author tomic
 *
 */

public class SimAnnealRunner {
	
	/**
	 * Initial acceptance probability.
	 */
	public double pInit; 
	
	/**
	 * Number of transitions at each temperature.
	 */
	public int nTrans;
	
	/**
	 * Number of acceptances at each temperature.
	 */
	public int nAccept;
	
	/**
	 * Number of temperature values.
	 */
	public int nT;
	
	/**
	 * Number of temperature steps that minimum can be constant before aborting.
	 */
	public int nUnchanged;
	
	/**
	 * Cooling rate. <em>T<sub>k</sub> = &alpha;*T<sub>k-1</sub></em>
	 */
	public double alpha;
	
	/**
	 * Number of random trials to estimate average energy difference.
	 */
	public int nTrial;  
	
	
	/**
	 * Sets all fields.
	 */
	public SimAnnealRunner (double pInit, int nTrans, 
			int nAccept, int nT, int nUnchanged, double alpha, int nTrial) {
		this.pInit = pInit;
		this.nTrans = nTrans;
		this.nAccept = nAccept;
		this.nT = nT;
		this.nUnchanged = nUnchanged;
		this.alpha = alpha;
		this.nTrial = nTrial;
	}

	
	/**
	 * Runs a simulated annealing procedure for combinatorial
	 * minimization, starting with some initial configuration.
	 * 
	 * @param initConfig initial configuration.
	 * @return minimizing configuration.
	 * 
	 * @author tomic
	 */
	public SimAnnealConfiguration run(SimAnnealConfiguration initConfig) {

        // initial temperature
        double Tinit = - initialDeltaE(initConfig, this.nTrial)/Math.log(pInit);
  
        double bestEnergy = initConfig.energyFunction();
        double currentEnergy = bestEnergy;
        double tempEnergy = bestEnergy;
        double testEnergy;
        
        SimAnnealConfiguration currentConfig = initConfig;
 
        try {
        	Date today = new Date();
        	String outputFile = String.format("SimAnneal-%tF-%tR.txt", today, today);
            File f = new File(outputFile);
        	FileWriter fw = new FileWriter(f);
        	PrintWriter pw = new PrintWriter(fw);
            
        	// annealing procedure
            int k = 0;
            int cntUnchanged = 0;
            while (k<nT && cntUnchanged < nUnchanged){
        		double T = Math.pow(alpha, k)*Tinit; // decrease temperature
         		int cntTrans = 0;
        		int cntAccept = 0;
        		while (cntTrans < nTrans && cntAccept < nAccept){
        			// make random change in current config
        			SimAnnealConfiguration testConfig = currentConfig.randomChange();
        			testEnergy = testConfig.energyFunction();
        			// monte carlo
        			if (acceptChange(currentEnergy, testEnergy, T)){
        				currentConfig = testConfig;
        				currentEnergy = testEnergy;
        				cntAccept += 1;
        				pw.println(currentEnergy);
        				// if current energy better than best energy, update best energy
        				if (currentEnergy < bestEnergy){
        					bestEnergy = currentEnergy;
        				}
        			}
        			cntTrans += 1;
        		}
        		// check if current energy is different from value at previous T step
        		if (currentEnergy == tempEnergy)
        			cntUnchanged += 1;
        		else
        			tempEnergy = currentEnergy;
        		k += 1;
         	}
        } catch (IOException e) {
            System.out.println("IOException: " + e);
        }
        System.out.println("Best energy: " + bestEnergy);
        return currentConfig;
	}
	
	/**
	 * Generates initial estimate for average energy difference. Used to find
	 * an initial temperature corresponding to {@link SimAnnealRunner#pInit}.
	 * 
	 * @author tomic
	 */
	public double initialDeltaE (SimAnnealConfiguration conf, int nTrial) {
		DoubleArrayList E = new DoubleArrayList(nTrial);
		ArrayList<SimAnnealConfiguration> configs = new ArrayList<SimAnnealConfiguration>(nTrial);
		for (int i=0; i<nTrial; i++){
			// generate nTrial random configs
			conf.randomInit();
			configs.add(conf);
			double en = conf.energyFunction();
			E.add(en);
		}
		// average increase in energy
		double inc = 0;
		if (nTrial != 1)
			inc = (Descriptive.sum(E) - nTrial*Descriptive.min(E))/(nTrial-1);
		return inc;
	}
	
	/**
	 * Monte Carlo. Accept a change with probability 1 if the new energy is lower 
	 * than the original energy, and with probability <em>e<sup>-(E<sub>test</sub>
	 * - E<sub>orig</sub>)/T</sup></em> otherwise.
	 * 
	 * @param orig energy value of original configuration.
	 * @param test test energy value of test configuration.
	 * @param T current temperature
	 * @return <code>true</code> if the test energy is accepted.
	 * 
	 * @author tomic
	 */
	public boolean acceptChange(double orig, double test, double T){
		boolean accept = false;
		if (Math.random() <= Math.exp(-(test-orig)/T) )
			accept = true;
		return accept;
	}
}
