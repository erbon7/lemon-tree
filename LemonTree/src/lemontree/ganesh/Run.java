/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.ganesh;

import java.util.*;
import java.io.*;

public class Run {
	public static void main(String[] args) {

		// get properties in a hashmap
		HashMap<String, String> hp = new HashMap<String, String>();

		if (args.length > 0)
			getProperties(args[0], hp);
		else
			getProperties("GaneSh.properties", hp);

		// Set data directory
		String dataDir = setProperty("dataDir", hp);

		// Set expression data file name
		String expressiondata = setProperty("expressiondata", hp);

		//Set output filename
		String output = setProperty("output", hp);

		//      Set normal-gamma prior paramaters
		//double lambda0 = Double.parseDouble(setProperty("lambda0", hp));
		//double mu0 = Double.parseDouble(setProperty("mu0", hp));
		//double alpha0 = Double.parseDouble(setProperty("alpha0", hp));
		//double beta0 = Double.parseDouble(setProperty("beta0", hp));

		int numRuns = Integer.parseInt(setProperty("numRuns", hp));
		int numSteps = Integer.parseInt(setProperty("numSteps", hp));
		int burnIn = Integer.parseInt(setProperty("burnIn", hp));
		int sampleStep = Integer.parseInt(setProperty("sampleStep", hp));
		int initClust = Integer.parseInt(setProperty("initClust", hp));

		Boolean oneway = Boolean.parseBoolean(setProperty("oneway", hp));
		Boolean outputMatrix = Boolean.parseBoolean(setProperty("outputMatrix",	hp));

		IO run = new IO(initClust, dataDir, expressiondata, output, numRuns, numSteps,
				burnIn, sampleStep, oneway, outputMatrix);
	}

	public static void getProperties(String filename, HashMap<String, String> h) {
		try {

			BufferedReader file = new BufferedReader(new FileReader(filename));
			String inputline;

			// Read data from file
			while ((inputline = file.readLine()) != null) {
				inputline = inputline.trim();
				if (!inputline.startsWith("#")) {
					Scanner s = new Scanner(inputline);
					s.useDelimiter("=");
					if (s.hasNext()) {
						String key = s.next().trim();
						String value = s.next().trim();
						//System.out.println(key + "<=>" + value);
						h.put(key, value);
					}
					s.close();
				}
			}
			file.close();

		} catch (IOException e) {
			System.out.println("Error: IOexception: " + e);
			System.exit(1);
		}
	}

	public static String setProperty(String name, HashMap<String, String> h) {
		String ret;
		if (h.get(name) == null)
			Die("property " + name + " not found in parameter file.");

		ret = h.get(name);
		return (ret);
	}

	public static void Die(String str) {
		System.out.println("Error: " + str);
		System.exit(1);
	}

}
