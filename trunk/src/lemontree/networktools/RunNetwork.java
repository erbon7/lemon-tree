/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.networktools;

public class RunNetwork {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String dir = args[0];
		String file = args[1];
		
		Network net = new Network(dir, file);
		
		System.out.println(net.numEdges + "\t" + net.numNodes);
		
		net.printEdges(args[0]+args[2], 1000);
		
	}

}
