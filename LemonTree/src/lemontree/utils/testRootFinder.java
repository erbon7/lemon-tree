/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.utils;

public class testRootFinder {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		RealFunction func = new testFunction();
		RootFinder rf = new RootFinder(func, -1.0, 2.0, 1E-16);
		try{
			rf.bisect();
			System.out.println(rf.xRoot);
		} catch (Exception e){
			System.out.println(e);
		}
	}

}
