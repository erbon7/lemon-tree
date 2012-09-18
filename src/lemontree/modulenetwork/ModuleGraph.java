/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.modulenetwork;

import java.util.*;

/**
 * A ModuleGraph is a directed graph between modules. An edge in the graph implies
 * that at least one of the parents of the target module belongs to the source
 * module.
 * 
 * Methods avoid any calling to Gene or Module fields since we are virtually
 * moving genes around without updating these fields.
 * 
 * @author tomic
 *
 */
public class ModuleGraph {

	/**
	 * The actual graph between modules.
	 */
	public HashMap<Module, HashSet<Module>> collapsedNetwork;

	/**
	 * The module network the graph belongs to
	 */
	public ModuleNetwork moduleNetwork;

	/**
	 * Edges from regulators to modules.
	 */
	public HashMap<Gene, HashSet<Module>> reg2Mod; 

	/**
	 * For each module: which regulators does it (virtually) contain
	 */
	public HashMap<Module, HashSet<Gene>> containsReg; 

	/**
	 * For each regulator: in which module is it (virtually).
	 */
	public HashMap<Gene, Module> inMod; 


	/**
	 * Empty constructor
	 */
	public ModuleGraph(){
		
	}

	/**
	 * Construct module graph for given module network.
	 * 
	 * @param moduleNetwork
	 */
	public ModuleGraph(ModuleNetwork moduleNetwork) {
		this.collapsedNetwork = new HashMap<Module, HashSet<Module>>();
		this.moduleNetwork = moduleNetwork;
		this.reg2Mod = new HashMap<Gene, HashSet<Module>>();
		this.containsReg = new HashMap<Module, HashSet<Gene>>();
		this.inMod = new HashMap<Gene, Module>();

		for (Module mod : moduleNetwork.moduleSet) {
			containsReg.put(mod, new HashSet<Gene>());
			collapsedNetwork.put(mod, new HashSet<Module>());
		}

		for (Gene reg : moduleNetwork.regulatorSet) {
			reg2Mod.put(reg, new HashSet<Module>());
			if (moduleNetwork.inModule.containsKey(reg.name)){
				Module mod = moduleNetwork.moduleSet.get(moduleNetwork.inModule.get(reg.name));
				inMod.put(reg, mod);
				containsReg.get(mod).add(reg);
			}
		}
	}

	/**
	 * Add a regulator (parent) of a module to the network.
	 * 
	 * @param regulator
	 * @param module
	 * 
	 * @author tomic
	 * 
	 */
	public void addRegulator(Gene regulator, Module module) {
		// check if regulator is in the regulatorset, if not, do nothing
		if (moduleNetwork.regulatorSet.contains(regulator)) {
			reg2Mod.get(regulator).add(module);
			if (inMod.containsKey(regulator))
				collapsedNetwork.get(inMod.get(regulator)).add(module);
		}
	}

	/**
	 * Remove a regulator (parent) of a module from the network.
	 * The reverse of {@link ModuleGraph#addRegulator(Gene, Module)}.
	 *
	 * @param regulator
	 * @param module
	 * 
	 * @author tomic
	 * 
	 */
	public void removeRegulator(Gene regulator, Module module) {
		// check if regulator is in the regulatorset, if not, do nothing
		if (moduleNetwork.regulatorSet.contains(regulator)) {
			reg2Mod.get(regulator).remove(module);
			// remove edge in collapsedNetwork, unless they are represented by
			// other reg's
			boolean b = inMod.containsKey(regulator);
			for (Gene reg : containsReg.get(inMod.get(regulator)))
				if (reg2Mod.get(reg).contains(module)) {
					b = false;
				}
			if (b)
				collapsedNetwork.get(inMod.get(regulator)).remove(module);
		}
	}

	/**
	 * Add a regulator gene to a module and update network. 
	 * Always apply {@link ModuleGraph#addGene(Gene, Module)} and
	 * {@link ModuleGraph#removeGene(Gene)} in pairs, 
	 * never one of those by itself.
	 * 
	 * @param gene
	 * @param module
	 * 
	 * @author tomic
	 * 
	 */
	public void addGene(Gene gene, Module module) {
		// check if gene is in the regulatorset, if not, do nothing
		if (moduleNetwork.regulatorSet.contains(gene)) {
			containsReg.get(module).add(gene);
			inMod.put(gene, module);
			for (Module mod : reg2Mod.get(gene)) {
				collapsedNetwork.get(module).add(mod);
			}
		}
	}

	/**
	 * Remove a regulator gene from the network. 
	 * Always apply {@link ModuleGraph#addGene(Gene, Module)} and
	 * {@link ModuleGraph#removeGene(Gene)} in pairs, 
	 * never one of those by itself.
	 * 
	 * @param gene
	 * 
	 * @author tomic
	 * 
	 */
	public void removeGene(Gene gene) {
		// check if gene is in the regulatorset and belongs to a module, if not, do nothing
		if (moduleNetwork.regulatorSet.contains(gene) && inMod.containsKey(gene)) {
			containsReg.get(inMod.get(gene)).remove(gene);
			// remove edges in collapsedNetwork, unless they are represented by
			// other reg's
			for (Module mod : reg2Mod.get(gene)) {
				boolean b = true;
				for (Gene reg : containsReg.get(inMod.get(gene)))
					if (reg2Mod.get(reg).contains(mod)) {
						b = false;
					}
				if (b)
					collapsedNetwork.get(inMod.get(gene)).remove(mod);
			}
			inMod.remove(gene);
		}
	}

	/**
	 * Test the cyclicity of the module graph if a regulator would be added 
	 * to the parents of a module. Changes to the graph are undone at the end.
	 * Cyclicity is tested by searching for a path from the module to itself.
	 * 
	 * @param regulator
	 * @param module
	 * 
	 * @author tomic
	 * 
	 */
	public boolean testCycle(Gene regulator, Module module) {
		// check if regulator is in the regulatorset, if not, return false (no cycle)
		if (moduleNetwork.regulatorSet.contains(regulator)) {
			// if regulator already regulates module, return false (no cycle)
			if (!reg2Mod.get(regulator).contains(module)) {
				addRegulator(regulator, module);
				boolean b = modPath(module, module);
				removeRegulator(regulator, module);
				return b;
			}
			return false;
		}
		return false;
	}

	/**
	 * Test the cyclicity of the module graph if a regulator gene would be moved 
	 * from the gene set of module 1 to the gene set of module 2. Changes to the 
	 * graph are undone at the end. Cyclicity is tested by searching for a path 
	 * from module 2 to itself.
	 * 
	 * @param gene
	 * @param module1
	 * @param module2
	 * 
	 * @author tomic
	 */
	public boolean testCycle(Gene gene, Module module1, Module module2) {
		// check if gene is in the regulatorset, if not, return false (no cycle)
		if (moduleNetwork.regulatorSet.contains(gene)) {
			removeGene(gene);
			addGene(gene, module2);
			boolean b = modPath(module2, module2);
			removeGene(gene);
			addGene(gene, module1);
			return b;
		}
		return false;
	}

	public boolean isCyclic(){
		boolean b = false;
		for (Module mod : this.moduleNetwork.moduleSet){
			b = modPath(mod, mod);
			if (b)
				break;
		}
		return b;
	}
	
	
	/**
	 * Test if there is a path between 2 modules in the module graph.
	 * This method works with a set of modules that is initially the set of
	 * targets of module 1 in the graph. Then it is updated to contain the targets 
	 * of the targets, etc., until the set is empty (no more paths continue beyond 
	 * this point), or until module 2 appears in the set (we found a path between 1 and 2).  
	 * 
	 * @param mod1
	 * @param mod2
	 * @return success or failure to find a path.
	 * 
	 * @author tomic
	 * 
	 */
	public boolean modPath(Module mod1, Module mod2) {
		
		HashSet<Module> reachableModules = new HashSet<Module>();
		
		// initially start with children of mod1
		for (Module mod: collapsedNetwork.get(mod1))
			reachableModules.add(mod);
		
		// loop until we find mod2 or can't continue
		while (!reachableModules.isEmpty() && !reachableModules.contains(mod2)) {
			HashSet<Module> nextReachable = new HashSet<Module>();
			for (Module modf: reachableModules) // current modules
				for (Module modr: collapsedNetwork.get(modf)) { // children of current modules
					nextReachable.add(modr);
				}
			reachableModules = nextReachable;
		}
		
		if (reachableModules.isEmpty()){
			return false;
		}
		else if (reachableModules.contains(mod2)){
			return true;
		}
		else {// java not smart enough to see that this wouldn't happen
			return false;
		}
		
	}
	
	
	/**
	 * Count the number of edges in the graph
	 * 
	 * @author erbon
	 * 
	 */
	public int CountEdges() {
		int ct = 0;
		for (Module mod : collapsedNetwork.keySet()) {
			ct += collapsedNetwork.get(mod).size();
		}
		return (ct);
	}

}
