/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.utils;

import BiNGO.BingoParameters;
import BiNGO.BingoAlgorithm;
import BiNGO.AnnotationParser;
import BiNGO.CalculateTestTask;
import BiNGO.CalculateCorrectionTask;

import java.util.Collections;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import cytoscape.data.annotation.OntologyTerm;


/**
 * Wrapper around the BiNGO library to generate GO statistics for each module of a 
 * ModuleNetwork object. 
 * 
 * @author stmae, erbon
 * 
 */
@SuppressWarnings("unchecked")
public class BiNGO {
	
	/**
	 * BingoParameters object.
	 */
    private BingoParameters bp;
    
    /**
     * BiNGO object constructor
     * 
     * @param annotationFile annotation file name (gene_association format) 
     * @param ontologyFile ontology file name (.obo file)
     * @param significance p-value threshold
     * @param namespace as defined in the .obo file, e.g. biological_process
     * @return BiNGO object
     * 
     * @author eric
     */
    public BiNGO(String annotationFile, String ontologyFile, String significance, String nameSpace) {

    	// create bingoparameters object
    	try {
    		this.bp = new BingoParameters("");
    	} catch(IOException e){
    		System.out.println("Error: Can't create the BingoParameters object.");
    		System.exit(1);
    	}

    	// Override standards parameters with 
    	bp.setReferenceSet(BingoAlgorithm.GRAPH);
    	bp.setCategory(BingoAlgorithm.CATEGORY_CORRECTION);
    	bp.setOverOrUnder("Overrepresentation") ;
    	bp.setSignificance(new BigDecimal(significance)) ;                
    	bp.setAnnotationFile(annotationFile);
    	bp.setSpecies(bp.getSpeciesNameFromFilename(annotationFile));
    	bp.setNameSpace(nameSpace);
    	bp.setOntologyFile(ontologyFile) ;
    	bp.setAnnotation_default(false) ;
    	bp.setOntology_default(false);
    	
    	System.out.print("Initialization of the BiNGO parser...");
    	AnnotationParser annParser = bp.initializeAnnotationParser();
    	annParser.calculate();
    	
    	if (annParser.getStatus()) {
    		bp.setAnnotation(annParser.getAnnotation());
    		bp.setOntology(annParser.getOntology());
    		bp.setAlias(annParser.getAlias());
    	}
    	System.out.println(" [ok]");
    	
    	if (bp.getAnnParser().getOrphans())
    		System.out.println("Warning: Some category labels in the annotation file are not defined in the ontology.");
    }

    /**
     * Load cluster file and calculate GO over-representation for each cluster.
     * Re-map gene names if a map file is given.
     * 
     * @param clusterFileName cluster file name.
     * @param referenceSetFileName gene reference set file name
     * @param outputFileName output file name
     * @param mapFile gene map file name
     * @throws IOException
     * 
     * @author eric
     */
    public void GOstats(String clusterFileName, String referenceSetFileName, String outputFileName, String mapFile) throws IOException{

    	
    	// if mapFile is defined, re-map gene names
    	HashMap<String, String> geneMap = new HashMap<String, String>();
    	if (mapFile != null) {
    		try {
    			System.out.print("Reading map file...");
    			BufferedReader buf = new BufferedReader(new FileReader(mapFile));
    			String line;
    			while ((line = buf.readLine()) != null) {
    				if (!line.startsWith("#")) {
    					String[] tk = line.trim().split("\\s+|\\t");
    					geneMap.put(tk[0], tk[1]);
    				}
    			}
    			buf.close();
    		}
    		catch (FileNotFoundException e) {
    			System.out.println("File "+mapFile+" not found.");
    			System.exit(1);
    		}
    		System.out.println(" [ok]");
    	}
    	
    	// read cluster file 
    	System.out.print("Reading cluster file...");
    	ArrayList<String> lines = new ArrayList<String>();
    	HashSet<String> clusterID = new HashSet<String>();
    	try {
    		BufferedReader buf = new BufferedReader(new FileReader(clusterFileName));
    		String line;
    		while((line = buf.readLine()) != null) {
    			if (!line.startsWith("#")) {
    				lines.add(line.trim());
    				String[] tk = line.trim().split("\\s+|\\t");
    				clusterID.add(tk[1]);
    			}
    		}
    		buf.close();
    	}
    	catch (FileNotFoundException f) {
    		System.out.println("Error: file "+clusterFileName+" not found.");
    		System.exit(1);
    	}

    	HashMap<String,HashSet<String>> clusters = new HashMap<String, HashSet<String>>();
    	ArrayList<Integer> orderedID = new ArrayList<Integer>();
    	for (String id : clusterID) {
    		clusters.put(id, new HashSet<String>());
    		orderedID.add(Integer.parseInt(id));
    	}
    	Collections.sort(orderedID);
    	
    	for (String l : lines) {
    		String[] tk = l.split("\\s+|\\t");
    		String clusterId = tk[1];
    		String geneName = tk[0];
    		if (!geneMap.isEmpty()) {
    			String newName = geneMap.get(geneName);
    			if (newName != null)
    				geneName = newName;
    		}
    		clusters.get(clusterId).add(geneName.toUpperCase());
    	}
    	System.out.println(" [ok]");
    	
    	// read gene reference set file
    	System.out.print("Reading gene reference set file...");
    	HashSet<String> allNodes =  new HashSet<String>();
    	try {
    		BufferedReader buf = new BufferedReader(new FileReader(referenceSetFileName));
    		String line;
    		while((line = buf.readLine()) != null) {
    			allNodes.add(line.trim().toUpperCase());
    		}
    		buf.close();
    	}
    	catch (FileNotFoundException f) {
    		System.out.println("Error: file "+referenceSetFileName+" not found.");
    		System.exit(1);
    	}
    	System.out.println(" [ok]");

    	// create output file
    	PrintWriter out = new PrintWriter(new FileWriter(outputFileName));
    	String newline = System.getProperty("line.separator");

    	bp.setFileoutput(false);
    	bp.setAllNodes(allNodes);
    	bp.setCluster_name("z");
    	
    	out.println("# "+ new java.util.Date().toString());
    	out.println("# Cluster file: "+clusterFileName);
    	out.println("# Reference set file: "+referenceSetFileName);
    	out.println("# Ontology file: "+bp.getOntologyFile());
    	out.println("# Annotation file: "+bp.getAnnotationFile());
    	out.println("# Namespace: "+bp.getNameSpace());
    	out.println("# Statistical test: "+bp.getOverOrUnder());
    	out.println("# Correction method: "+bp.getCorrection());
    	out.println("# Significance level: "+bp.getSignificance());
    	out.println("#");
    	out.println("# prop_clust: proportion of genes having the GO category in the cluster");
    	out.println("# prop_ref: proportion of genes having the GO category in the reference set");
    	out.println("#");
    	out.println("# odds-ratio: prop_clust / prop_ref");
    	out.println("#");
    	out.println("# cluster_id GO_code prop_clust prop_ref odds-ratio p-value corrected_p-value go_category");
    	
    	System.out.print("Calculating p-values...");
    	for (int numID : orderedID) {
    		
    		String id = Integer.toString(numID);
    		
    		bp.setSelectedNodes(clusters.get(id));

    		HashMap testMap = null;
    		HashMap correctionMap = null;

    		BingoAlgorithm algorithm = new BingoAlgorithm(bp);
    		CalculateTestTask test = algorithm.calculate_distribution();
    		test.run();
    		testMap = test.getTestMap();
    		CalculateCorrectionTask correction = algorithm.calculate_corrections(testMap);
    		if (correction != null) {
    			correction.run();
    			correctionMap = correction.getCorrectionMap();
    		}
    		
    		HashMap mapSmallX = test.getMapSmallX();
    		HashMap mapSmallN = test.getMapSmallN();
    		HashMap mapBigX = test.getMapBigX();
    		HashMap mapBigN = test.getMapBigN();

    		// order the GO terms according to p-values
    		String[] ordenedKeySet = ordenKeys(correctionMap);
    		
    		// formatter for decimal formats
    		NumberFormat ft = new DecimalFormat("0.00E0");
    		NumberFormat ft2 = new DecimalFormat("0.00");
    		
    		for(String s: ordenedKeySet) {
    			double correctedPval = Double.parseDouble(correctionMap.get(s).toString()); 
    			if( correctedPval < bp.getSignificance().doubleValue()){
    				double pval = Double.parseDouble(testMap.get(new Integer(s)).toString());
    				OntologyTerm o = bp.getOntology().getTerm(Integer.parseInt(s));
    				int smallX = (Integer) mapSmallX.get(new Integer(s));
    				int bigX = (Integer) mapBigX.get(new Integer(s));
    				int smallN = (Integer) mapSmallN.get(new Integer(s));
    				int bigN = (Integer) mapBigN.get(new Integer(s));
    				double den = (double)smallN/(double)bigN;
    				double oddsRatio=0.0;
    				if (den>0)
    					oddsRatio = ((double)smallX/(double)bigX)/den;
    				out.write(id + "\t");
    				out.write(formatGOcat(s) + "\t");
    				out.write(smallX+"/"+bigX+"\t");
    				out.write(smallN+"/"+bigN+"\t");
    				out.write(ft2.format(oddsRatio)+"\t");
    				out.write(ft.format(pval)+"\t");
    				out.write(ft.format(correctedPval)+"\t");
    				out.write(o.getName()+"\t");
    				out.write(newline);
    			}
    		}
    	}
    	out.close();
    	System.out.println(" [ok]");
    }
    
    /**
     * Order hashmap according to p-values
     * 
     * @param map
     * @return string array
     * 
     * @author stmae
     */
    public String [] ordenKeys(HashMap map) {
        
        String[] keySet = new String[map.keySet().size()];
        int count = 0 ;
        for(Object k: map.keySet()){
            keySet[count] = k.toString();
            count++;
        }
        
        for (int i = 1; i < keySet.length; i++) {
            int j = i;
            // get the first unsorted value ...
            String insert_key = keySet[i];
            double val = Double.parseDouble(map.get(keySet[i]).toString());
            // ... and insert it among the sorted
            while ((j > 0) && (val < Double.parseDouble(map.get(keySet[j-1]).toString()))) {
                keySet[j] = keySet[j - 1];
                j--;
            }
            // reinsert value
            keySet[j] = insert_key;
        }
        return keySet;
    }
    
    /**
     * Format integer to GO category string.
     * 
     * @param id
     * @return
     * 
     * @author eric
     */
    private String formatGOcat (String id) {
    	String ret = "GO:";
    	for (int i=0;i<7-id.length();i++)
    		ret += "0";
    	ret += id;
    	return ret;
    }
}

