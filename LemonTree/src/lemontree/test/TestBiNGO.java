package lemontree.test;

import java.io.*;

import lemontree.utils.BiNGO;

public class TestBiNGO {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String go_annot_file = "/Users/eric/work/GO/gene_association.goa_human";
		//String go_annot_file = "/Users/eric/work/GO/gene_association.mgi";
		String go_ontology_file = "/Users/eric/work/GO/gene_ontology_ext.obo";
		
		String go_p_value ="0.05";
		
		BiNGO b = new BiNGO(go_annot_file, go_ontology_file, go_p_value, "biological_process");
		
		try {
			//b.GOstats("/Users/eric/cluster.txt", "/Users/eric/gene_ref", "/Users/eric/go.txt", null);
			b.GOstats("/Users/eric/work/lemone_test/data/tc4.txt", 
					"/Users/eric/work/lemone_test/go_ref.txt", 
					"/Users/eric/work/lemone_test/go.txt", 
					"/Users/eric/work/lemone_test/data/all_genes_map");
//			b.GOstats("/Users/eric/work/baker_data/baker_modules.txt", 
//					"/Users/eric/work/baker_data/go_ref.txt", 
//					"/Users/eric/work/baker_data/go.txt",
//					"/Users/eric/work/baker_data/all_gene_map");
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

}
