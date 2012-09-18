package lemontree.test;
import lemontree.modulenetwork.ModuleNetwork;

public class TestReadData {
	public static void main(String[] args) {
		String file = "/Users/eric/mat.dat";
		ModuleNetwork m = new ModuleNetwork();
		m.setNormalGammaPriors(0.1F, 0.0F, 0.1F, 0.1F);
		m.readExpressionMatrix(file,null);
		
		for (int i=0;i<m.numGenes;i++) {
			System.out.print(m.geneSet.get(i).name+" ");
			for (int j=0;j<m.numCond;j++) {
				System.out.print(m.data[i][j]+" ");
			}
			System.out.println();
		}
	}
}
