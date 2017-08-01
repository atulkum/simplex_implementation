import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;


public class RunSimplex {
	public static void main(String[] args) throws IOException{
		
		String exprefix1 = "/Users/atulk/Downloads/course/LP/part1TestCases/unitTests/dict%d.dict";
		String testprefix1 = "/Users/atulk/Downloads/course/LP/part1TestCases/assignmentParts/part%d";
		String exprefix2 = "/Users/atulk/Downloads/course/LP/part2TestCases/unitTests/dict%d";
		String testprefix2 = "/Users/atulk/Downloads/course/LP/part2TestCases/assignmentParts/part%d.dict";
		String exprefix3 = "/Users/atulk/Downloads/course/LP/initializationTests/unitTests/idict%d";
		String testprefix3 = "/Users/atulk/Downloads/course/LP/initializationTests/assignmentTests/part%d.dict";
		String exprefix4 = "/Users/atulk/Downloads/course/LP/ilpTests/assignmentTests/part%d.dict";
		String testprefix4 = "/Users/atulk/Downloads/course/LP/ilpTests/unitTests/ilpTest%d";
		
		runCompleteILPSimplex(4,4, exprefix4, false);
	}
	
	static void runCompleteILPSimplex(int start, int end, String prefix, boolean isWrite) throws IOException{
		for(int i =start;i <= end; i++){
			String filename = String.format(prefix, i);
			Simplex smpx = new Simplex();
			smpx.initialize(filename);
			List<Integer>  fracIndex;
			System.out.println("initial");
			System.out.println(smpx.toString());
			if(smpx.isInfeasible()){
				smpx.runIntializationSteps();
				System.out.println("after AUX initialization");
				System.out.println(smpx.toString());
			}
			smpx.runOptimizationStep(false, false);
			
			for(int k = 0; k < 100; ++k){
			//while(true){
				System.out.println("optimized");
				System.out.println(smpx.toString());
				fracIndex = smpx.getAllFractional();
				if(fracIndex.isEmpty()){
					System.out.println("final");
					System.out.println(smpx.toString());
					break;
				}else{
					smpx.addGOMORY_CHVATAL_CUTS(fracIndex);
					System.out.println("GOMORY_CHVATAL_CUTS");
					System.out.println(smpx.toString());
					Simplex dualsmpx = smpx.getDualSimplex();
					System.out.println("dual");
					System.out.println(dualsmpx.toString());
					try{
						dualsmpx.runOptimizationStep(false, true);
					}catch(UnboundedException ue){
						System.out.println("infeasible");
						break;
					}
					System.out.println("dual optimized");
					System.out.println(dualsmpx.toString());
					smpx = dualsmpx.getDualSimplex();
				}
			}
		}
	}
	static void runCompleteSimplex(int n, String prefix, boolean isWrite) throws IOException{
		for(int i =0;i < n; i++){
			String filename = String.format(prefix,  (i +1));
			Simplex smpx = new Simplex();
			smpx.initialize(filename);
			if(smpx.isInfeasible()){
				smpx.runIntializationSteps();
			}
			smpx.runOptimizationStep(true, false);
			System.out.println(smpx.getOptimum());
			System.out.println(smpx.getCount());
			if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "" + 
					smpx.getOptimum() + "\n"+
					smpx.getCount());
		}
	}
	static void part3(int n, String prefix, boolean isWrite) throws IOException{
		for(int i =0;i < n; i++){
			String filename = String.format(prefix,  (i +1));
			Simplex smpx = new Simplex();
			smpx.initialize(filename);
			//System.out.println(smpx);
			Simplex auxCmpx = smpx.getAuxSimplex();
			auxCmpx.magicPivoting();
			System.out.println(auxCmpx.getOptimum());
			if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "" + 
					smpx.getOptimum() + "\n");
		}
	}
	static void part2(int n, String prefix, boolean isWrite) throws IOException{
		for(int i =0;i < n; i++){
			String filename = String.format(prefix,  (i +1));
			try{
				Simplex smpx = new Simplex();
				smpx.initialize(filename);
				//System.out.println(smpx);
				smpx.runOptimizationStep(true, false);
				
				System.out.println(smpx.getOptimum());
				System.out.println(smpx.getCount());
				if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "" + 
						smpx.getOptimum() + "\n"+
						smpx.getCount());
			}catch(DegeneratedException de){
				System.out.println("DEGENERATE");
				if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "DEGENERATE");
			}catch(UnboundedException de){
				System.out.println("UNBOUNDED");
				if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "UNBOUNDED");
			}
		}
	}
	static void part1(int n, String prefix, boolean isWrite) throws IOException{
		for(int i =0;i < n; i++){
			Simplex smpx = new Simplex();
			String filename = String.format(prefix,  (i +1));
			smpx.initialize(filename);
			//System.out.println(smpx);
			int entering = smpx.searchEnteringVarIndex();
			if(entering < 0){
				System.out.println("UNBOUNDED");
				if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "UNBOUNDED");
				continue;
			}
			int leaving = smpx.searchLeavingVarIndex(entering);
			if(leaving < 0){
				System.out.println("UNBOUNDED");
				if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), "UNBOUNDED");
				continue;
			}
			System.out.println(smpx.getNonbasis()[entering]);
			System.out.println(smpx.getBasis()[leaving]);
			double objective = smpx.fixDictionary(entering, leaving);
			if(isWrite) FileUtils.writeStringToFile(new File(filename + ".output"), 
					"" + smpx.getNonbasis()[entering] + "\n" +
							smpx.getBasis()[leaving] + "\n" +
							objective);
			
			System.out.println(objective);
		}
	}
}
