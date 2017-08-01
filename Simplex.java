import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;


public class Simplex {
	public static double epsilone = 1e-10;
	
	private int basisCount;
	private int nonbasisCount;
	private int[] basis;
	private int[] nonbasis;
	private double[] b;
	private double[][] a;
	private double z;
	private double[] c;
	private int count;
	private double optimum;
	
	void runIntializationSteps() throws IOException{
		Simplex auxCmpx = getAuxSimplex();
		auxCmpx.magicPivoting();
		fixInfeasibleDict(auxCmpx);
		System.out.println("1");
		System.out.println(toString());
	}
	void runOptimizationStep(boolean exceptionOnDegenerate,
			boolean isDual) throws IOException {
		optimum = -1;
		count = 0;
		while(true){
			if(exceptionOnDegenerate & isDegenrate()){
				throw new DegeneratedException();
			}
			int entering = searchEnteringVarIndex();
			if(entering < 0){
				break;
			}
			int leaving = searchLeavingVarIndex(entering);
			if(leaving < 0){
				throw new UnboundedException();
			}
			optimum = fixDictionary(entering, leaving);
			count++;
		}
	}
	
	double fixDictionary(int entering, int leaving){
		double aij = a[leaving][entering];
		b[leaving] = b[leaving]/(-aij);
		
		for(int i = 0; i < nonbasisCount; ++i){
			if(i == entering){
				a[leaving][entering] = (-1)/(-aij);
			}else{
				a[leaving][i] /= (-aij); 
			}
		}
		for(int i = 0; i < basisCount; ++i){
			if(i == leaving) continue;
			double ai_enter = a[i][entering];
			a[i][entering] = 0;
			for(int j = 0; j < nonbasisCount; ++j){
				a[i][j] += ai_enter*a[leaving][j]; 
			}
			b[i] += ai_enter*b[leaving];
		}
		double c_enter =c[entering];
		c[entering] = 0;
		for(int j = 0; j < nonbasisCount; ++j){
			c[j] += c_enter*a[leaving][j]; 
		}
		z += c_enter*b[leaving];
		
		int enterId = nonbasis[entering];
		nonbasis[entering] = basis[leaving];
		basis[leaving] = enterId;
		
		return z;
	}
	int searchEnteringVarIndex(){
		int entering = -1;
		boolean hasEntering = false;
		for(int i = 0; i < nonbasisCount; ++i){
			if(c[i] > epsilone){
				if(!hasEntering){
					entering = i;
					hasEntering = true;
				}else if(nonbasis[entering] > nonbasis[i]){
					entering = i;
				}
			}
		}
		return entering;
		
	}
	int searchLeavingVarIndex(int nonbasis){
		int leaving = -1;
		double minBound = Float.MAX_VALUE;
		for(int i = 0; i < basisCount; ++i){
			double bi = b[i];
			double aij = a[i][nonbasis];
			if(aij < epsilone){
				double bound = bi/(-aij);
				if(minBound > bound){
					minBound= bound;
					leaving = i;
				}else if(minBound == bound && basis[i] < basis[leaving]){
					leaving = i;
				}
			}
		}
		return leaving;
	}
	boolean isDegenrate(){
		for(int i = 0; i < basisCount; ++i){
			if(b[i]  == 0) return true;
		}
		return false;
	}
	boolean isInfeasible(){
		for(int i = 0; i < basisCount; ++i){
			if(b[i]  < 0) return true;
		}
		return false;
	}
	
	List<Integer> getAllFractional(){
		List<Integer> allfrac = new ArrayList<Integer>();
		for(int i = 0; i < basisCount; ++i){
			double floorDiff = b[i] - Math.floor(b[i]);
			double ceilDiff = Math.ceil(b[i]) - b[i];
			//System.out.println(floorDiff + ":" + ceilDiff + ":" + b[i] + ":" + (floorDiff > epsilone)+":" + ( ceilDiff > epsilone ) );
			if(floorDiff > epsilone && ceilDiff > epsilone ) {
				allfrac.add(i);
			}
		}
		return allfrac;
	}
	void fixInfeasibleDict(Simplex auxsmpx){
		System.out.println(auxsmpx.toString());
		if(auxsmpx.optimum == 0){
			z = 0;
			double[] newc = new double[nonbasisCount];
			for(int j = 0; j < nonbasisCount; ++j){
				int nonbasicid = nonbasis[j];
				int auxasisid=-1;
				for(int i = 0; i < auxsmpx.basisCount; ++i){
					//System.out.println(auxsmpx.basis[i] + "==" + nonbasicid);
					if(auxsmpx.basis[i] == nonbasicid) {
						auxasisid = i;
						break;
					}
				}
				if(auxasisid == -1) {
					newc[j] += c[j];
					continue;
				}
				int ithis = 0;
				for(int i = 0; i < auxsmpx.nonbasisCount; ++i){
					if(auxsmpx.nonbasis[i] == 0) continue;
					newc[ithis++] += c[j]*auxsmpx.a[auxasisid][i];
				}
				
				z += c[j]*auxsmpx.b[auxasisid];
			}
			c = newc;
			basisCount = auxsmpx.basisCount;
			nonbasisCount = auxsmpx.nonbasisCount -1;
			basis = new int[basisCount];
			for(int i = 0; i < basisCount; ++i){
				basis[i] = auxsmpx.basis[i];
			}
			b = new double[basisCount];
			for(int i = 0; i < basisCount; ++i){
				b[i] = auxsmpx.b[i];
			}
			
			nonbasis = new int[nonbasisCount];
			int ithis = 0;
			for(int i = 0; i < auxsmpx.nonbasisCount; ++i){
				if(auxsmpx.nonbasis[i] == 0) continue;
				nonbasis[ithis++] = auxsmpx.nonbasis[i];
			}
			
			a = new double[basisCount][nonbasisCount];
			for(int i = 0; i < basisCount; ++i){
				int jthis = 0;
				for(int j = 0; j < auxsmpx.nonbasisCount; ++j){
					if(auxsmpx.nonbasis[j] == 0) continue;
					a[i][jthis++] = auxsmpx.a[i][j];
				}
			}
		}else{
			double x0 = - auxsmpx.optimum;
			for(int i = 0; i < basisCount; ++i){
				b[i] += x0;
			}
		}
	}
	void initialize(String filename) throws FileNotFoundException{
		Scanner sc = new Scanner(new File(filename));
		
		basisCount = sc.nextInt();
		nonbasisCount = sc.nextInt();
		basis = new int[basisCount];
		for(int i = 0; i < basisCount; ++i){
			basis[i] = sc.nextInt();
		}
		nonbasis = new int[nonbasisCount];
		for(int i = 0; i < nonbasisCount; ++i){
			nonbasis[i] = sc.nextInt();
		}
		b = new double[basisCount];
		for(int i = 0; i < basisCount; ++i){
			b[i] = sc.nextFloat();
		}
		a = new double[basisCount][nonbasisCount];
		for(int i = 0; i < basisCount; ++i){
			for(int j = 0; j < nonbasisCount; ++j){
				a[i][j] = sc.nextFloat();
			}
		}
		z = sc.nextFloat();
		c = new double[nonbasisCount];
		for(int i = 0; i < nonbasisCount; ++i){
			c[i] = sc.nextFloat();
		}
	}
	Simplex getAuxSimplex(){
		Simplex smpx = new Simplex();
		
		smpx.nonbasisCount  = nonbasisCount+1;
		smpx.nonbasis = new int[smpx.nonbasisCount];
		smpx.nonbasis[0] = 0;
		for(int i = 1; i < smpx.nonbasisCount; ++i){
			smpx.nonbasis[i] = nonbasis[i-1];
		}
		
		smpx.basisCount = basisCount;
		smpx.basis = new int[smpx.basisCount];
		for(int i = 0; i < smpx.basisCount; ++i){
			smpx.basis[i] = basis[i];
		}
		
		smpx.b = new double[smpx.basisCount];
		for(int i = 0; i < smpx.basisCount; ++i){
			smpx.b[i] = b[i];
		}
		
		smpx.a = new double[smpx.basisCount][smpx.nonbasisCount];
		for(int i = 0; i < smpx.basisCount; ++i){
			smpx.a[i][0] = 1;
			for(int j = 1; j < smpx.nonbasisCount; ++j){
				smpx.a[i][j] = a[i][j-1];
			}
		}
		smpx.z = 0;
		smpx.c = new double[smpx.nonbasisCount];
		smpx.c[0] = -1;
		for(int i = 1; i < smpx.nonbasisCount; ++i){
			smpx.c[i] = 0;
		}
		return smpx;
	}
	void addGOMORY_CHVATAL_CUTS(List<Integer> index){
		int newbasisCount = basisCount + index.size();
		int[] temp = basis;
		basis = new int[newbasisCount];
		System.arraycopy(temp, 0, basis, 0, basisCount);
		
		double[][] tempa = a;
		a = new double[newbasisCount][nonbasisCount];
		for(int i = 0; i < basisCount; ++i){
			System.arraycopy(tempa[i], 0, a[i], 0, nonbasisCount);
		}
		double[] tempb = b;
		b = new double[newbasisCount];
		System.arraycopy(tempb, 0, b, 0, basisCount);
		
		for(int id:index){
			basisCount++;
			basis[basisCount-1] = basisCount + nonbasisCount;
			double[] a_index = a[id];
			for(int i = 0; i < nonbasisCount; ++i){
				a[basisCount-1][i] = (double)((-a_index[i]) - Math.floor(-a_index[i]));
			}
			b[basisCount-1] = -(double)(b[id] - Math.floor(b[id]));
		}
	}
	Simplex getDualSimplex(){
		Simplex smpx = new Simplex();
		
		smpx.basisCount  = nonbasisCount;
		smpx.basis = new int[smpx.basisCount];
		for(int i = 0; i < smpx.basisCount; ++i){
			smpx.basis[i] = nonbasis[i];
		}
		
		smpx.nonbasisCount = basisCount;
		smpx.nonbasis = new int[smpx.nonbasisCount];
		for(int i = 0; i < smpx.nonbasisCount; ++i){
			smpx.nonbasis[i] = basis[i];
		}
		
		smpx.b = new double[smpx.basisCount];
		for(int i = 0; i < smpx.basisCount; ++i){
			smpx.b[i] = -c[i];
		}
		smpx.c = new double[smpx.nonbasisCount];
		for(int i = 0; i < smpx.nonbasisCount; ++i){
			smpx.c[i] = -b[i];
		}
		smpx.a = new double[smpx.basisCount][smpx.nonbasisCount];
		for(int i = 0; i < smpx.basisCount; ++i){
			for(int j = 0; j < smpx.nonbasisCount; ++j){
				smpx.a[i][j] = -a[j][i];
			}
		}
		smpx.z = -z;
		
		return smpx;
	}
	void magicPivoting() {
		int entering = 0;
		int leaving = -1;
		double minBound = Float.MAX_VALUE;
		for(int i = 0; i < basisCount; ++i){
			double bi = b[i];
			if(minBound > bi){
				minBound= bi;
				leaving = i;
			}else if(minBound == bi && basis[i] < basis[leaving]){
				leaving = i;
			}
		}
		count = 1;
		optimum = fixDictionary(entering, leaving);
		while(true){
			entering = searchEnteringVarIndex();
			if(entering < 0){
				break;
			}
			leaving = searchLeavingVarIndex(entering);
			if(leaving < 0){
				break;
			}
			optimum = fixDictionary(entering, leaving);
			count++;
		}
	}
	
	public String toString(){
		StringBuffer sb =  new StringBuffer();
		for(int i = 0; i < basisCount; ++i){
			sb.append("X" + basis[i] + " = " + b[i]);
			for(int j = 0; j < nonbasisCount; ++j){
				sb.append(" + " + a[i][j]  + "*X" + nonbasis[j]);
			}
			sb.append("\n");
		}
		sb.append("\n");
		sb.append("Z = " + z);
		for(int j = 0; j < nonbasisCount; ++j){
			sb.append(" + " + c[j]  + "*X" + nonbasis[j]);
		}
		sb.append("\n");
		
		return sb.toString();
	}
	public int[] getBasis() {
		return basis;
	}
	public int[] getNonbasis() {
		return nonbasis;
	}
	public int getCount() {
		return count;
	}
	public double getOptimum() {
		return optimum;
	}
}


