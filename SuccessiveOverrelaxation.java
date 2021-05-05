import java.util.Random;

public class SuccessiveOverrelaxation {
	static final Random RAND = new Random();
	static final double THRESHOLD = 0.0000000001;
	
//	Some sample times (averages over 20 trials of the times it takes 
//	for the fastest-converging omega to finish)
//	
//	Note that the fastest times are highly variable across different arrays
//	(as are the times for the non-optimal omega)
//	Hence sometimes the average time goes down for larger n
	
//	Also another thing that is noticeable is that the optimal omegas tend to be higher for larger n
// 
//	n = 1: 281.25 nanoseconds (0.28 milliseconds) -- as a benchmark
//	n = 2: 4805.3 nanoseconds (4.8 milliseconds)
//	n = 3: 44527.3 nanoseconds (44.5 milliseconds)
//	n = 4: 59373.55 nanoseconds (59.4 milliseconds)
//	n = 5: 121993.45 nanoseconds (122.0 milliseconds)
//	n = 6: 181313.65 nanoseconds (181.3 milliseconds)
//	n = 7: 132130.35 nanoseconds (132.1 milliseconds)
//	n = 8: 228510.60 nanoseconds (228.5 milliseconds)
//	n = 9: 151463.00 nanoseconds (151.5 milliseconds)
//	n = 10: 195300.00 nanoseconds (195.3 milliseconds)

	public static void main (String args[]) {
		double sum = 0;
		for (int i=1; i<=20; i++) {
			double time = sorTiming(2);
			if (time != 1000000000) {
				sum = sum + time;
				 System.out.println("Current total: " + sum);
				 System.out.println();
			}
			else {
				// System.out.println("Whoops didn't converge");
			}
		}
		System.out.println(sum/20);
	}
	
	// Returns the amount of time (in nanoseconds) the fastest omega took
	// for a random nxn matrix and random nxn vector
	public static double sorTiming (int n) {
		double[][] testA = generateRandMat(n, n);
		// printTwoDArr(testA);
		testA = multiplyMatrices(transpose(testA), testA);
		double[] testB = generateRandMat(n, 1)[0];
		int bestIters = 1000000000;
		double bestTime = 1000000000;
		double bestOmegaIters = 0;
		double bestOmegaTime = 0;
		for (int i=1; i<199; i++) {
			double omega = 0.01*i;
			double startTime = System.nanoTime();
			double[][] result = sor(testA, testB, omega);
			double endTime = System.nanoTime();
			if (result[0][0]!=-1) {
//				System.out.print("omega = " + omega + " required " + result[1][0] + " iterations and returned ");
//				printArr(result[0]);
//				System.out.print("and took " + (endTime-startTime) + " nanoseconds");
//				System.out.println();
				
				if (result[1][0]<=bestIters) {
					bestIters = (int) result[1][0];
					bestOmegaIters = omega;
				}
				if (endTime-startTime<bestTime) {
					bestTime = endTime - startTime;
					bestOmegaTime = omega;
				}
			}
			else {
//				System.out.println("Whoops, " + omega + " became too big");
			}
			
			// Creates a histogram showing the number of iterations required for each omega
//			for (int j=0; j<result[1][0]/10; j++) {
//				System.out.print('*');
//			}
			// System.out.println();
			
		}
//		System.out.println("Solution: ");
//		printArr(sor(testA, testB, bestOmegaIters)[0]);
//		System.out.println("Best omega is " + bestOmegaIters + " which takes " + bestIters + " iterations");
//		System.out.println("Fastest omega is " + bestOmegaTime + " which takes " + bestTime + " nanoseconds");
		return bestTime;
	}
	
	// Performs sor using inputs A, b, omega until the residual is less than THRESHOLD
	// Returns a 2D array composed of the result of sor (1D array) and a 1D array with
	// the number of iterations it took
	public static double[][] sor (double[][] A, double[] b, double omega) {
		double[] phi = getGuess(b.length);
		double residual = computeResidual(A, b, phi);
		int iteration = 0;
		while (residual > THRESHOLD) {
			for (int i=0; i<A[0].length; i++) {
				double sigma = 0;
				for (int j=0; j<A[0].length; j++) {
					if (j!=i) {
						sigma = sigma + A[i][j]*phi[j];
					}
				}
				phi[i] = (1 - omega) * phi[i] + (omega/A[i][i])*(b[i] - sigma);
			}
			if (!isValid(phi)) {
				return new double[][] {{-1}, {-1}};
			}
			residual = computeResidual(A, b, phi);
			iteration++;
		}
		return new double[][] {phi, {iteration}};
	}
	
	// Returns an initial guess (0 vector) for sor
	public static double[] getGuess (int len) {
		double[] guess = new double[len];
		for (int i=0; i<guess.length; i++) {
			guess[i] = 0;
		}
		return guess;
	}
	
	// Computes the residual (a measure of "how good" of a guess phi is)
	// by computing the euclidean distance between A*phi and b (||a*phi-b||)
	public static double computeResidual (double[][] A, double[] b, double[] phi) {
		double[] product = multiplyMatrixVector(A, phi);
		double[] error = new double[phi.length];
		for (int i=0; i<phi.length; i++) {
			error[i] = product[i] - b[i];
		}
		double distance = 0;
		for (int i=0; i<error.length; i++) {
			distance = distance + Math.pow(error[i], 2);
		}
		return Math.pow(distance, 0.5);
	}
	
	// Generates a nxm matrix with random values between 0 and 1
	public static double[][] generateRandMat (int n, int m){
		double[][] randomMat = new double[n][n];
		for (int i=0; i<n; i++) {
			for (int j=0; j<m; j++) {
				randomMat[i][j] = Math.random();
			}
		}
		return randomMat;
	}
	
	// Returns the transpose of A (A^T)
	public static double[][] transpose (double[][] A){
		for (int i=0; i<A.length; i++) {
			for (int j=0; j<A.length; j++) {
				if (j>i) {
					double temp = A[i][j];
					A[i][j] = A[j][i];
					A[j][i] = temp;
				}
			}
		}
		return A;
	}
	
	// Multiplies A*B assuming A and B are both nxn
	public static double[][] multiplyMatrices (double[][] A, double[][] B){
		double[][] result = new double[A.length][A.length];
		for (int i=0; i<A.length; i++) {
			for (int j=0; j<A.length; j++) {
				double sum = 0;
				for (int k=0; k<A.length; k++) {
					sum = sum + A[i][k]*B[k][j];
				}
				result[i][j] = sum;
			}
		}
		return result;
		
	}
	
	// Returns A*phi
	public static double[] multiplyMatrixVector (double[][] A, double[] phi) {
		double[] result = new double[phi.length];
		for (int i=0; i<result.length; i++) {
			double sum = 0;
			for (int j=0; j<phi.length; j++) {
				sum = sum + phi[j]*A[i][j];
			}
			result[i] = sum;
		}
		return result;
	}
	
	public static boolean isValid(double[] arr) {
		for (int j=0; j<arr.length; j++) {
			if (Double.isNaN(arr[j]) || Double.isInfinite(arr[j])) {
				return false;
			}
			if (arr[j]>Math.pow(10, 305)) {
				return false;
			}
		}
		return true;
	}
	
	public static void printTwoDArr (double[][] arr) {
		for (int i=0; i<arr.length; i++) {
			for (int j=0; j<arr[0].length; j++) {
				System.out.print(arr[i][j] + " ");
			}
			System.out.println();
		}
	}
	
	public static void printArr (double[] arr) {
		for (int i=0; i<arr.length; i++) {
			System.out.print(arr[i] + " ");
		}
	}
	
}
