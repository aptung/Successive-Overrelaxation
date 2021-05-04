import java.util.Random;

public class SuccessiveOverrelaxation {
	static final Random RAND = new Random();
	static final double THRESHOLD = 0.0000000001;
	
	public static void main (String args[]) {
		double[][] testA = generateRandMat(2, 2);
		printTwoDArr(testA);
		testA = multiplyMatrices(transpose(testA), testA);
		double[] testB = generateRandMat(2, 1)[0];
		int bestIters = 1000000000;
		double bestTime = 1000000000;
		double bestOmegaIters = 0;
		double bestOmegaTime = 0;
		for (int i=1; i<1999; i++) {
			double omega = 0.001*i;
			double startTime = System.nanoTime();
			double[][] result = sor(testA, testB, omega);
			double endTime = System.nanoTime();
			System.out.print("omega = " + omega + " required " + result[1][0] + " iterations and returned ");
			printArr(result[0]);
			
			// Creates a histogram showing the number of iterations required for each omega
			for (int j=0; j<result[1][0]/10; j++) {
				System.out.print('*');
			}
			System.out.println();
			if (result[1][0]<=bestIters) {
				bestIters = (int) result[1][0];
				bestOmegaIters = omega;
			}
			if (endTime-startTime<bestTime) {
				bestTime = endTime - startTime;
				bestOmegaTime = omega;
			}
		}
		System.out.println("Solution: ");
		printArr(sor(testA, testB, bestOmegaIters)[0]);
		System.out.println("Best omega is " + bestOmegaIters + " which takes " + bestIters + " iterations");
		System.out.println("Fastest omega is " + bestOmegaTime + " which takes " + bestTime + " nanoseconds");
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
		System.out.println();
	}
	
}
