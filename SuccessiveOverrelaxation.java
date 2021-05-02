
public class SuccessiveOverrelaxation {
	static final double THRESHOLD = 0.0000001;
	
	public static void main (String args[]) {
		double[][] testA = {{4, -1, -6, 0},
							{-5, -4, 10, 8},
							{0, 9, 4, -2},
							{1, 0, -7, 5}};
		printArr(multiplyMatrix(testA, new double[] {3, -2, 2, 1}));
		double[] testB = {2, 21, -12, -6};
		printArr(sor(testA, testB, 0.5)[0]);
		int bestIters = 1000000000;
		double bestTime = 1000000000;
		double bestOmegaIters = 0;
		double bestOmegaTime = 0;
		for (int i=1; i<570; i++) {
			double omega = 0.001*i;
			double startTime = System.nanoTime();
			double[][] result = sor(testA, testB, omega);
			double endTime = System.nanoTime();
			 System.out.print("omega = " + omega + " required " + result[1][0] + " iterations and returned ");
			 printArr(result[0]);
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
			// System.out.print(iteration + ": ");
			// printArr(phi);
		}
		return new double[][] {phi, {iteration}};
	}
	
	public static double[] getGuess (int len) {
		double[] guess = new double[len];
		for (int i=0; i<guess.length; i++) {
			guess[i] = 0;
		}
		return guess;
	}
	
	public static double computeResidual (double[][] A, double[] b, double[] phi) {
		double[] product = multiplyMatrix(A, phi);
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
	
	// Returns A*phi
	public static double[] multiplyMatrix (double[][] A, double[] phi) {
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
	
	
	public static void printArr (double[] arr) {
		for (int i=0; i<arr.length; i++) {
			System.out.print(arr[i] + " ");
		}
		System.out.println();
	}
	
	public static void printTwoDArr (double[][] arr) {
		for (int i=0; i<arr.length; i++) {
			for (int j=0; j<arr[0].length; j++) {
				System.out.print(arr[i][j] + " ");
			}
			System.out.println();
		}
	}
}
