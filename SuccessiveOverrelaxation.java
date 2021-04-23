
public class SuccessiveOverrelaxation {
	static final double THRESHOLD = 0.0000001;
	
	public static void main (String args[]) {
		double[][] testA = {{4, -1, -6, 0}, 
							{-5, -4, 10, 8},
							{0, 9, 4, -2},
							{1, 0, -7, 5}};
		double[] testB = {2, 21, -12, -6};
		int bestIters = 10000;
		double bestOmega = 0;
		for (int i=0; i<1900; i++) {
			double omega = 0.1 + 0.001*i;
			double[][] result = sor(testA, testB, omega);
			System.out.print("omega = " + omega + " required " + result[1][0] + " iterations and returned ");
			printArr(result[0]);
			if (result[1][0]<=bestIters) {
				bestIters = (int) result[1][0];
				bestOmega = omega;
			}
		}
		System.out.println("Best omega is " + bestOmega + " which takes " + bestIters + " iterations");
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
