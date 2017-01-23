#include "Complex/complex.h" //including Complex math library because Arduino doesn't have native complex number compatability
// http://playground.arduino.cc/Main/ComplexMath

const double pi = 3.14159265359;
const Complex i(0,1);

// estimate_coeffs(c, xs, lam, Omega, k, ats, reps3, N)
//       c: complex[] of coefficients
//      xs: samples
//     lam: freq,coeff pair struct in current approximation
//   Omega: freq,coeff pair struct in this iteration
//       k: length of short filter
//     ats: list of offset,dilation paits for random arrithmetic progreesions
//   reps3: constant set in driver
//       N: int signal length (power of 2)

void estimate_coeffs(Complex &c[], xs, Lambda lam, Lambda Omega, int k, tspair ats[], int reps, int N){
	reps -= 1;
	int L = k; // k equivalent to length(Omega) ?
	Complex c[reps][L];
	
	for(int j = 0 ; j < reps ; j++){
		
		t   = ats[j][0];
		sig = ats[j][1];
		Complex u[];
		sample_residual(u, xs(j:1), lam, t, sig, N); //xs(j:1) ?
		
		for(int l = 0 ; l < L ; l++){
			// c[j][l] = 
			// c[j][l] = 	
		}
	}

	return;
}