#include <StandardCplusplus.h>
#include <vector>

// sample_residual(r, samples, Lambda, t, sig, N)
//       r: complex[] that the function sets
// samples: sampling points of form t+sj (arithmetic progression)
//  lam: struct containing frequencies and their corresponding coefficients
//       t: double t value
//     sig: double sigma value
//       N: int signal length (power of 2)

std::vector <Complex> sample_residual(std::vector <Complex> samples, Lambda lam, double t, double sig, int N){
	std::vector <Complex> r;
	int k = samples.size();

	for(int i = 0 ; i < k ; i++){ //initializing r as a vector of length k
		r.push_back(Complex(0,0));
	}

	if(sizeof(lam.freq) > 0){
		for (int q = 0 ; q < k ; q++){
			Complex vq(0,0);
			for(int j = 0 ; j < (lam.freq).size() ; j++){
				vq += Complex(lam.coef[j],0) * (Complex(2.0,0)*Complex(PI,0)*i*Complex(lam.freq[q]-1,0)).c_exp() * Complex(sig*(q-1)/N,0);
				// * in the Complex library needs all factors to be of Complex datatype
				// complex.c_exp() is equivalent to exp(complex) -- see Arduino Playground link for more information
			}
			r[q] = samples[q] - Complex(1.0/sqrt(N),0)*vq;
		}
		return r;
	}
	else{
		return samples;
	}
}
