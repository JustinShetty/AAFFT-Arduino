#include "Complex/complex.h" //including Complex math library because Arduino doesn't have native complex number compatability
// http://playground.arduino.cc/Main/ComplexMath

// sample_residual(r, samples, Lambda, t, sig, N)
//       r: complex[] that the function sets
// samples: sampling points of form t+sj (arithmetic progression)
//  lam: struct containing frequencies and their corresponding coefficients
//       t: double t value
//     sig: double sigma value
//       N: int signal length (power of 2)

const double pi = 3.14159265359; //what precision do we want for this?
const Complex i(0,1);

void sample_residual(Complex &r, Complex samples[], Lambda lam, double t, double sig, int N){
	if(sizeof(lam.freq) > 0){
		int lamlen = sizeof(lam.freq) / sizeof(double); //length of both arrays in struct Lambda
		double freq[lamlen]; //length dependent on lamlen rather than a literal so the elements
		double coef[lamlen]; //must be set with a for loop rather than at declaration
		for(int i = 0 ; i < lamlen ; i++){
			freq[i] = lam.freq[i];
			coef[i] = lam.coef[i];
		}

		int k = sizeof(samples) / sizeof(double); //k is the length of samples
		Complex r_temp[k]; //r_temp is an empty array of length k

		for (int q = 0 ; q < k ; q++){
			Complex vq(0,0);
			for(int j = 0 ; j < lamlen ; j++){
				vq += Complex(coef[j],0) * (Complex(2.0,0)*Complex(pi,0)*i*Complex(freq[q]-1,0)).c_exp() * Complex(sig*(q-1)/N,0);
				// * in the Complex library needs all factors to be of Complex datatype
				// complex.c_exp() is equivalent to exp(complex) -- see Arduino Playground link for more information
			}
			r_temp[q] = samples[q] - Complex(1.0/sqrt(N),0)*vq; //subtract off the normalized vq from the existing sample
		}

		r = *r_temp; //setting the value of our return array
	}
	else{
		r = *samples; //if lambda was empty, there are no values to compare samples to, so return samples
	}

	return;
}
