#include <fix_fft.h>
#include <vector>

// sample_shattering(z, samples, Lambda, t, sig, N)
//       z: complex[] that the function sets
// samples: sampling points of form t+sj (arithmetic progression)
//  Lambda: struct containing frequencies and their corresponding coefficients
//       t: double t value
//     sig: double sigma value
//       N: int signal length (power of 2)
//
//  Result: z is set to an array of the FFTs of residual samples

void sample_shattering(Complex z[], Complex samples[], Lambda lam, double t, double sig, int N){
  
	Complex temp[];
	sample_residual(temp, samples, lam, t, sig, N);
	double n = sizeof(temp) / sizeof(Complex);
  
  int i = 0;
  int j = 0;
  while(i < n*2){
    fft_input[i] = temp[j].real(); //store data to even indices
    fft_input[i+1] = temp[j].imag();  //blank odd indices
    i+=2;
    j++;
  }
  
  for(int i = 0 ; i < n ; i++){
    z[i] = Complex(1/sqrt(n),0) * fft_lin_out[i]; //store the fft results
  }
}
