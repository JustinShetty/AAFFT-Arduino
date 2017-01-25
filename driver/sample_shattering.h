#define LIN_OUT 1
#include <FFT.h>

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
	Complex temp(0,0);
	sample_residual(temp, samples, lam, t, sig, N);
	double n = sizeof(temp) / sizeof(Complex);
  int i = 0;
  int j = 0;
  while(i < n*2){
    fft_input[i] = temp[j];
    fft_input[i+1] = 0;
    i+=2;
    j++;
  }
  fft_window();
  fft_reorder();
  fft_run();
  fft_mag_lin();
	z = Complex(1/sqrt(n),0) * Complex(0,0);
  for(int i = 0 ; i < n ; i++){
    z[i] = Complex(1/sqrt(n),0) * fft_lin_out[i];
  }
}
