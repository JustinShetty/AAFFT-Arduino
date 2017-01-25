#define LIN_OUT 1 //specifies that the linear fft output will be used
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
    fft_input[i] = temp[j]; //store data to even indices
    fft_input[i+1] = 0; //blank odd indices
    i+=2;
    j++;
  }
  //the setup of fft_input[] a little weird but it is how the library authors used it
  fft_window();
  fft_reorder();
  fft_run(); //run the fft operations
  fft_mag_lin(); //store the linear output
  
  for(int i = 0 ; i < n ; i++){
    z[i] = Complex(1/sqrt(n),0) * fft_lin_out[i]; //store the fft results
  }
}
