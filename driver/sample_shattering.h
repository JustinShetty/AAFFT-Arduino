#include <fix_fft.h>
#include <vector>

// sample_shattering(samples, Lambda, t, sig, N)
// samples: sampling points of form t+sj (arithmetic progression)
//  Lambda: struct containing frequencies and their corresponding coefficients
//       t: double t value
//     sig: double sigma value
//       N: int signal length (power of 2)
//
//  Result: z is set to an array of the FFTs of residual samples

std::vector <Complex> sample_shattering(std::vector <Complex> samples, Lambda lam, double t, double sig, int N){
	std::vector <Complex> z = sample_residual(samples, lam, t, sig, N);
	int n = z.size();
  
  char re[128]; //samples limited to length 128?
  char im[128];

  for(int i = 0 ; i < n ; i++){
    re[i] = (char) z[i].real();
    im[i]   = (char) z[i].imag();
  }

  fix_fft(re, im, 7, 0); //operations occur in-place
  // re: real component
  // im: imaginary component
  //  m: 0 <= n < 2^m
  //  0: specifies FFT rather than iFFT (1)
  
  for(int i = 0 ; i < n ; i++){
    z[i] = Complex(1/sqrt(n),0) * Complex(re[i], im[i]) ; //store the fft results
  }

  return z;
}
