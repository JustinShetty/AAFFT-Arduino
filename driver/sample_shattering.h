// sample_shattering(z, samples, Lambda, t, sig, N)
//       z: complex[] that the function sets
// samples: sampling points of form t+sj (arithmetic progression)
//  Lambda: struct containing frequencies and their corresponding coefficients
//       t: double t value
//     sig: double sigma value
//       N: int signal length (power of 2)
//
//  Result: z is set to an array of the FFTs of residual samples

void sample_shattering(Complex &z, Complex samples[], Lambda lam, double t, double sig, int N){
	Complex temp(0,0);
	sample_residual(temp, samples, lam, t, sig, N);
	double n = sizeof(temp) / sizeof(Complex); 
	z = Complex(1/sqrt(n),0) * Complex(0,0);
	// Complex(0,0) in line 17 should be the FFT of temp!
	// possible FFT implementation:
	// https://github.com/kosme/arduinoFFT/blob/master/Examples/FFT_01/FFT_01.ino
}

//Currently, this function will set the input z to 0+0i ( Complex(0,0) ), but deciding on an FFT implementation to use will resolve this issue
