// estimate_coeffs(c, xs, lam, Omega, k, ats, reps3, N)
//       c: complex[] of coefficients
//      xs: samples
//     lam: freq,coeff pair struct in current approximation
//   Omega: freq,coeff pair struct in this iteration
//       k: length of short filter
//     ats: list of offset,dilation paits for random arrithmetic progreesions
//   reps3: constant set in driver
//       N: int signal length (power of 2)

void estimate_coeffs(Complex &c, Complex xs[], Lambda lam, Lambda Omega, int k, tspair ats[], int reps, int N){
	reps -= 1;
	int L = k; // k equivalent to length(Omega) ?
	Complex temp[reps][L];
	
	for(int j = 0 ; j < reps ; j++){
		
		int t   = ats[j].t;
		int sig = ats[j].s;
		//Complex u[];
		//sample_residual(u, xs(j:1), lam, t, sig, N); //xs(j:1) ?
		
		for(int l = 0 ; l < L ; l++){
			// c[j][l] = 
			// c[j][l] = 	
		}
	}
  for(int i = 0 ; i < reps ; i++){
    for(int j = 0 ; j < L ; i++){
      c[i][j] = temp[i][j];
    }
  }
	return;
}
