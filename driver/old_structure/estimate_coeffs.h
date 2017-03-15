// estimate_coeffs(c, xs, lam, Omega, k, ats, reps3, N)
//       c: complex[] of coefficients
//      xs: samples
//     lam: freq,coeff pair struct in current approximation
//   Omega: freq,coeff pair struct in this iteration
//       k: length of short filter
//     ats: list of offset,dilation pairs for random arrithmetic progressions
//   reps3: constant set in driver
//       N: int signal length (power of 2)

std::vector <Complex> median(std::vector < std::vector <Complex> > m);
Complex sum(std::vector <Complex> vec);
bool complexComp(Complex c1, Complex c2);

std::vector <Complex> estimate_coeffs(std::vector < std::vector <Complex> > xs, Lambda lam, std::vector <double> Omega, 
										   int k, std::vector <tspair> ats, int N){
	int reps = ats.size();
	int L = Omega.size();

	std::vector <Complex> row(L, Complex(0,0) );
	std::vector < std::vector <Complex> > c(reps, row); //c is a 2d vector with width: L and height: reps
	
	for(int j = 1 ; j <= reps ; j++){ //indexing from 1 or 0?
		int t   = ats[j].t;
		int sig = ats[j].s;
		std::vector <Complex> u = sample_residual( xs[j], lam, t, sig, N);

		for(int l = 1 ; l <= L ; l++){

			std::vector <Complex> tempVec;

			for(int w = 0 ; w < u.size() ; w++){
				Complex tempComp(0,0);
				tempComp = u[w] * ( Complex(-2,0)*Complex(PI,0)*i*Complex(Omega[l],0)
								    *Complex(sig,0)*Complex(1/N,0) ).c_exp();
				tempVec.push_back(tempComp);
			}

			c[j][l] = sum(tempVec);
			c[j][l] = Complex(sqrt(N)/k, 0) * c[j][l] * ( Complex(-2,0) * Complex(PI,0) * i * Complex(Omega[l],0) * (t/N)).c_exp() * Complex(t/N,0); 
		}
	}

	return median(c); // NOT THAT EASY
}

std::vector <Complex> median(std::vector < std::vector <Complex> > m){
	std::vector <Complex> res( m[0].size() );
	Complex val(0,0);
	for(int j = 0 ; j < m.size() ; j++){
		std::sort(m[j].begin(), m[j].end(), complexComp);
		int n = m[j].size();
		if(m.size() % 2 == 0){
			Complex a = m[j][n/2];
			Complex b = m[j][n/2 + 1];
			val = Complex( ( a.real() + b.real() )/2 , ( a.imag() + b.imag() )/2 );
		}
		else{
			val = m[j][n/2];
		}
		res.push_back(val);
	}
	return res;
}

Complex sum(std::vector <Complex> vec){
	Complex a(0,0);
	for(int i = 0 ; i < vec.size() ; i++ ){
		a += vec[i];
	}
	return a;
}

bool complexComp(Complex c1, Complex c2){
	double a = sqrt( pow(c1.real(),2) + pow(c1.imag(),2) );
	double b = sqrt( pow(c2.real(),2) + pow(c2.imag(),2) );
	return a < b;
}