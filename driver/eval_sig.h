#include <vector>

std::vector <Complex> eval_sig(sig_struct x, double pts[], int ptsLen, int N){
	// s has already been initializied to the proper size (ptslen)
	std::vector <Complex> s;
	for(int i = 0 ; i < ptsLen ; i++){
		s.push_back( Complex(0,0) );
	}

	for(int j = 0 ; j < ptsLen ; j++){
		for(int l = 0 ; l < x.len ; l++){
			s[j] += Complex(x.spx[l],0) * (Complex(2,0) * Complex(pi,0) * i * pts[j] * (x.inds[l]-1) / N ).c_exp();
		}
		double randomDecimal = (double) random(1000)/1000; //gaussian distribution number pending
		s[j] = Complex((1/sqrt(N))*(s[j],0) + (x.nu *randomDecimal)); //in MATLAB, randomDecimal would be randn(1)
	}
	return s;
}
