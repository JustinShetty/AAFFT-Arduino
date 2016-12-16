const double pi = 3.14159265359;
const Complex i(0,1);

void eval_sig(Complex s[], sig_struct x, double pts[], int ptsLen, int N){
	// s has already been initializied to the proper size (ptslen)
	for(int j = 0 ; j < ptslen ; j++){
		for(int l = 0 ; l < x.len ; l++){
			s[j] += x.spx[l] * exp(2 * pi * i * pts[j] * (x.inds[l]-1) / N );
		}
		double randomDecimal = (double) random(1000)/1000; //gaussian distribution number pending
		s[j] = (1/sqrt(N))*(s[j] + (x.nu *randomDecimal)); //in MATLAB, randomDecimal would be randn(1)
	}
	return;
}