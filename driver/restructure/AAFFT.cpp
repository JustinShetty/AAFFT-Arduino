#include "AAFFT"

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

std::vector <Complex> estimate_coeffs(std::vector < std::vector <Complex> > xs, Lambda lam, std::vector <double> Omega, int k, std::vector <tspair> ats, int N){
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

std::vector <Complex> eval_sig(sig_struct x, std::vector <int> pts, int N){
	std::vector <Complex> s;
	for(int i = 0 ; i < pts.size() ; i++){
		s.push_back( Complex(0,0) );
	}

	for(int j = 0 ; j < pts.size() ; j++){
		for(int l = 0 ; l < x.inds.size(); l++){
			s[j] += Complex(x.spx[l],0) * (Complex(2,0) * Complex(PI,0) * i * pts[j] * (x.inds[l]-1) / N ).c_exp();
		}
		double randomDecimal = (double) random(1000)/1000; //gaussian distribution number pending
		s[j] = Complex((1/sqrt(N))*(s[j],0) + (x.nu *randomDecimal)); //in MATLAB, randomDecimal would be randn(1)
	}
	return s;
}

// fourier_sampling();

std::vector<int> flatten(std::vector< std::vector< std::vector<int> > > x){
	std::vector<int> points;
	for(int j = 0 ; j < x.size() ; j++){
		for(int k = 0 ; k < x[j].size() ; k++){
			for(int l = 0 ; l < x[j][k].size(); l++){
				if( std::find(points.begin(), points.end(), x[j][k][l]) == points.end() ){
					points.push_back(x[j][k][l]);
				}
			}
		}
	}
}

std::vector<int> flatten(std::vector< std::vector<int> > x){
	std::vector<int> points;
	for(int j = 0 ; j < x.size() ; j++){
		for(int k = 0 ; k < x[j].size() ; k++){
			if( std::find(points.begin(), points.end(), x[j][k]) == points.end() ){
					points.push_back(x[j][k]);
				}
		}
	}
}

void generate_sample_set(std::vector< std::vector< std::vector<Complex> > > &xs1,
						 std::vector< std::vector<Complex> > &xs2,
						 std::vector< std::vector< std::vector<int> > > &samp1,
						 std::vector< std::vector<int> > &samp2, 
						 sig_struct x, int N, int m, 
						 std::vector <tspair> ats1, std::vector <tspair> ats2, 
						 int width, int input_type){
	int K = width*m;

	// initialization
	xs1.clear();
	xs2.clear();
	samp1.clear();
	samp2.clear();

	std::vector<Complex> row1(K, Complex(0,0));
	std::vector< std::vector<Complex> > mat1(log(N)/log(2), row1);
	std::vector< std::vector< std::vector<Complex> > > cube(ats1.size(), mat1);
	xs1 = cube;

	std::vector<int> row2(K, 0);
	std::vector< std::vector<int> > mat2(log(N)/log(2), row2);
	std::vector< std::vector< std::vector<int> > > cube2(ats1.size(), mat2);
	samp1 = cube2;

	std::vector<Complex> row3(K, Complex(0,0));
	std::vector< std::vector<Complex> > mat3(ats2.size(), row3);
	xs2 = mat3;

	std::vector<int> row4(K, 0);
	std::vector< std::vector<int> > mat4(ats2.size(), row4);
	samp2 = mat4;


	// xs1 and samp1
	int nr1 = ats1.size();
	for(int j = 0 ; j < nr1 ; j++){
		int t     = ats1[j].t;
		int sig   = ats1[j].s;
		int final = t + sig*(K-1);

		std::vector<int> aprog;
		for(int q = t ; q <= final ; q+=sig){
			aprog.push_back( q );
		}

		int geo_end = log(N)/log(2); // log2(N)
		for(int b = 0 ; b < geo_end ; b++){
			
			std::vector<int> geoprog;
			for(int q = 0 ; q < aprog.size() ; q++){
				geoprog.push_back( ( aprog[q] + N/((int)pow(2,b)) ) % N);
			}

			if(input_type){
				std::vector<Complex> temp = eval_sig(x, geoprog, N);
				for(int w = 0 ; w < geoprog.size() ; w++){
					xs1[b][ w ][j] = temp[w];
				}
			}
			else{
				//xs2[j] = eval_sig_vect(x, aprog, N);
			}

			for(int w = 0 ; w < geoprog.size() ; w++){
				samp1[b][ w ][j] = geoprog[w];
			}
		}
	}

	// xs2 and samp2
	int nr2 = ats2.size();
	for(int j = 0 ; j < nr2 ; j++){

		int t     = ats2[j].t;
		int sig   = ats2[j].s;
		int final = t + sig*(K-1);

		std::vector<int> aprog;
		for(int q = t ; q <= final ; q+=sig){
			aprog.push_back( q % N );
		}

		if(input_type){
			xs2[j] = eval_sig(x, aprog, N);
		}
		else{
			//xs2[j] = eval_sig_vect(x, aprog, N);
		}

		for(int q = 0 ; q < samp2[j].size() ; q++){
			samp2[j][q] = aprog[q];
		}
	}

}

void generate_signal(sig_struct &x, int sigsize, int sparsity, double noise){ //function generate_signal analagous to the Matlab namesake
  for(int i = 0 ; i < sparsity ; i++){ //iterating 'sparsity' number of times
    double randomDecimal = (double) random(1000)/1000;
    int a = (int) sigsize*randomDecimal; //random integer on [0,sigsize)
    bool contains = false;
    for(int j = 0 ; j < i ; j++){
      if(x.inds[i] == a){
        contains = true; //if a appears in x.inds, contains indicates that it should not be included
      }
    }
    if(!contains){
       x.inds.push_back(a); //if x.inds does not contain a, a is added to x.inds (uniqueness check)
    }

    double s = (double) random(-100,101) / 100.0; //the ith element of spx is loaded by random() with upper bound sigsize
    x.spx.push_back(s);
  }
  //x.len = sparsity;
  x.nu = noise; //loading nu to the passed value of noise
  
  return; //return the loaded structure
}

void generate_tspairs(std::vector <tspair> ats1, std::vector <tspair> ats2, int N, int reps1, int reps2, int reps3){

	int alpha = log(N)/log(2); //equivalent to log2(N)
  
	for(int j =  0 ; j < reps1 ; j++){ 
    double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
		int r = (int) pow(2,alpha-1)*randomDecimal + 1;
		int s = 2*r - 1; //random odd integer on the interval [1,N)
		
		for(int n = 0 ; n < reps2 ; n++){
			double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
			int t =  (int) N*randomDecimal; //integer on [0,N)
      		tspair temp = {t,s};
			ats1[(j*reps2)+n] = temp;
		}
		
		for(int n = 0 ; n < reps3 ; n++){
      		double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
			int r = (int) pow(2,alpha-1)*randomDecimal + 1;
      		int s = 2*r - 1; //random odd integer on the interval [1,N)
      
		    randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
		    int t =  (int) N*randomDecimal;
		    tspair temp = {t,s};
      
			ats2[(j*reps3)+n] = temp;
		}
	}
	return;
}

// identify_frequencies(xs, Lambda, k, ats, N);

std::vector <Complex> sample_residual(std::vector <Complex> samples, Lambda lam, double t, double sig, int N){
	std::vector <Complex> r;
	int k = samples.size();

	for(int i = 0 ; i < k ; i++){ //initializing r as a vector of length k
		r.push_back(Complex(0,0));
	}

	if(sizeof(lam.freq) > 0){
		for (int q = 0 ; q < k ; q++){
			Complex vq(0,0);
			for(int j = 0 ; j < (lam.freq).size() ; j++){
				vq += Complex(lam.coef[j],0) * (Complex(2.0,0)*Complex(PI,0)*i*Complex(lam.freq[q]-1,0)).c_exp() * Complex(sig*(q-1)/N,0);
				// * in the Complex library needs all factors to be of Complex datatype
				// complex.c_exp() is equivalent to exp(complex) -- see Arduino Playground link for more information
			}
			r[q] = samples[q] - Complex(1.0/sqrt(N),0)*vq;
		}
		return r;
	}
	else{
		return samples;
	}
}

std::vector <Complex> sample_shattering(std::vector <Complex> samples, Lambda lam, double t, double sig, int N){
	std::vector <Complex> z = sample_residual(samples, lam, t, sig, N);
	// int n = z.size();
  
  // char re[128]; //samples limited to length 128?
  // char im[128];

  // for(int i = 0 ; i < n ; i++){
  //   re[i] = (char) z[i].real();
  //   im[i]   = (char) z[i].imag();
  // }

  //fix_fft(re, im, 7, 0); //operations occur in-place
  // re: real component
  // im: imaginary component
  //  m: 0 <= n < 2^m
  //  0: specifies FFT rather than iFFT (1)
  
  // for(int i = 0 ; i < n ; i++){
  //   z[i] = Complex(1/sqrt(n),0) * Complex(re[i], im[i]) ; //store the fft results
  // }

  return z;
}