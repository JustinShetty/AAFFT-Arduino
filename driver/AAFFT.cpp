/*
	AAFFT.cpp - Ann Arbor Fast Fourier Transform on the Arduino
	Written by Justin Shetty, Oct 2016 - April 2017
	Adapted from MATLAB implementation (Anna C. Gilbert, University of Michigan Dept. of Mathematics)
	URL: https://github.com/JustinShetty/AAFFT-Arduino
	Released into the public domain
*/

#include "AAFFT.h"

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
	for(int b = 0 ; b < vec.size() ; b++ ){
		a += vec[b];
	}
	return a;
}

bool complexComp(Complex c1, Complex c2){ //returns true if c1 < c2
	double a = sqrt( pow(c1.real(),2) + pow(c1.imag(),2) );
	double b = sqrt( pow(c2.real(),2) + pow(c2.imag(),2) );
	return a < b;
}

double abs(Complex c){
	return sqrt( pow(c.real(),2) + pow(c.imag(),2));
}

std::vector <Complex> estimate_coeffs(std::vector < std::vector <Complex> > xs, lam Lambda, std::vector <double> Omega, int k, std::vector <tspair> ats, int N){
	int reps = ats.size();
	int L = Omega.size();

	std::vector <Complex> row(L, Complex(0,0) );
	std::vector < std::vector <Complex> > c(reps, row); //c is a 2d vector with width: L and height: reps
	
	for(int j = 1 ; j <= reps ; j++){ //indexing from 1 or 0?
		int t   = ats[j].t;
		int sig = ats[j].s;
		std::vector <Complex> u = sample_residual( xs[j], Lambda, t, sig, N);

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

	return median(c); // MORE DIFFICULT THAN EXPECTED
}

std::vector <Complex> eval_sig(sig_struct x, std::vector <int> pts, int N){
	std::vector <Complex> s( pts.size() );

	for(int j = 0 ; j < pts.size() ; j++){
		for(int l = 0 ; l < x.inds.size(); l++){
			s[j] += Complex(x.spx[l],0) * (Complex(2,0) * Complex(PI,0) * i * pts[j] * (x.inds[l]-1) / N ).c_exp();
		}
		double randomDecimal = (double) random(1000)/1000; //gaussian distribution number pending
		s[j] = Complex((1/sqrt(N))*(s[j],0) + (x.nu *randomDecimal)); //in MATLAB, randomDecimal would be randn(1)
	}
	return s;
}

void fourier_sampling(lam &Lambda,
	                  std::vector< std::vector< std::vector<Complex> > > &xs1,
					  std::vector< std::vector<Complex> > &xs2,
					  int m, std::vector <tspair> ats1, std::vector <tspair> ats2,
					  int reps1, int reps2, int reps3, int N, int width){
	int k = width*m;
	std::vector<double> Omega(k);
	std::vector<Complex> c(k);
	int list_length = 0;

	for(int j = 0 ; j < reps1 ; j++){
		std::vector<Complex> row1(xs1.size(), Complex(0,0));
		std::vector< std::vector<Complex> > mat1(xs1[0].size(), row1);
		std::vector< std::vector< std::vector<Complex> > > temp1(reps2, mat1);
		for(int x = 0 ; x < xs1.size() ; x++){
			for(int y = 0 ; y < xs1[0].size() ; y++){
				for(int z = reps2*(j-1) ; z <= reps2*(j-1)+reps2 ; z++){
					temp1[x][y][z] = xs1[x][y][z];
				}
			}
		}

		std::vector<tspair> temp2;
		for(int c = reps2*(j-1) ; c <= reps2*(j-1)+reps2 ; j++){
			temp2.push_back(ats1[c]);
		}

		Omega = identify_frequencies(temp1, Lambda, k, temp2, N);

		std::vector<Complex> row2(reps3, Complex(0,0));
		std::vector< std::vector<Complex> > temp3(xs2[0].size(), row2);
		for(int x = 0 ; x < reps3 ; x++){
			for(int y = 0 ; y < xs2[0].size() ; y++){
				temp3[x][y] = xs2[x][y];
			}
		}

		std::vector<tspair> temp4;
		for(int c = reps3*(j-1) ; c <= reps3*(j-1)+reps3 ; j++){
			temp4.push_back(ats2[c]);
		}

		c = estimate_coeffs(temp3, Lambda, Omega, k, temp4, N);

		for(int q = 0 ; q < Omega.size() ; q++){
			if(Lambda.freq.size() > 0){
				int I = 0;
				for(int b = 0 ; b < Lambda.coef.size() ; b++){
					if(Lambda.coef[b] == c[q]){
						I = b;
					}
				}
				if(I > 0){
					Lambda.coef[I] += c[q];
				}
				else{
					list_length++;
					Lambda.freq.push_back(Omega[q]);
					Lambda.coef.push_back(c[q]);
				}
			}
			else{
				list_length++;
				Lambda.freq.push_back(Omega[q]);
				Lambda.coef.push_back(c[q]);
			}
		}

		int p;
		if(k < list_length){
			p = k;
		}
		else{
			p = list_length;
		}
		lam temp = Lambda;
		Lambda.freq.clear();
		Lambda.coef.clear();
		for(int b = 0 ; b < p ; b++){ //retaining the top k frequencies
			Complex max = temp.coef[0];
			int maxInd = 0;
			for(int c = 0 ; c < temp.freq.size() ; c++){
				if( complexComp(max,temp.coef[c]) ){
					max = temp.coef[c];
					maxInd = c;
				}
			}
			Lambda.freq.push_back(temp.freq[maxInd]);
			Lambda.coef.push_back(temp.coef[maxInd]);
			temp.freq.erase(temp.freq.begin() + maxInd);
			temp.coef.erase(temp.coef.begin() + maxInd);
		}
	}
}

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

	std::vector<Complex> row1(log(N)/log(2)+1, Complex(0,0));
	std::vector< std::vector<Complex> > mat1(K, row1);
	std::vector< std::vector< std::vector<Complex> > > cube(ats1.size(), mat1);
	xs1 = cube;

	std::vector<int> row2(log(N)/log(2)+1, 0);
	std::vector< std::vector<int> > mat2(K, row2);
	std::vector< std::vector< std::vector<int> > > cube2(ats1.size(), mat2);
	samp1 = cube2;

	std::vector<Complex> row3(ats2.size(), Complex(0,0));
	std::vector< std::vector<Complex> > mat3(K, row3);
	xs2 = mat3;

	std::vector<int> row4(ats2.size(), 0);
	std::vector< std::vector<int> > mat4(K, row4);
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
  for(int b = 0 ; b < sparsity ; b++){ //iterating 'sparsity' number of times
    double randomDecimal = (double) random(1000)/1000;
    int a = (int) sigsize*randomDecimal; //random integer on [0,sigsize)
    bool contains = false;
    for(int j = 0 ; j < b ; j++){
      if(x.inds[j] == a){
        contains = true; //if a appears in x.inds, contains indicates that it should not be included
      }
    }
    if(!contains){
       x.inds.push_back(a); //if x.inds does not contain a, a is added to x.inds (uniqueness check)
    }

    double s = (double) random(-100,101) / 100.0; //the ith element of spx is loaded by random() with upper bound sigsize
    x.spx.push_back(s);
  }
  x.nu = noise; //loading nu to the passed value of noise
  
  return;
}

void generate_tspairs(std::vector <tspair> &ats1, std::vector <tspair> &ats2, int N, int reps1, int reps2, int reps3){

	int alpha = log(N)/log(2); //equivalent to log2(N)
  
	for(int j =  0 ; j < reps1 ; j++){ 
    double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
		int r = (int) pow(2,alpha-1)*randomDecimal + 1;
		int s = 2*r - 1; //random odd integer on the interval [1,N)
		
		for(int n = 0 ; n < reps2 ; n++){
			double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
			int t =  (int) N*randomDecimal; //integer on [0,N)
      		tspair temp = {t,s};
//          Serial.println(t);
//          Serial.println(s);
//          Serial.println(" ");
//          Serial.println((j*reps2)+n);
//          Serial.println(" ");
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

std::vector <double> identify_frequencies(std::vector< std::vector< std::vector<Complex> > > xs, lam Lambda, int k, std::vector <tspair> ats, int N){
	int reps = ats.size();
	std::vector<double> Omega(k);
	int alpha = log(N)/log(2);
	int sig = ats[0].s;
	
	for(int b = 0 ; b < alpha ; b++){
		std::vector<int> vote(k);
		
		for(int j = 0 ; j < reps ; j++){
			int t = ats[j].t;
			std::vector<Complex> samples;
			for(int f = 0 ; f < xs[1].size() ; f++){
				samples.push_back(xs[1][f][j]);
			}
			std::vector<Complex> u = sample_shattering(samples, Lambda, t, sig, N);
			std::vector<Complex> v = sample_shattering(samples, Lambda, t+(N/pow(2,b+1)), sig, N);

			for(int s = 0 ; s < k ; s++){
				Complex E0 = u[s] + v[s]*(-i * Complex(PI,0) * Complex( Omega[s]/pow(2,b) ,0)).c_exp();
				Complex E1 = u[s] - v[s]*(-i * Complex(PI,0) * Complex( Omega[s]/pow(2,b) ,0)).c_exp();
				if(abs(E1) > abs(E0)){
					vote[s]++;
				}
			}
		}

		for(int s = 0 ; s < k ; s++){
			if(vote[s] > reps/2){
				Omega[s] = Omega[s] + pow(2,b);
			}
		}

	}
	std::vector <double> temp;
	for(int b = 0 ; b < Omega.size() ; b++){
		bool present = false;
		for(int c = 0 ; c < temp.size() ; c++){
			if(Omega[b] == temp[c]){
				present = true;
			}
		}
		if(present){
			temp.push_back(Omega[b]);
		}
	}
	return temp;
}


std::vector <Complex> sample_residual(std::vector <Complex> samples, lam Lambda, double t, double sig, int N){
	int k = samples.size();
	std::vector <Complex> r(k);

	if(sizeof(Lambda.freq) > 0){
		for (int q = 0 ; q < k ; q++){
			Complex vq(0,0);
			for(int j = 0 ; j < (Lambda.freq).size() ; j++){
				vq += Lambda.coef[j] * (Complex(2.0,0) * Complex(PI,0) * i * (Lambda.freq[q] - Complex(1,0) ) * Complex( sig*(q-1)/N , 0) ).c_exp();
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

std::vector <Complex> sample_shattering(std::vector <Complex> samples, lam Lambda, double t, double sig, int N){
	
	std::vector <Complex> z = sample_residual(samples, Lambda, t, sig, N);
	int n = z.size();
  	if(n & (n-1) != 0){ //if n is not a power of 2, the FFT will not work
  		return z;
  	}
	
	PlainFFT FFT = PlainFFT();

	double re[n]; 
	double im[n];
	
	for(int b = 0 ; b < n ; b++){
		re[b] =  z[b].real();
		im[b] =  z[b].imag();
	}

  	FFT.Compute(re,im,n,FFT_FORWARD);
  
	for(int b = 0 ; b < n ; b++){
		z[b] = Complex(1/sqrt(n),0) * Complex(re[b], im[b]) ; //store the fft results
	}

	return z;
}
