/*
	AAFFT.cpp - Ann Arbor Fast Fourier Transform on the Arduino
	Written by Justin Shetty, Oct 2016 - April 2017
	Adapted from MATLAB implementation (Anna C. Gilbert, University of Michigan Dept. of Mathematics)
	URL: https://github.com/JustinShetty/AAFFT-Arduino
	Released into the public domain
*/

#include "AAFFT.h"

std::vector <Complex> median(std::vector < std::vector <Complex> > m){
	std::vector <Complex> res;
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
	return c_abs(c1) < c_abs(c2);
}

double c_abs(Complex c){
	return sqrt( pow(c.real(),2) + pow(c.imag(),2));
}

std::vector <Complex> estimate_coeffs(Complex xs[][WIDTH*M], lam Lambda, std::vector <int> Omega, int k, tspair ats[], int N){
	int reps = REPS3;
	int L = Omega.size();

	std::vector <Complex> row(L, Complex(0,0) );
	std::vector < std::vector <Complex> > c(reps, row); //c is a 2d vector with width: L and height: reps
	
	for(int j = 0 ; j < reps ; j++){
		int t   = ats[j].t;
		int sig = ats[j].s;
		std::vector <Complex> u = sample_residual( xs[j], Lambda, t, sig, N);

		for(int l = 0 ; l < L ; l++){

			std::vector <Complex> tempVec;

			for(int w = 0 ; w < u.size() ; w++){
				Complex tempComp(0,0);
				tempComp = u[w] * ( Complex(-2,0)*Complex(PI,0)*i*Complex(Omega[l],0)
								* Complex(sig,0)* Complex(w,0) * Complex(1/(double)N,0) ).c_exp();
				tempVec.push_back(tempComp);
			}
			c[j][l] = sum(tempVec);
			c[j][l] = Complex(sqrt(N)/k, 0) * c[j][l] * ( Complex(-2,0) * Complex(PI,0) * i * Complex(Omega[l],0) * Complex(t/(double)N,0)).c_exp();
    		// c[j][l] *= Complex(sqrt(N)/k,0);
    		// c[j][l] *= ( Complex(-2,0) * Complex(PI,0) * i * Complex(Omega[l],0) * Complex(t/(double)N,0)).c_exp();
		}
	}
  for(int x = 0 ; x < c.size() ; x++){
    for(int y = 0 ; y < c[0].size(); y++){
      // Serial.println(c[x][y]);
    }
  }
	return median(c); // MORE DIFFICULT THAN EXPECTED
}

std::vector <Complex> eval_sig(sig_struct x, std::vector <int> pts, int N){
	std::vector <Complex> s( pts.size() );

	for(int j = 0 ; j < pts.size() ; j++){
		for(int l = 0 ; l < x.inds.size(); l++){
			s[j] += Complex(x.spx[l],0) * (Complex(2,0) * Complex(PI,0) * i * Complex(pts[j],0) * Complex((x.inds[l]-1) / (double)N,0) ).c_exp();
		}
		double randomDecimal = (double) random(1000)/1000; //gaussian distribution number pending
		s[j] = Complex(1/sqrt((double)N),0) * s[j] + Complex( (x.nu *randomDecimal), 0); //in MATLAB, randomDecimal would be randn(1)
	}
	return s;
}

void fourier_sampling(lam &Lambda, Complex xs1[][WIDTH*M][REPS1*REPS2], Complex xs2[][WIDTH*M],
					            int m, std::vector <tspair> ats1, std::vector <tspair> ats2,
					            int reps1, int reps2, int reps3, int N, int width){
	int k = width*m;
	std::vector<int> Omega(k);
	std::vector<Complex> c(k);
	int list_length = 0;

	 for(int j = 0 ; j < REPS1 ; j++){
//	 for(int j = 0 ; j < 1 ; j++){
		//////////////////////////////////////////////////////////////////////////////////////////////
		Complex temp1[log2N+1][WIDTH*M][REPS2];
		for(int x = 0 ; x < log2N+1 ; x++){
			for(int y = 0 ; y < WIDTH*M ; y++){
				for(int z = 0 ; z < REPS2 ; z++){
					temp1[x][y][z] = xs1[x][y][z+(reps2*j)];
					//Serial.println(temp1[x][y][z]);
				}
			}
		}
    
		tspair temp2[REPS2];
		for(int c = 0 ; c < reps2 ; c++){
			temp2[c] = ats1[(reps2*j)+c];
		}
		Omega = identify_frequencies(temp1, Lambda, k, temp2, N);

   		//////////////////////////////////////////////////////////////////////////////////////////////
		Complex temp3[REPS3*REPS1][WIDTH*M];
		for(int x = 0 ; x < reps3 ; x++){
			for(int y = 0 ; y < WIDTH*M ; y++){
				temp3[x][y] = xs2[x][y];
			}
		}

		tspair temp4[REPS3];
		for(int c = 0 ; c <= reps3 ; c++){
			temp4[c] = ats2[c+(reps3*j)];
		}

		c = estimate_coeffs(temp3, Lambda, Omega, k, temp4, N);

		for(int q = 0 ; q < Omega.size() ; q++){
       		Omega[q] += 1;
//       		Serial.println(Omega[q]);
    	}
    	//////////////////////////////////////////////////////////////////////////////////////////////
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
    for(int b = 0 ; b < Lambda.freq.size() ; b++){
     Serial.print(Lambda.freq[b]);
     Serial.print(" ");
     Serial.println(Lambda.coef[b]);
     Serial.println();
   }
   Serial.println();
		// int p;
		// if(k < list_length){
		// 	p = k;
		// }
		// else{
		// 	p = list_length;
		// }
		// lam temp = Lambda;
		// Lambda.freq.clear();
		// Lambda.coef.clear();
		// for(int b = 0 ; b < p ; b++){ //retaining the top k frequencies
		// 	Complex max = temp.coef[0];
		// 	int maxInd = 0;
		// 	for(int c = 0 ; c < temp.freq.size() ; c++){
		// 		if( complexComp(max,temp.coef[c]) ){
		// 			max = temp.coef[c];
		// 			maxInd = c;
		// 		}
		// 	}
		// 	Lambda.freq.push_back(temp.freq[maxInd]);
		// 	Lambda.coef.push_back(temp.coef[maxInd]);
		// 	temp.freq.erase(temp.freq.begin() + maxInd);
		// 	temp.coef.erase(temp.coef.begin() + maxInd);
		// }
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

void generate_sample_set(Complex xs1[][WIDTH*M][REPS1*REPS2], Complex xs2[][WIDTH*M], int samp1[][WIDTH*M][REPS1*REPS2], int samp2[][WIDTH*M], 
						 sig_struct x, int N, int m, std::vector <tspair> ats1, std::vector <tspair> ats2, 
						 int width, int input_type){
	int K = width*m;

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

		int geo_end = log2N+1;
		for(int b = 0 ; b < geo_end ; b++){
			
			std::vector<int> geoprog;
			for(int q = 0 ; q < aprog.size() ; q++){
				geoprog.push_back( ( aprog[q] + N/round(pow(2,b)) ) % N);
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
      std::vector<Complex> temp = eval_sig(x, aprog, N);
      for(int q = 0 ; q < K ; q++){
        xs2[j][q] = temp[q];
      }
		}
		else{
			//xs2[j] = eval_sig_vect(x, aprog, N);
		}

		for(int q = 0 ; q < K ; q++){
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

std::vector <int> identify_frequencies(Complex xs[][WIDTH*M][REPS2], lam Lambda, int k, tspair ats[], int N){
	// for(int xq = 0 ; xq < log2N+1 ; xq++){
	//     for(int y = 0 ; y < WIDTH*M ; y++){
	//       for(int z = 0 ; z < REPS2*REPS1 ; z++){
	//         Serial.println(xs[xq][y][z]);
	//       }
	//     }
 //  	}
	int reps = REPS2;
	std::vector<int> Omega(k,0);
	int alpha = log(N)/log(2);
	int sig = ats[0].s;
	
	for(int b = 0 ; b < alpha ; b++){
//    for(int b = 0 ; b < 1 ; b++){
		std::vector<int> vote(k,0);
		
		for(int j = 0 ; j < reps ; j++){
			int t = ats[j].t;
			std::vector<Complex> samplesU;
			std::vector<Complex> samplesV;
			for(int f = 0 ; f < WIDTH*M ; f++){
				samplesU.push_back(xs[0][f][j]);
			}
			for(int f = 0 ; f < WIDTH*M ; f++){
				samplesV.push_back(xs[b+1][f][j]);
			}
			std::vector<Complex> u = sample_shattering(samplesU, Lambda, t, sig, N);
			std::vector<Complex> v = sample_shattering(samplesV, Lambda, t+(N/round(pow(2,b+1))), sig, N);

			for(int s = 0 ; s < k ; s++){
//        Serial.println(Omega[s]);
        Complex f = (Complex(-1,0) * i * Complex(PI,0) * Complex( Omega[s]/pow(2,b) ,0)).c_exp();
//        Serial.println(f);
				Complex E0 = u[s] + v[s]*f;
				Complex E1 = u[s] - v[s]*f;
//         Serial.print(u[s]);
//         Serial.print(" ");
//         Serial.println(v[s]);
//         Serial.println();
//				 Serial.print(E0);
//				 Serial.print(" ");
//				 Serial.println(E1);
//				 Serial.println();
				if(c_abs(E1) >= c_abs(E0)){
					vote[s]++;
				}
			}
		}
    for(int s = 0 ; s < k ; s++){
//      Serial.println(vote[s]);
//        Serial.println(Omega[s]);
    }
		for(int s = 0 ; s < k ; s++){
			if(vote[s] > reps/2){
				 Omega[s] += round(pow(2,b));   
			}
		}

	}
	std::vector <int> temp;
	for(int b = 0 ; b < Omega.size() ; b++){
		bool present = false;
		for(int c = 0 ; c < temp.size() ; c++){
			if(Omega[b] == temp[c]){
				present = true;
			}
		}
		if(!present){
			temp.push_back(Omega[b]);
		}
	}
	return temp;
}


std::vector <Complex> sample_residual(Complex samples[], lam Lambda, double t, double sig, int N){
	int k = WIDTH*M;
	std::vector <Complex> r(k);
	// for(int q = 0 ; q < k ; q++){
	// 	Serial.println(samples[q]);
	// }
	if(Lambda.freq.size() > 0){
		for (int q = 0 ; q < k ; q++){
			Complex vq(0,0);
			for(int j = 0 ; j < (Lambda.freq).size() ; j++){
				vq += Lambda.coef[j] * (Complex(2.0,0) * Complex(PI,0) * i * (Lambda.freq[q] - Complex(1,0) ) * Complex(t + (sig*q/(double)N) , 0) ).c_exp();
				// * in the Complex library needs all factors to be of Complex datatype
				// complex.c_exp() is equivalent to exp(complex) -- see Arduino Playground link for more information
			}
      Serial.print("v[q]: ");
      Serial.println((1/sqrt(N))*vq);
      Serial.print("samples[q]: ");
      Serial.println(samples[q]);
			r[q] = samples[q] - Complex(1.0/sqrt(N),0)*vq;
      
		}
		return r;
	}
	else{
		std::vector<Complex> temp;
		for(int q = 0 ; q < k ; q++){
			temp.push_back(samples[q]);
		}
		return temp;
	}
}

std::vector <Complex> sample_shattering(std::vector <Complex> samples, lam Lambda, double t, double sig, int N){
	
	Complex temp[WIDTH*M];
	for(int q = 0 ; q < WIDTH*M ; q++){
		temp[q] = samples[q];
	}
	std::vector <Complex> z = sample_residual(temp, Lambda, t, sig, N);
	// for(int q = 0 ; q < z.size() ; q++){
	// 	Serial.println(z[q]);
	// }
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

  	// for(int q = 0 ; q < n ; q++){
  	// 	Serial.println(Complex(re[q],im[q]));
  	// }
  
	for(int b = 0 ; b < n ; b++){
		z[b] = Complex(1/sqrt(n),0) * Complex(re[b], im[b]) ; //store the fft results
	}
	
	// for(int q = 0 ; q < z.size() ; q++){
	// 	Serial.println(z[q]);
	// }

	return z;
}

int getFreeRam()
{
  extern int __heap_start, *__brkval; 
  int v;

  v = (int) &v - (__brkval == 0 ? (int) &__heap_start : (int) __brkval);

  Serial.print("Free RAM = ");
  Serial.println(v, DEC);

  return v;
}

void printDouble( double val, unsigned int precision){
// prints val with number of decimal places determine by precision
// NOTE: precision is 1 followed by the number of zeros for the desired number of decimial places
// example: printDouble( 3.1415, 100); // prints 3.14 (two decimal places)
   Serial.print (int(val));  //prints the int part
   Serial.print("."); // print the decimal point
   unsigned int frac;
   if(val >= 0)
       frac = (val - int(val)) * precision;
   else
       frac = (int(val)- val ) * precision;
   Serial.println(frac,DEC) ;
} 
