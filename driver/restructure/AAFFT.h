#include "Arduino.h"
#ifndef AAFFT
#define AAFFT

#include "Arduino.h"

#include <StandardCplusplus.h>
#include <vector>
#include <math.h>
#include <complex.h>

struct sig_struct{
	std::vector<int> inds;
	std::vector<double> spx;
	double nu;
};

struct tspair{
	int t;
	int s;
};

struct Lambda{
	std::vector<double> freq;
	std::vector<double> coef;
};

std::vector <Complex> median(std::vector < std::vector <Complex> > m);

Complex sum(std::vector <Complex> vec);

bool complexComp(Complex c1, Complex c2);

std::vector <Complex> estimate_coeffs(std::vector < std::vector <Complex> > xs, 
									  Lambda lam, std::vector <double> Omega, 
									  int k, std::vector <tspair> ats, int N);

std::vector <Complex> eval_sig(sig_struct x, std::vector <int> pts, int N);

// fourier_sampling();

std::vector<int> flatten(std::vector< std::vector< std::vector<Complex> > > x);

std::vector<int> flatten(std::vector< std::vector<Complex> > x);

void generate_sample_set(std::vector< std::vector< std::vector<Complex> > > &xs1,
						 std::vector< std::vector<Complex> > &xs2,
						 std::vector< std::vector< std::vector<int> > > &samp1,
						 std::vector< std::vector<int> > &samp2, 
						 sig_struct x, int N, int m, 
						 std::vector <tspair> ats1, std::vector <tspair> ats2, 
						 int width, int input_type);

void generate_signal(sig_struct &x, int sigsize, int sparsity, double noise);

void generate_tspairs(std::vector <tspair> ats1, std::vector <tspair> ats2, int N, int reps1, int reps2, int reps3);

// identify_frequencies(xs, Lambda, k, ats, N);

std::vector <Complex> sample_residual(std::vector <Complex> samples, Lambda lam, double t, double sig, int N);

std::vector <Complex> sample_shattering(std::vector <Complex> samples, Lambda lam, double t, double sig, int N);



#endif