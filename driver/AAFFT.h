#ifndef AAFFT
#define AAFFT

#include "Arduino.h"

#include <StandardCplusplus.h>
#include <vector>
#include <math.h>
#include <complex.h>
#include <PlainFFT.h>

const int N     =  round(pow(2,8));
const int log2N = 8;
const int M     = 1;
const int REPS1 = 3;
const int REPS2 = 3;
const int REPS3 = 3;
const int WIDTH = 2;
const Complex i(0,1);

struct sig_struct{
	std::vector<int> inds;
	std::vector<double> spx;
	double nu;
};

struct tspair{
	int t;
	int s;
};

struct lam{
	std::vector<Complex> freq;
	std::vector<Complex> coef;
};

std::vector <Complex> median(std::vector < std::vector <Complex> > m);

Complex sum(std::vector <Complex> vec);

bool complexComp(Complex c1, Complex c2);

double c_abs(Complex c);

std::vector <Complex> estimate_coeffs(Complex xs[][WIDTH*M], lam Lambda, std::vector <int> Omega, int k, tspair ats[], int N);

std::vector <Complex> eval_sig(sig_struct x, std::vector <int> pts, int N);

 void fourier_sampling(lam &Lambda, Complex xs1[][WIDTH*M][REPS1*REPS2], Complex xs2[][WIDTH*M],
 					  int m, std::vector <tspair> ats1, std::vector <tspair> ats2,
 					  int reps1, int reps2, int reps3, int N, int width);


std::vector<int> flatten(std::vector< std::vector< std::vector<Complex> > > x);

std::vector<int> flatten(std::vector< std::vector<Complex> > x);

void generate_sample_set(Complex xs1[][WIDTH*M][REPS1*REPS2], Complex xs2[][WIDTH*M], int samp1[][WIDTH*M][REPS1*REPS2], int samp2[][WIDTH*M], 
						 sig_struct x, int N, int m, std::vector <tspair> ats1, std::vector <tspair> ats2, 
						 int width, int input_type);

void generate_signal(sig_struct &x, int sigsize, int sparsity, double noise);

void generate_tspairs(std::vector <tspair> &ats1, std::vector <tspair> &ats2, int N, int reps1, int reps2, int reps3);

std::vector <int> identify_frequencies(Complex xs[][WIDTH*M][REPS2], lam Lambda, int k, tspair ats[], int N);

std::vector <Complex> sample_residual(Complex samples[], lam Lambda, double t, double sig, int N);

std::vector <Complex> sample_shattering(std::vector <Complex> samples, lam Lambda, double t, double sig, int N);

int getFreeRam();

void printDouble( double val, unsigned int precision);

#endif
