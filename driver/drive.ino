#include "AAFFT.h"

void setup() {
  Serial.begin(9600);
  randomSeed(analogRead(0));

  int N = (int) pow(2,15);
  int m = 2;
  double nu = 0.01;

  sig_struct x;
  generate_signal(x,N,m,nu);

  int input_type = 1;
  int reps1 = 3;
  int reps2 = 5;
  int reps3 = 11;
  int width = 15;
  
  //generate the tspairs
  std::vector <tspair> ats1(reps1*reps2);
  std::vector <tspair> ats2(reps1*reps3); 
  generate_tspairs(ats1, ats2, N, reps1, reps2, reps3);
 
  //generate sample set
  std::vector< std::vector< std::vector<Complex> > > xs1;
  std::vector< std::vector<Complex> > xs2;
  std::vector< std::vector< std::vector<int> > > samp1;
  std::vector< std::vector<int> > samp2;
  generate_sample_set(xs1, xs2, samp1, samp2, x, N, m, ats1, ats2, width, input_type);

  //do the sampling
  //fourier_sampling()
}

void loop() {
  delay(1);
}
