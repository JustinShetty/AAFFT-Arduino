#include "AAFFT.h"

void setup() {
  Serial.begin(9600);
  delay(1);
  Serial.println("Ready");

  randomSeed(analogRead(0));

  int N = (int) pow(2,15);
  int m = 2;
  double nu = 0.01;

  sig_struct x;
  generate_signal(x,N,m,nu);

  int input_type = 1;
  int reps1 = 3;
  int reps2 = 3;
  int reps3 = 3;
  int width = 16;
  
  //generate the tspairs
  std::vector <tspair> ats1(reps1*reps2);
  std::vector <tspair> ats2(reps1*reps3); 
  generate_tspairs(ats1, ats2, N, reps1, reps2, reps3);
  
  Serial.println("wow"); 
  //generate sample set
  std::vector< std::vector< std::vector<Complex> > > xs1;
  std::vector< std::vector<Complex> > xs2;
  std::vector< std::vector< std::vector<int> > > samp1;
  std::vector< std::vector<int> > samp2;
  generate_sample_set(xs1, xs2, samp1, samp2, x, N, m, ats1, ats2, width, input_type);
  
  Serial.println("wow");
  //do the sampling
  lam Lambda;
  fourier_sampling(Lambda, xs1, xs2, m, ats1, ats2, reps1, reps2, reps3, N, width);
  
  Serial.println("frequency, coefficient");
  for(int b = 0 ; b < Lambda.freq.size() ; b++){
    Serial.print(Lambda.freq[b].real(),6);
    Serial.print(" ");
    Serial.print(Lambda.freq[b].imag(),6);
    Serial.print("i");
    Serial.print(", ");
    Serial.println(Lambda.coef[b].real(),6);
    Serial.print(" ");
    Serial.print(Lambda.coef[b].imag(),6);
    Serial.print("i");
  }
}

void loop() {
  delay(1);
}
