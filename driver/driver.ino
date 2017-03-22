#include "AAFFT.h"

void setup() {
  Serial.begin(9600);
  delay(1);
  Serial.println("Ready");

  randomSeed(analogRead(0));

  int N = (int) pow(2,8);
  int m = 1;
  double nu = 0.01;
  
  sig_struct x;
  generate_signal(x,N,m,nu);
  
  int input_type = 1;
  int reps1 = 3;
  int reps2 = 3;
  int reps3 = 3;
  int width = 2;
  
  //generate the tspairs
  std::vector <tspair> ats1(reps1*reps2);
  std::vector <tspair> ats2(reps1*reps3); 
  generate_tspairs(ats1, ats2, N, reps1, reps2, reps3);
   
  //generate sample set
  int K = width*m;
  int log2N = log(N)/log(2);
  
  Complex xs1[log2N+1][WIDTH*M][REPS1*REPS3];
  Complex xs2[ats2.size()][WIDTH*M];
  int samp1[log2N+1][WIDTH*M][REPS1*REPS3];
  int samp2[ats2.size()][WIDTH*M];
  
  Complex zeroC(0,0);
  for(int x = 0 ; x < log2N+1 ; x++){
    for(int y = 0 ; y < K ; y++){
      for(int z = 0 ; z < ats1.size() ; z++){
        xs1[x][y][z] = zeroC;
        samp1[x][y][z] = 0;
      }
    }
  }
  for(int x = 0 ; x < ats2.size() ; x++){
    for(int y = 0 ; y < K ; y++){
      xs2[x][y] = zeroC;
      samp2[x][y] = 0;
    }
  }
  
  generate_sample_set(xs1, xs2, samp1, samp2, x, N, m, ats1, ats2, width, input_type);
  
  //do the sampling
  lam Lambda;
  fourier_sampling(Lambda, xs1, xs2, m, ats1, ats2, reps1, reps2, reps3, N, width);
  delay(1);
  getFreeRam();
  // Serial.println("frequency, coefficient");
  //Serial.println(Lambda.freq.size());
  // for(int b = 0 ; b < Lambda.freq.size() ; b++){
  //   Serial.print(Lambda.freq[b].real(),6);
  //   Serial.print(" ");
  //   Serial.print(Lambda.freq[b].imag(),6);
  //   Serial.print("i");
  //   Serial.print(", ");
  //   Serial.println(Lambda.coef[b].real(),6);
  //   Serial.print(" ");
  //   Serial.print(Lambda.coef[b].imag(),6);
  //   Serial.print("i");
  // }
  Serial.println("END"); 
}

void loop() {
  delay(1);
}
