#include "AAFFT.h"

void setup() {
  Serial.begin(9600);
  delay(1);
  Serial.println("READY");

  randomSeed(analogRead(0));

  int N = (int) pow(2,8);
  int m = 1;
  double nu = 0;
  
  sig_struct x;
  generate_signal(x,N,m,nu);
  x.inds[0] = 3;
  x.spx[0] = 1.0;
//  Serial.println(x.inds[0]);
//  Serial.println(x.spx[0]);
//  Serial.println(x.nu);
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
  
  Complex czero(0,0);
  for(int xq = 0 ; xq < log2N+1 ; xq++){
    for(int y = 0 ; y < K ; y++){
      for(int z = 0 ; z < ats1.size() ; z++){
        xs1[xq][y][z] = czero;
        samp1[xq][y][z] = 0;
      }
    }
  }
  for(int xq = 0 ; xq < ats2.size() ; xq++){
    for(int y = 0 ; y < K ; y++){
      xs2[xq][y] = czero;
      samp2[xq][y] = 0;
    }
  }
  
  generate_sample_set(xs1, xs2, samp1, samp2, x, N, m, ats1, ats2, width, input_type);
  // for(int xq = 0 ; xq < log2N+1; xq++){
  //   for(int y = 0 ; y < WIDTH*M ; y++){
  //     for(int z = 0 ; z < REPS1*REPS3 ; z++){
  //       Serial.println(xs1[xq][y][z]);
  //     }
  //   }
  // }
  
  //do the sampling
  lam Lambda;
  fourier_sampling(Lambda, xs1, xs2, m, ats1, ats2, reps1, reps2, reps3, N, width);
  delay(1);

//   for(int b = 0 ; b < Lambda.freq.size() ; b++){
//     Serial.print(Lambda.freq[b]);
//     Serial.print(" ");
//     Serial.println(Lambda.coef[b]);
//     Serial.println();
//   }
  Serial.println("END"); 
}

void loop() {
  delay(1);
}
