#include <StandardCplusplus.h>
#include <vector>
#include <math.h>
#include <complex.h>

const double pi = 3.14159265359;
const Complex i(0,1);

#include "structs.h" //including self-defined data structures
#include "generate_signal.h" //including the function generate_signal
#include "generate_tspairs.h" //including the function generate_tspairs
#include "sample_residual.h"  //including the function sample_residual 
#include "eval_sig.h"
//#include "sample_shattering.h" //including the function sample_shattering

using namespace std;

void setup(){
  Serial.begin(9600);
  randomSeed(analogRead(0));
}

void loop(){
  
  int N = (int) pow(2,15);
  int m = 2;
  double nu = 0.01;

  // generate signal
  sig_struct x;
  generate_signal(x,N,m,nu);
  
//  for(int i = 0 ; i < m ; i++){
//    Serial.print(x.inds[i]);
//    Serial.print(" ");
//    Serial.println(x.spx[i]);
//  }
//  Serial.print("\n\n");
  
  int input_type = 1; //generated signal rather than collected via sensor

  int reps1 = 3;
  int reps2 = 5;
  int reps3 = 11;
  int width = 15;
  
  //generate the tspairs
  tspair ats1[reps1*reps2];
  tspair ats2[reps1*reps3];
  
  generate_tspairs(ats1, ats2, N, reps1, reps2, reps3);
  
  for(int i = 0 ; i < reps1*reps2 ; i++){
    //Serial.print("t: ");
    Serial.println(ats1[i].t);
    //Serial.print(" s: ");
    //Serial.println(ats1[i].s);
  }
  Serial.print("\n\n");
  
  //generate the sample set

  //do the fourier sampling

  delay(1000);
}
