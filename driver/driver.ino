#include <math.h>

#include "structs.h" //including self-defined data structures
#include "generate_tspairs.h"
//#include "generate_signal.h" //including sig_struct and the function generate_signal
//#include "sample_residual.h" //including Lambda struct and the function sample_residual
//#include "sample_shattering.h" //including the function sample_shattering

using namespace std;

void setup(){
  Serial.begin(9600);
}

void loop(){
  randomSeed(analogRead(0));
  int N = pow(2,10);
  int m = 2;
  double nu = 0.0;

  // generate signal
  // set input type

  int reps1 = 3;
  int reps2 = 5;
  int reps3 = 5;
  int width = 5;
  
  //generate the tspairs
  tspair ats1[reps1*reps2];
  tspair ats2[reps1*reps3];
  
  generate_tspairs(ats1, ats2, N, reps1, reps2, reps3);
  
  for(int i = 0 ; i < reps1*reps2 ; i++){
    Serial.print("t: ");
    Serial.print(ats1[i].t);
    Serial.print(" s: ");
    Serial.println(ats1[i].s);
  }

  Serial.print("\n\n");
  
  //generate the sample set

  //do the fourier sampling

  delay(1000);
}
