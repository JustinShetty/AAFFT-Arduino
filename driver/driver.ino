#include <math.h>
#include "Complex/complex.h" //including Complex math library because Arduino doesn't have native complex number compatability

#include "structs.h" //including self-defined data structures
#include "generate_signal.h" //including sig_struct and the function generate_signal
//#include "sample_residual.h" //including Lambda struct and the function sample_residual
//#include "sample_shattering.h" //including the function sample_shattering

void setup(){
  Serial.begin(9600);
	sig_struct x; //initializing the sig_struct
  x = generate_signal(10, 2, 1.1); //generating a signal with arbitrary parameters for testing
  
  //printing the values returned by generate_signal();
  Serial.print("   Indices: "); //intentionally pre-spaced the word so it will align with Amplitudes
  for(int i = 0 ; i < sizeof(x.inds)/sizeof(int) ; i++){
    Serial.print(x.inds[i]);
    Serial.print(" ");
  }
  Serial.print('\n');
  
  Serial.print("Amplitudes: ");
  for(int i = 0 ; i < sizeof(x.spx)/sizeof(int) ; i++){
    Serial.print(x.spx[i]);
    Serial.print(" ");
  }
  Serial.print('\n');

  Serial.print("Noise: ");
  Serial.println(x.nu);
}

void loop(){ 
}
