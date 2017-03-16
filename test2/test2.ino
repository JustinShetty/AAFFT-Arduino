#include "PlainFFT.h"

PlainFFT FFT = PlainFFT(); // Create FFT object
// These values can be changed in order to evaluate the functions
const uint16_t samples = 8;

// These are input and output vectors
double vReal[samples]; 
double vImag[samples];
uint8_t runOnce = 0x00;

void setup(){
	Serial.begin(9600);
	Serial.println("Ready");
}

void loop() {
	if (runOnce == 0x00) {
		runOnce = 0x01;
    if(samples & (samples-1) != 0){
      Serial.println("FUGG");
      return;
    }
		for (uint8_t i = 0; i < samples; i++) {
      vReal[i] = 1.0;
      vImag[i] = 1.0;
      Serial.print(vReal[i],6);
      Serial.print(" ");
      Serial.print(vImag[i],6);
      Serial.println("i");
		}
   
    //FFT.windowing(vReal, samples); 
		FFT.Compute(vReal, vImag, samples, FFT_FORWARD); // Compute FFT
		//FFT.complexToMagnitude(vReal, vImag, samples); // Compute magnitudes
   
   Serial.println(" ");
   for(uint8_t i = 0 ; i < samples ; i++){
     Serial.print(vReal[i],6);
     Serial.print(" ");
     Serial.print(vImag[i],6);
     Serial.println("i");
   }
   Serial.println(" ");
	}
}
