#include <math.h> //included for pow()

void generate_tspairs(tspair ats1[], tspair ats2[], int N, int reps1, int reps2, int reps3){

	int alpha = log(N)/log(2); //equivalent to log2(N)
  
	for(int j =  0 ; j < reps1 ; j++){ 
    double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
		int r = (int) pow(2,alpha-1)*randomDecimal + 1;
		int s = 2*r - 1; //random odd integer on the interval [1,N)
		
		for(int n = 0 ; n < reps2 ; n++){
			double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
			int t =  (int) N*randomDecimal; //integer on [0,N)
      tspair temp = {t,s};
			ats1[(j*reps2)+n] = temp;
		}
		
		for(int n = 0 ; n < reps3 ; n++){
      double randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
			int r = (int) pow(2,alpha-1)*randomDecimal + 1;
      int s = 2*r - 1; //random odd integer on the interval [1,N)
      
      randomDecimal = (double) random(1000)/1000; //pseudorandom number on [0,1) to 4 decimals
      int t =  (int) N*randomDecimal;
      tspair temp = {t,s};
      
			ats2[(j*reps3)+n] = temp;
		}
	}
	return;
}
