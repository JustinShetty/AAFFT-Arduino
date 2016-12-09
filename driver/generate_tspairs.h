#include <math.h>

void generate_tspairs(tspair ats1[], tspair ats2[], int N, int reps1, int reps2, int reps3){

	int alpha = log(N)/log(2); //equivalent to log2(N)
  
	for(int i = 0 ; i < reps1*reps2 ; i++){
    tspair temp = {0,0};
    ats1[i] = temp;
	}
	for(int i = 0 ; i < reps1*reps3 ; i++){
		tspair temp = {0,0};
    ats2[i] = temp;
	}


	for(int j =  0 ; j < reps1 ; j++){
    double randomDecimal = (double) random(1000)/1000;
		int r = (int) pow(2,alpha-1)*randomDecimal + 1;
		int s = 2*r - 1; //s is a random odd integer on the interval [1,N-1]
		
		for(int n = 0 ; n < reps2 ; n++){
			double randomDecimal = (double) random(1000)/1000;
			int t =  (int) N*randomDecimal;
      tspair temp = {t,s};
			ats1[(j*reps2)+n] = temp;
		}
		
		for(int n = 0 ; n < reps3 ; n++){
      double randomDecimal = (double) random(1000)/1000;
			int r = (int) pow(2,alpha-1)*randomDecimal + 1;
      int s = 2*r - 1; //s is a random odd integer on the interval [1,N-1]
      
      randomDecimal = (double) random(1000)/1000;
      int t =  (int) N*randomDecimal;
      tspair temp = {0,s};
      
			ats2[(j*reps3)+n] = temp;
		}
	}


	return;
}
