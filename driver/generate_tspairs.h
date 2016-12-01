#include <math.h>

void generate_tspairs(double &ats1, double &ats2, int N, int reps1, int reps2, int reps3){
	randomSeed(analogRead(0));
	double temp1[reps1*reps2][2];
	double temp2[reps1*reps3][2];
	int alpha = log(N)/log(2); //equivalent to log2(N)
	for(int i = 0 ; i < reps1*reps2 ; i++){
		temp1[i][0] = 0;
		temp1[i][1] = 0;
	}
	for(int i = 0 ; i < reps1*reps3 ; i++){
		temp2[i][0] = 0;
		temp2[i][1] = 0;
	}

	for(int j =  0 ; j < (reps1-1) ; j++){
		int r = (int) pow(2,alpha-1)*random(0,1) + 1;
		int s = 2*r - 1; //s is a uniformly random integer on the interval [1,N-1]
		
		for(int n = 0 ; n < reps2 ; n++){
			int t = (int) N*random(1);
			temp1[(j*reps2)+n][1] = t;
			temp1[(j*reps2)+n][2] = s;
		}
		
		for(int n = 0 ; n < reps3 ; n++){
			int r = (int) pow(2,alpha-1)*random(0,1) + 1;

			int s = 2*r - 1; //s is a uniformly random integer on the interval [1,N-1]
			int t = (int) N*random(1);

			temp2[(j*reps2)+n][1] = t;
			temp2[(j*reps2)+n][2] = s;
		}
		
	}

	ats1 = temp1;
	ats2 = temp2;
	return;
}