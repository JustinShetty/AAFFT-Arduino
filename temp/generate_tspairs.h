

void generate_tspairs(double &ats1, double &ats2, int N, int reps1, int reps2, int reps3){
	ats1 = double temp[reps1*reps2][2];
	ats2 = double temp[reps1*reps3][2];
	for(int i = 0 ; i < reps1*reps2 ; i++){
		ats1[i][0] = 0;
		ats1[i][1] = 0;
	}
	for(int i = 0 ; i < reps1*reps3 ; i++){
		ats2[i][0] = 0;
		ats2[i][1] = 0;
	}
}