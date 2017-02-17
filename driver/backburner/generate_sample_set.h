void generate_sample_set(double &xs1, double &xs2, double &samp1, double &samp2, sig_struct x, int N, int m, double ats1[][], double ats2[][], int width, int input_type ){
	int K = width * m;
	int nr1 = sizeof(ats1)/sizeof(ats1[0]); //number of rows in ats1
	double temp[(log(N)/log(2))+1][K][nr1];
	samp1 = temp;
	double temp[(log(N)/log(2))+1][K][nr1];
	xs1 = temp;
	// TODO: ATS1 is complicated!

	int nr2 = sizeof(ats2)/sizeof(ats2[0]); //number of rows in ats2
	double temp[nr2][K];
	samp2 = temp;
	double temp[nr2][K];
	xs2 = temp;
	for(int j = 0 ; j < nr2 ; j++){
		int t = ats2[j].t;
		int s = ats2[j].s;
		//int final = t + s*(K-1);
		//aprog = t:s:final
		if(input_type == 1){
			//xs2[j] = eval_sig(x, mod(aprog,N), N);
		}
		else{
			//xs2[j] = eval_sig_vect(x, mod(aprog,N), N);
		}
		//samp2[j] = mod(aprog,N);
	}
}