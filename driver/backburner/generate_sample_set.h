void generate_sample_set(std::vector< std::vector< std::vector<Complex> > > &xs1,
						 std::vector< std::vector<Complex> > &xs2,
						 std::vector< std::vector< std::vector<Complex> > > &samp1,
						 std::vector< std::vector<Complex> > &samp2, 
						 sig_struct x, int N, int m, 
						 std::vector <tspair> ats1, std::vector <tspair> ats2, 
						 int width, int input_type){
	// initialization
	xs1.clear();
	xs2.clear();
	samp1.clear();
	samp2.clear();



	int K = width*m;
	// xs1 and samp1
	int nr1 = ats1.size();
	for(int j = 0 ; j < nr1 ; j++){
		int t     = ats1[j].t;
		int sig   = ats1[j].s;
		int final = t + sig*(K-1);

		std::vector<int> aprog;
		for(int q = t ; q <= final ; q+=sig){
			aprog.push_back( q );
		}

		int geo_end = log(N)/log(2); // log2(N)
		for(int b = 0 ; b < geo_end ; b++){
			
			std::vector<int> geoprog;
			for(int q = 0 ; q < aprog.size() ; q++){
				geoprog.push_back( ( aprog[q] + N/((int)pow(2,b)) ) % N);
			}

			if(input_type){
				std::vector<Complex> temp = eval_sig(x, geoprog, N);
				for(int w = 0 ; w < geoprog.size() ; w++){
					xs1[b][ w ][j] = temp[w];
				}
			}
			else{
				//xs2[j] = eval_sig_vect(x, aprog, N);
			}

			for(int w = 0 ; w < geoprog.size() ; w++){
				samp1[b][ w ][j] = geoprog[w];
			}
		}
	}

	// xs2 and samp2
	int nr2 = ats2.size();
	for(int j = 0 ; j < nr2 ; j++){

		int t     = ats2[j].t;
		int sig   = ats2[j].s;
		int final = t + sig*(K-1);

		std::vector<int> aprog;
		for(int q = t ; q <= final ; q+=sig){
			aprog.push_back( q % N );
		}

		if(input_type){
			xs2[j] = eval_sig(x, aprog, N);
		}
		else{
			//xs2[j] = eval_sig_vect(x, aprog, N);
		}

		for(int q = 0 ; q < samp2[j].size() ; q++){
			samp2[j][q] = Complex(aprog[q], 0);
		}
	}

}