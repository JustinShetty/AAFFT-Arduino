std::vector<int> flatten(std::vector< std::vector< std::vector<Complex> > > x);
std::vector<int> flatten(std::vector< std::vector<Complex> > x);

void generate_sample_set(std::vector< std::vector< std::vector<Complex> > > &xs1,
						 std::vector< std::vector<Complex> > &xs2,
						 std::vector< std::vector< std::vector<int> > > &samp1,
						 std::vector< std::vector<int> > &samp2, 
						 sig_struct x, int N, int m, 
						 std::vector <tspair> ats1, std::vector <tspair> ats2, 
						 int width, int input_type){
	int K = width*m;

	// initialization
	xs1.clear();
	xs2.clear();
	samp1.clear();
	samp2.clear();

	std::vector<Complex> row1(K, Complex(0,0));
	std::vector< std::vector<Complex> > mat1(log(N)/log(2), row1);
	std::vector< std::vector< std::vector<Complex> > > cube(ats1.size(), mat1);
	xs1 = cube;

	std::vector<int> row2(K, 0);
	std::vector< std::vector<int> > mat2(log(N)/log(2), row2);
	std::vector< std::vector< std::vector<int> > > cube2(ats1.size(), mat2);
	samp1 = cube2;

	std::vector<Complex> row3(K, Complex(0,0));
	std::vector< std::vector<Complex> > mat3(ats2.size(), row3);
	xs2 = mat3;

	std::vector<int> row4(K, 0);
	std::vector< std::vector<int> > mat4(ats2.size(), row4);
	samp2 = mat4;


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
			samp2[j][q] = aprog[q];
		}
	}

}

std::vector<int> flatten(std::vector< std::vector< std::vector<int> > > x){
	std::vector<int> points;
	for(int j = 0 ; j < x.size() ; j++){
		for(int k = 0 ; k < x[j].size() ; k++){
			for(int l = 0 ; l < x[j][k].size(); l++){
				if( std::find(points.begin(), points.end(), x[j][k][l]) == points.end() ){
					points.push_back(x[j][k][l]);
				}
			}
		}
	}
}

std::vector<int> flatten(std::vector< std::vector<int> > x){
	std::vector<int> points;
	for(int j = 0 ; j < x.size() ; j++){
		for(int k = 0 ; k < x[j].size() ; k++){
			if( std::find(points.begin(), points.end(), x[j][k]) == points.end() ){
					points.push_back(x[j][k]);
				}
		}
	}
}