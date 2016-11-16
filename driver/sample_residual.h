// sample_residual(r, samples, Lambda, t, sig, N)
//       r: double[] that the function sets
// samples: ???????????????????? 
//  Lambda: struct containing frequencies and their corresponding coefficients
//       t: int t value
//     sig: int sigma value

const double pi = 3.14159265359; //what precision do we want for this?
void sample_residual(double &r, double samples[], Lambda lam, double t, double sig, int N){
	if(sizeof(lam.freq) > 0){
		int lamlen = sizeof(lam.freq) / sizeof(double); //length of both arrays in struct Lambda
		double freq[lamlen]; //length dependent on lamlen rather than a literal so the elements
		double coef[lamlen]; //must be set with a for loop rather than at declaration
		for(int i = 0 ; i < lamlen ; i++){
			freq[i] = lam.freq[i];
			coef[i] = lam.coef[i];
		}

		int k = sizeof(samples) / sizeof(double); //k is the length of samples
		double r_temp[k]; //r_temp is an empty array of length k

		for (int q = 0 ; q < k ; q++){
			double vq = 0;
			for(int j = 0 ; j < lamlen ; j++){
				vq += coef[j] * exp(2.0*pi*(freq[q]-1)) * (sig*(q-1)/N); 
				//there should be an i in the above expression but I'm not sure how to handle it
			}
			r_temp[q] = samples[q] - (1.0/sqrt(N))*vq; //subtract off the normalized vq from the existing sample
		}

		r = *r_temp; //setting the value of our return array (include * to derefrence pointer)
	}
	else{
		r = *samples; //if lambda was empty, there are no values to compare samples to, so return samples
	}

	return;
}
