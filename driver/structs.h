struct sig_struct { //struct to store the generated signal's information
  int inds[]; //indices
  double spx[]; //spectrum of x
  double nu; //noise
};

struct Lambda{ //struct to store frequencies and their corresponding coefficients
  double freq[]; //frequencies
  double coef[]; //coefficients
};
