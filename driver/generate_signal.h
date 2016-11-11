struct sig_struct { //defining a structure to hold the generated signal's data
  int inds[];
  double spx[];
  double nu;
};

sig_struct generate_signal(int sigsize, int sparsity, double noise){ //function generate_signal analagous to the Matlab namesake
  sig_struct x; //defining a sig_struct x which we will load data into;
  randomSeed(analogRead(0)); //setting the random num generators seed based on the analogRead of a disconnected pin

  for(int i = 0 ; i < sparsity ; i++){ //iterating 'sparsity' number of times
    x.inds[i] = random(sigsize); //the ith element of inds is loaded by random() with upper bound sigsize
    x.spx[i] = random(sigsize); //the ith eleemnt of spx is loaded by random() with upper bound sigsize
  }
  x.nu = noise; //loading nu to the passed value of noise
  return x; //return the loaded structure
}
