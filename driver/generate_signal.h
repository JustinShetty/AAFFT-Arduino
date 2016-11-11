struct sig_struct { //defining a structure to hold the generated signal's data
  int inds[]; //indices
  double spx[]; //spectrum of x
  double nu; //noise
};

sig_struct generate_signal(int sigsize, int sparsity, double noise){ //function generate_signal analagous to the Matlab namesake
  sig_struct x; //defining a sig_struct x which we will load data into;
  randomSeed(analogRead(0)); //setting the random num generators seed based on the analogRead of a disconnected pin

  for(int i = 0 ; i < sparsity ; i++){ //iterating 'sparsity' number of times
    double a = (double) random(0,sigsize); //random double on the interval (0,sigsize)
    bool contains;
    for(int j = 0 ; j < i ; j++){
      if(x.inds[i] == a){
        contains = true; //if a appears in x.inds, contains indicates that it should not be included
      }
    }
    if(!contains){
       x.inds[i] = a; //if x.inds does not contain a, x.inds[i] = a (uniqueness check)
    }

    x.spx[i] = (double)random(-1,1); //the ith element of spx is loaded by random() with upper bound sigsize
  }
  x.nu = noise; //loading nu to the passed value of noise
  return x; //return the loaded structure
}
