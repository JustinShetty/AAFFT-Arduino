void generate_signal(sig_struct &x, int sigsize, int sparsity, double noise){ //function generate_signal analagous to the Matlab namesake
  for(int i = 0 ; i < sparsity ; i++){ //iterating 'sparsity' number of times
    double randomDecimal = (double) random(1000)/1000;
    int a = (int) sigsize*randomDecimal; //random integer on [0,sigsize)
    bool contains = false;
    for(int j = 0 ; j < i ; j++){
      if(x.inds[i] == a){
        contains = true; //if a appears in x.inds, contains indicates that it should not be included
      }
    }
    if(!contains){
       x.inds.push_back(a); //if x.inds does not contain a, a is added to x.inds (uniqueness check)
    }

    double s = (double) random(-100,101) / 100.0; //the ith element of spx is loaded by random() with upper bound sigsize
    x.spx.push_back(s);
  }
  //x.len = sparsity;
  x.nu = noise; //loading nu to the passed value of noise
  
  return; //return the loaded structure
}
