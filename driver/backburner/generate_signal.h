sig_struct generate_signal(int sigsize, int sparsity, double noise){ //function generate_signal analagous to the Matlab namesake
  sig_struct x; //defining a sig_struct x which we will load data into;
  int temp1[sparsity];
  x.inds = temp1;
  double temp2[sparsity];
  x.spx = temp2;
  //randomSeed(analogRead(0)); //setting the random num generators seed based on the analogRead of a disconnected pin
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

    x.spx[i] = (double)random(-100,101)/100.0; //the ith element of spx is loaded by random() with upper bound sigsize
  }
  x.nu = noise; //loading nu to the passed value of noise
  for(int i = 0 ; i < sparsity ; i++){
    Serial.print(x.spx[i]);
    Serial.print(" ");
  }
  Serial.print("\n");
  return x; //return the loaded structure
}
