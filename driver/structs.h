#include <StandardCplusplus.h>
#include <vector>

struct sig_struct{
	int len;
	std::vector<int> inds;
 	std::vector<double> spx;
  	double nu;
};

struct tspair{
	int t;
  int s;
};

struct Lambda{ //struct to store frequencies and their corresponding coefficients
  double freq[]; //frequencies
  double coef[]; //coefficients
};