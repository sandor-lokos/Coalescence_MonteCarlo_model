#include "includes/model_functions.h"


int main(int argc, char *argv[])
{
  if(argc < 2)
	{
		cerr << endl << argv[0] << " need an argument |-->  N_iter  <--| (number of iteration) !!!" << endl << endl;
		return 1;
	}

  long int N_iter = atoi(argv[1]);

  double m  = 0.14;
  double m1 = 0.03 ;
  double m2 = 0.03 ;
  double dy1 = 1.5;
  double dy2 = 1.5;
  double tau = 7.1;
  double T = .170;
  double sigma = 0.2;
  
  vector<double> par_array = { m, m1, m2, dy1, dy2, tau, T, sigma } ;
  vector<double> func_vars_array = { 0. , 0. } ; // Y, pT
  
  ofstream outfile;

  for(double Y = 0.0 ; Y < 0.10 ; Y++)
  {
    stringstream outputfilename("");
    outputfilename << "test_mom_rap_dist_" << Y << ".txt";
    outfile.open(outputfilename.str().c_str());
    for(double pT = 0. ; pT < 2.; pT += 0.1)
    {
      func_vars_array.at(0) = Y ;
      func_vars_array.at(1) = pT ;
      double result = rapidity_momentum_distribution( func_vars_array , par_array, N_iter ) ;
      outfile << Y << "\t" << pT << "\t" << result << endl;
      cerr << Y << "\t" << pT << "\t" << result << endl;
    }
    outfile.close();
  }


	return 0;
}