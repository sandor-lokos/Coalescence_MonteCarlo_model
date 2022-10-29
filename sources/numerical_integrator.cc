#include "../includes/model_funcs.h"

int main(int argc, char *argv[])
{
  if(argc < 2)
	{
		cerr << endl << argv[0] << " need an argument |-->  N_iter  <--| (number of iteration) !!!" << endl << endl;
		return 1;
	}

  long int N_iter = atoi(argv[1]);

  double m  = 3.1;
  double m1 = 1.5;
  double m2 = 1.5;
  double dy1 = 1.5;
  double dy2 = 1.5;
  double tau = 7.1;
  double T = .170;
  double sigma = 3.2;
  
  vector<double> par_array = { m, m1, m2, dy1, dy2, tau, T, sigma } ;
  vector<double> func_vars_array = { 0. , 0. } ; // Y, pT
  
  ofstream outfile;

  model_funcs model;

  for(double Y = 0.0 ; Y < 1.0 ; Y++)
  {
    stringstream outputfilename("");
    outputfilename << "test_mom_rap_dist_" << Y << ".txt";
    outfile.open(outputfilename.str().c_str());
    for( double pT = 0.5 ; pT < 8.5 ; pT += 0.5 )
    {
      func_vars_array.at(0) = Y ;
      func_vars_array.at(1) = pT ;
      double general_result = model.general_rapidity_momentum_distribution( func_vars_array , par_array, N_iter ) ;
      double Bjorken_result = model.GKL_rapidity_momentum_distribution( func_vars_array , par_array, N_iter ) ;
      outfile << Y << "\t" << pT << "\t" << general_result / ( 2. * M_PI * pT ) << "\t" << Bjorken_result / ( 2. * M_PI * pT ) << endl;
      cerr    << Y << "\t" << pT << "\t" << general_result / ( 2. * M_PI * pT ) << "\t" << Bjorken_result / ( 2. * M_PI * pT ) << endl;
    }
    outfile.close();
  }

	return 0;
}