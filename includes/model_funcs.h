#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;

class model_funcs
{
  public:
    double eta_low = -5. ;
    double eta_high = 5. ;

    double phi_low = 0. ;
    double phi_high = 2.* M_PI ;

    double rT_low = 0.0 ;
    double rT_high = 7. ;

    double pT_low = 0.0 ;
    double pT_high = 2.0 ;

    double y_low = -5. ;
    double y_high = 5. ;

    bool Theta(double x, double limit) { return ( x < limit ); }
    double DiracDelta(double x, double eps)
    {
      if( abs(x) > eps)
        return x;
      else
        return 0;
    }
    // double DiracDelta(double x, double eps) { return exp( -x * x / eps / eps ); }

    double dNi_GKL(vector<double> par_array);
    double CoalescenceDistribution_GKL(vector<double> par_array);
    double integrand_GKL( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_vars_array );
    double GKL_rapidity_momentum_distribution( vector<double> func_vars_array , vector<double> par_array , long int N_iter );

    double dNi_general(vector<double> par_array);
    double CoalescenceDistribution_general(vector<double> par_array);
    double integrand_general( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_vars_array );
    double general_rapidity_momentum_distribution( vector<double> func_vars_array , vector<double> par_array , long int N_iter );
};