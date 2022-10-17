#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
// #include <gsl/gsl_sf_bessel.h>
#include <vector>

using namespace std;

double CoalescenceDistribution(vector<double> par_array);
double dNi(vector<double> par_array);
double integrand( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_vars_array );
double rapidity_momentum_distribution( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_lims_array , long int N_iter );