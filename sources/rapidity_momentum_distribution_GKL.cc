#include "../includes/model_funcs.h"

double model_funcs::CoalescenceDistribution_GKL(vector<double> par_array)
{
  double rT1   = par_array.at(0);
  double rT2   = par_array.at(1);
  double phi1  = par_array.at(2);
  double phi2  = par_array.at(3);
  double Phi1  = par_array.at(4);
  double Phi2  = par_array.at(5);
  double pT1   = par_array.at(6);
  double pT2   = par_array.at(7);
  double m1    = par_array.at(8);
  double m2    = par_array.at(9);
  double tau   = par_array.at(10);
  double T     = par_array.at(11);
  double sigma = par_array.at(12);

  double mT1 = sqrt(m1*m1 + pT1*pT1);
  double mT2 = sqrt(m2*m2 + pT2*pT2);
  double x2_exponent = rT1*rT1 + rT2*rT2 - 2.*(rT1*rT2*cos(phi1-phi2));  // (x1-x2)^2
  double p2_exponent = mT1*mT1 + mT2*mT2 + pT1*pT1 + pT2*pT2 - 2.*pT1*pT2*cos(Phi1-Phi2);  // (p1-p2)^2
  
  // cerr << x2_exponent << " " << sigma*sigma << endl;
  // cerr << p2_exponent << " " << sigma*sigma << endl;
  
  return 9 * M_PI / pow( sigma , 6. ) * Theta( x2_exponent , sigma*sigma ) * Theta( p2_exponent + (m1-m2)*(m1-m2) , sigma * sigma ) ;
  // return 1. / pow( sigma , 6. ) * exp( - x2_exponent / sigma / sigma - p2_exponent * sigma * sigma ) ;
}

double model_funcs::dNi_GKL(vector<double> par_array)
{
  // Implementation of the GKL version
  double rTi  = par_array.at(0);
  double phii = par_array.at(1);
  double Phii = par_array.at(2);
  double pTi  = par_array.at(3);
  double mi   = par_array.at(4);
  double tau  = par_array.at(5);
  double T    = par_array.at(6);
  double dyi  = par_array.at(7);

  double mTi = sqrt(mi*mi + pTi*pTi);
  double Cyi = 1.;//exp( - yi * yi / dyi / dyi );
  double pui = mTi ;
  double gammaTi = 1. / sqrt( 1 - rTi * rTi / tau / tau / 4. ) ;
  
  return Cyi * tau * rTi * pui * exp( - gammaTi / T * ( pui - pTi * rTi / tau /2. * cos( Phii - phii ) ) ) ;
}

double model_funcs::integrand_GKL( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_vars_array )
{
  double Y     = func_vars_array.at(0);
  double pT    = func_vars_array.at(1);

  double m     = par_array.at(0);
  double m1    = par_array.at(1);
  double m2    = par_array.at(2);
  double dy1   = par_array.at(3);
  double dy2   = par_array.at(4);
  double tau   = par_array.at(5);
  double T     = par_array.at(6);
  double sigma = par_array.at(7);

  double Phi   = integ_vars_array.at(0);
  double rT    = integ_vars_array.at(1);
  double pT1   = integ_vars_array.at(2);
  double pT2   = integ_vars_array.at(3);
  double Phi1  = integ_vars_array.at(4);
  double Phi2  = integ_vars_array.at(5);
  double phi1  = integ_vars_array.at(6);
  double phi2  = integ_vars_array.at(7);
  double rT1   = integ_vars_array.at(8);
  double rT2   = integ_vars_array.at(9);

  double mT1 = sqrt(m1*m1 + pT1*pT1);
  double mT2 = sqrt(m2*m2 + pT2*pT2);
  double mT = sqrt(m*m + pT*pT);

  vector<double> coal_pars = {rT1,rT2,phi1,phi2,Phi1,Phi2,pT1,pT2,m1,m2,tau,T,sigma};
  vector<double> momdist1_pars = {rT1,phi1,Phi1,pT1,m1,tau,T,dy1} ;
  vector<double> momdist2_pars = {rT2,phi2,Phi2,pT2,m2,tau,T,dy2} ;

  double result = 10.0; 
  // if ( DiracDelta( pT * cos(Phi) - ( pT1 * cos(Phi1) + pT2 * cos(Phi2) ) , 0.001) > 0. )
    // if ( DiracDelta( pT * sin(Phi) - ( pT1 * sin(Phi1) + pT2 * sin(Phi2) ) , 0.001) > 0. )
      result = dNi_GKL( momdist1_pars ) * dNi_GKL( momdist2_pars ) * CoalescenceDistribution_GKL( coal_pars ) ;
    // if ( (y1 > 0 && y2 > 0) || (y1 < 0 && y2 < 0) )
    // if ( (eta1 > 0 && eta2 > 0) || (eta1 < 0 && eta2 < 0) )
      
  return result ;
}

double model_funcs::GKL_rapidity_momentum_distribution( vector<double> func_vars_array , vector<double> par_array , long int N_iter )
{
  double totalSum = 0 ;
  double functionVal = 0 ;
	int iter = 0;

  vector<double> integ_vars_array = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } ;

	while( iter < N_iter - 1 )
	{
    integ_vars_array.at(0)  = phi_low + ( double( rand() ) / RAND_MAX ) * ( phi_high - phi_low );  // Phi
    integ_vars_array.at(1)  = rT_low +  ( double( rand() ) / RAND_MAX ) * ( rT_high  - rT_low  );  // rT
    integ_vars_array.at(2)  = pT_low +  ( double( rand() ) / RAND_MAX ) * ( pT_high  - pT_low  );  // pT1
    integ_vars_array.at(3)  = pT_low +  ( double( rand() ) / RAND_MAX ) * ( pT_high  - pT_low  );  // pT2
    integ_vars_array.at(4)  = phi_low + ( double( rand() ) / RAND_MAX ) * ( phi_high - phi_low );  // Phi1
    integ_vars_array.at(5)  = phi_low + ( double( rand() ) / RAND_MAX ) * ( phi_high - phi_low );  // Phi2
    integ_vars_array.at(6) = phi_low + ( double( rand() ) / RAND_MAX ) * ( phi_high - phi_low );  // phi1
    integ_vars_array.at(7) = phi_low + ( double( rand() ) / RAND_MAX ) * ( phi_high - phi_low );  // phi2
    integ_vars_array.at(8) = rT_low +  ( double( rand() ) / RAND_MAX ) * ( rT_high  - rT_low  );  // rT1
    integ_vars_array.at(9) = rT_low +  ( double( rand() ) / RAND_MAX ) * ( rT_high  - rT_low  );  // rT2

		functionVal = integrand_GKL( func_vars_array , par_array , integ_vars_array ) ;

    if(iter%10000000==0) cerr << "\t --> Y = " << func_vars_array.at(0) << " , pT = " << func_vars_array.at(1) << " iter = " << iter/1e6 << "M" << endl;

		totalSum += functionVal;

		iter++;
	}

	double estimate = totalSum / N_iter * pow( ( phi_high - phi_low ) , 5 ) * pow( ( rT_high - rT_low ) , 3 ) * pow( ( pT_high - pT_low ) , 2 ) * pow( ( y_high - y_low ) , 2 ) * pow( ( eta_high - eta_low ) , 2 ) ;

	return estimate;
}