double CoalescenceDistribution(vector<double> par_array)
{
  double rT1   = par_array.at(0);
  double rT2   = par_array.at(1);
  double phi1  = par_array.at(2);
  double phi2  = par_array.at(3);
  double eta1  = par_array.at(4);
  double eta2  = par_array.at(5);
  double Phi1  = par_array.at(6);
  double Phi2  = par_array.at(7);
  double y1    = par_array.at(8);
  double y2    = par_array.at(9);
  double pT1   = par_array.at(10);
  double pT2   = par_array.at(11);
  double m1    = par_array.at(12);
  double m2    = par_array.at(13);
  double tau   = par_array.at(14);
  double T     = par_array.at(15);
  double sigma = par_array.at(16);

  double mT1 = sqrt(m1*m1 + pT1*pT1);
  double mT2 = sqrt(m2*m2 + pT2*pT2);
  double x2_exponent = rT1*rT1 + rT2*rT2 - 2.*(rT1*rT2*cos(phi1-phi2) + tau*tau*(cosh(eta1-eta2)-1.) );  // (x1-x2)^2
  double p2_exponent = mT1*mT1 + mT2*mT2 + pT1*pT1 + pT2*pT2 + (m1-m2)*(m1-m2) - 2.*(mT1*mT2*cosh(y1-y2) + pT1*pT2*cos(Phi1-Phi2));  // (p1-p2)^2

  return 9 * M_PI / pow( sigma , 6. ) * Theta( x2_exponent , sigma*sigma ) * Theta( p2_exponent + (m1-m2)*(m1-m2) , sigma * sigma ) ;
  // return 1. / pow( sigma , 6. ) * exp( - x2_exponent / sigma / sigma - p2_exponent * sigma * sigma ) ;
}

double dNi(vector<double> par_array)
{
  double rTi  = par_array.at(0);
  double phii = par_array.at(1);
  double etai = par_array.at(2);
  double Phii = par_array.at(3);
  double yi   = par_array.at(4);
  double pTi  = par_array.at(5);
  double mi   = par_array.at(6);
  double tau  = par_array.at(7);
  double T    = par_array.at(8);
  double dyi  = par_array.at(9);

  double mTi = sqrt(mi*mi + pTi*pTi);
  double Cyi = exp( - yi * yi / dyi / dyi );
  double pui = mTi * cosh( yi - etai ) ;
  double gammaTi = 1. / sqrt( 1 - rTi * rTi / tau / tau ) ;
  
  return Cyi * pui * exp( - gammaTi / T * ( pui - pTi * rTi / tau * cos( Phii - phii ) ) ) ;
}

double integrand( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_vars_array )
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
  double y1    = integ_vars_array.at(4);
  double y2    = integ_vars_array.at(5);
  double Phi1  = integ_vars_array.at(6);
  double Phi2  = integ_vars_array.at(7);
  double eta1  = integ_vars_array.at(8);
  double eta2  = integ_vars_array.at(9);
  double phi1  = integ_vars_array.at(10);
  double phi2  = integ_vars_array.at(11);
  double rT1   = integ_vars_array.at(12);
  double rT2   = integ_vars_array.at(13);

  double mT1 = sqrt(m1*m1 + pT1*pT1);
  double mT2 = sqrt(m2*m2 + pT2*pT2);
  double mT = sqrt(m*m + pT*pT);

  vector<double> coal_pars = {rT1,rT2,phi1,phi2,eta1,eta2,Phi1,Phi2,y1,y2,pT1,pT2,m1,m2,tau,T,sigma};
  vector<double> momdist1_pars = {rT1,phi1,eta1,Phi1,y1,pT1,m1,tau,T,dy1} ;
  vector<double> momdist2_pars = {rT2,phi2,eta2,Phi2,y2,pT2,m2,tau,T,dy2} ;

  double result = 0.0; 
  if ( DiracDelta( pT * sin(Phi) - ( pT1 * sin(Phi1) + pT2 * sin(Phi2) ) , 0.001) > 0. )
    if ( DiracDelta( mT * sinh(Y) - ( mT1 * sinh(y1) + mT2 * sinh(y2) ) , 0.001) > 0. ) 
      result = dNi( momdist1_pars ) * dNi( momdist2_pars ) * CoalescenceDistribution( coal_pars ) ;
      // if ( (y1 > 0 && y2 > 0) || (y1 < 0 && y2 < 0) )
        // if ( (eta1 > 0 && eta2 > 0) || (eta1 < 0 && eta2 < 0) )
      
  return result ;
}


double rapidity_momentum_distribution( vector<double> func_vars_array , vector<double> par_array , vector<double> integ_lims_array , long int N_iter )
{
	double totalSum, functionVal = 0;
	int iter = 0;

  double low_deg = integ_lims_array.at(0) ;
  double up_deg = integ_lims_array.at(1) ;

  double low_cor = integ_lims_array.at(2) ;
  double up_cor = integ_lims_array.at(3) ;

  double low_mom = integ_lims_array.at(4) ;
  double up_mom = integ_lims_array.at(5) ;

  double low_rap = integ_lims_array.at(6) ;
  double up_rap = integ_lims_array.at(7) ;

  double low_eta = integ_lims_array.at(8) ;
  double up_eta = integ_lims_array.at(9) ;

  vector<double> integ_vars_array = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. } ;

	while( iter < N_iter - 1 )
	{
    integ_vars_array.at(0)  = low_deg + ( double( rand() ) / RAND_MAX ) * ( up_deg - low_deg );  // Phi
    integ_vars_array.at(1)  = low_cor + ( double( rand() ) / RAND_MAX ) * ( up_cor - low_cor );  // rT
    integ_vars_array.at(2)  = low_mom + ( double( rand() ) / RAND_MAX ) * ( up_mom - low_mom );  // pT1
    integ_vars_array.at(3)  = low_mom + ( double( rand() ) / RAND_MAX ) * ( up_mom - low_mom );  // pT2
    integ_vars_array.at(4)  = low_rap + ( double( rand() ) / RAND_MAX ) * ( up_rap - low_rap );  // y1
    integ_vars_array.at(5)  = low_rap + ( double( rand() ) / RAND_MAX ) * ( up_rap - low_rap );  // y2
    integ_vars_array.at(6)  = low_deg + ( double( rand() ) / RAND_MAX ) * ( up_deg - low_deg );  // Phi1
    integ_vars_array.at(7)  = low_deg + ( double( rand() ) / RAND_MAX ) * ( up_deg - low_deg );  // Phi2
    integ_vars_array.at(8)  = low_eta + ( double( rand() ) / RAND_MAX ) * ( up_eta - low_eta );  // eta1
    integ_vars_array.at(9)  = low_eta + ( double( rand() ) / RAND_MAX ) * ( up_eta - low_eta );  // eta2
    integ_vars_array.at(10) = low_deg + ( double( rand() ) / RAND_MAX ) * ( up_deg - low_deg );  // phi1
    integ_vars_array.at(11) = low_deg + ( double( rand() ) / RAND_MAX ) * ( up_deg - low_deg );  // phi2
    integ_vars_array.at(12) = low_cor + ( double( rand() ) / RAND_MAX ) * ( up_cor - low_cor );  // rT1
    integ_vars_array.at(13) = low_cor + ( double( rand() ) / RAND_MAX ) * ( up_cor - low_cor );  // rT2

		functionVal = integrand( func_vars_array , par_array , integ_vars_array ) ;

    if(iter%10000000==0) cerr << "\t --> Y = " << func_vars_array.at(0) << " , pT = " << func_vars_array.at(1) << " iter = " << iter/1e6 << "M" << endl;

		totalSum += functionVal;

		iter++;
	}

	double estimate = totalSum / N_iter * pow((up_deg - low_deg),5) * pow((up_cor - low_cor),3) * pow((up_mom - low_mom),2) * pow((up_rap - low_rap),2) * pow((up_eta - low_eta),2) ;

	return estimate;

}