#include "../includes/model_functions.h"

bool Theta(double x, double limit)
{
  return ( x < limit );
}

double DiracDelta(double x, double eps)
{
  // return exp( -x * x / eps / eps );
  if ( abs(x) > eps )
    return x;
  else
    return 0;
}