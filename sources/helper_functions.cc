#include "../includes/model_funcs.h"

bool model_funcs::Theta(double x, double y)
{
  return ( x > y );
}

double model_funcs::DiracDelta(double x, double eps)
{
  if( abs(x) > eps)
    return x;
  else
    return 0;
}

// double DiracDelta(double x, double eps)
// {
  // return exp( -x * x / eps / eps );
// }