#include <R.h>
#include <Rmath.h>
/* Fortran function for standard normal PDF, CDF */
double F77_SUB(fdnorm)(double *x)
{
  return dnorm(*x, 0.0, 1.0, 0); 
}
