#include <R.h>
#include <Rmath.h>
/* Fortran function for standard normal PDF, CDF */
double F77_SUB(fpnorm)(double *x)
{
  return pnorm(*x, 0.0, 1.0, 1, 0); 
}

