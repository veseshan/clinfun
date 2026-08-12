# Stopping rule for toxicity/futility monitoring

Computes a stopping rule and its operating characteristics for toxicity
monitoring based repeated significance testing.

## Usage

``` r
toxbdry(pLo, pHi, n, cP0=0.1, cP1=0.9, ngrid=6, niter=10, delta=0,
          priority=c("null","alt"))
futilbdry(rLo, rHi, n, size=0.1, power=0.9, ngrid=3, niter=10, delta=0.5)
bdrycross.prob(n, r, ptox)
# S3 method for class 'toxbdry'
print(x, ...)
# S3 method for class 'futilbdry'
print(x, ...)
```

## Arguments

- pLo:

  the toxicity rate that is acceptable.

- rLo:

  baseline (null) response rate.

- pHi:

  the toxicity rate that is too high and hence unacceptable.

- rHi:

  desirable response rate. Stop when it is too unlikely.

- n:

  vector of times (sample size) when toxicty/response is monitored.

- r:

  vector of maximum acceptable toxicities (non-responders for futility)
  corresponding to n.

- ptox:

  the toxicity rates for which the operating characteristics are
  calculated. For futility this is the non-response rate.

- cP0:

  boundary crossing probability under pLo i.e. type I error or the
  probability of declaring a treatment with toxicity rate pLo
  unacceptable.

- cP1:

  boundary crossing probability under pHi i.e. power or the probability
  of declaring a treatment with toxicity rate pHi unacceptable.

- size:

  probability of calling drug effective if response rate is rLo.

- power:

  probability of calling drug effective if response rate is rHi.

- ngrid:

  the number of toxicity rates from pLo to pHi for which the operating
  characteristics are computed.

- niter:

  the number of iterations run to obtain the boundary.

- delta:

  power determining the shape of the boundary. Should be between 0
  (default) and 0.5.

- priority:

  the error threshold to prioritize when the max sample size is too
  small to have both error thresholds satisfied. Default is the null
  i.e. error under pLo.

- x:

  object returned by the function toxbdry.

- ...:

  additional arguments to print.

## Value

the function returns a list with:

- looks:

  when toxicty is monitored - same as input n.

- lo.bdry:

  lower boundary is a vector of maximum acceptable number of toxicities
  corresponding the number of subjects in n. The boundary crossing
  probability for this is slightly above cP0.

- hi.bdry:

  upper boundary is a vector of maximum acceptable number of toxicities
  corresponding the number of subjects in n. The boundary crossing
  probability for this is slightly below cP0.

- bdry.oc:

  the operating characteristics i.e the toxicity rate, the probability
  of crossing, stopping (i.e. cross before the last observation) and the
  expected sample size for both the low (lo) and high (hi) boundaries.

- bdry.alpha:

  the alpha levels for testing at each look for the two boundaries.

stopping for toxicity is done when the number of toxicities exceeds the
boundary i.e. the boundary gives the maximum acceptable number.

## Details

The shape parameter delta is used to determine the stopping boundary
with 0 corresponding to the Pocock boundary where the same significance
level is used for all looks and 0.5 corresponding to the O'Brien-Fleming
boundary which has smaller probability of stopping at early looks.

Default value of delta for toxicity monitoring is 0; value between 0.1
and 0.2 is a reasonable choice to make it less likely to stop early.
Default values of delta for futility stopping is 0.5.

For toxicity monitoring two sets of probabilities - pstop and pcross -
are given which correspond to probability of stopping early and
probability of declaring the treatment too toxic with the full
complement of study subjects accrued and treated.

For futility monitoring instead two sets of probabilities - pstop and
peffective - are given corresponding to the probability of stopping
early for futility and probability of finishing the trial and declaring
it a success. Note that peffective is the complement of pcross in
toxicity monitoring.

The futility boundary can have a -1 in earlier looks which means that
even zero responses is not sufficient for stopping at that look.

The exact calculations in this function are done along the lines of the
method in Chapter 12 of Jennison and Turnbull (2000). Ivanova, Qaqish
and Schell (2005) have an illustrative paper.

## References

Jennison C and Turnbull BW. (2000). Group Sequential Methods with
Applications to Clinical Trials. *Chapman and Hall/CRC*

Ivanova A, Qaqish BF and Schell MJ. (2005). Continuous Toxicity
Monitoring in Phase II Trials in Oncology. *Biometrics* 61, 540-545.

## Examples

``` r
  toxbdry(0.2, 0.35, c(20,40,60,75))
#> 
#>  Toxicity boundary based on repeated significance testing 
#>     Boundary shape parameter delta = 0 
#> 
#>  ******************************************************************
#>  * Stop if the number of toxicities exceeds (i.e. >) the boundary *
#>  ******************************************************************
#> 
#>   Low boundary:   7/20   12/40   17/60   20/75   
#>  High boundary:   7/20   12/40   17/60   21/75   
#> 
#>  Operating Characteristics: 
#> 
#>      ptox pcross.lo pstop.lo ess.lo pcross.hi pstop.hi ess.hi
#> [1,] 0.20     0.101    0.079   72.0     0.088    0.079   72.0
#> [2,] 0.23     0.247    0.188   68.1     0.214    0.188   68.1
#> [3,] 0.26     0.456    0.354   62.1     0.406    0.354   62.1
#> [4,] 0.29     0.672    0.548   54.7     0.619    0.548   54.7
#> [5,] 0.32     0.838    0.727   47.0     0.797    0.727   47.0
#> [6,] 0.35     0.935    0.859   39.9     0.911    0.859   39.9
  toxbdry(0.2, 0.3, c(20,40,60,75), cP0=0.15, cP1=0.8)
#> Warning: Max sample size 75 may be small for the specified stopping probabilities
#> 
#>  Toxicity boundary based on repeated significance testing 
#>     Boundary shape parameter delta = 0 
#> 
#>  ******************************************************************
#>  * Stop if the number of toxicities exceeds (i.e. >) the boundary *
#>  ******************************************************************
#> 
#>   Low boundary:   6/20   12/40   16/60   20/75   
#>  High boundary:   7/20   12/40   16/60   20/75   
#> 
#>  Operating Characteristics: 
#> 
#>      ptox pcross.lo pstop.lo ess.lo pcross.hi pstop.hi ess.hi
#> [1,] 0.20     0.152    0.140   69.1     0.117    0.103   71.6
#> [2,] 0.22     0.253    0.229   65.7     0.213    0.186   69.0
#> [3,] 0.24     0.381    0.342   61.4     0.341    0.297   65.4
#> [4,] 0.26     0.522    0.469   56.5     0.487    0.427   61.0
#> [5,] 0.28     0.659    0.598   51.2     0.632    0.562   56.1
#> [6,] 0.30     0.777    0.715   45.9     0.758    0.688   51.0
  # continuous monitoring
  toxbdry(0.1, 0.3, 2:30)
#> 
#>  Toxicity boundary based on repeated significance testing 
#>     Boundary shape parameter delta = 0 
#> 
#>  ******************************************************************
#>  * Stop if the number of toxicities exceeds (i.e. >) the boundary *
#>  ******************************************************************
#> 
#>   Low boundary:   1/3   2/8   3/13   4/20   5/26   6/30   
#>  High boundary:   1/3   2/8   3/13   4/19   5/26   6/30   
#> 
#>  Operating Characteristics: 
#> 
#>      ptox pcross.lo pstop.lo ess.lo pcross.hi pstop.hi ess.hi
#> [1,] 0.10     0.100    0.099   28.1     0.097    0.096   28.1
#> [2,] 0.14     0.250    0.245   25.7     0.245    0.239   25.8
#> [3,] 0.18     0.448    0.437   22.4     0.441    0.430   22.5
#> [4,] 0.22     0.644    0.631   18.8     0.638    0.624   19.0
#> [5,] 0.26     0.799    0.786   15.5     0.795    0.781   15.6
#> [6,] 0.30     0.901    0.891   12.6     0.898    0.888   12.7
  # prioritize cP1 error threshold
  toxbdry(0.1, 0.3, 2:25, priority="alt")
#> Warning: Max sample size 25 may be small for the specified stopping probabilities
#> 
#>  Toxicity boundary based on repeated significance testing 
#>     Boundary shape parameter delta = 0 
#> 
#>  ******************************************************************
#>  * Stop if the number of toxicities exceeds (i.e. >) the boundary *
#>  ******************************************************************
#> 
#>   Low boundary:   1/4   2/10   3/16   4/23   5/25   
#>  High boundary:   1/4   2/10   3/16   4/22   5/25   
#> 
#>  Operating Characteristics: 
#> 
#>      ptox pcross.lo pstop.lo ess.lo pcross.hi pstop.hi ess.hi
#> [1,] 0.10     0.146    0.145   22.8     0.140    0.139   22.8
#> [2,] 0.14     0.312    0.310   20.5     0.300    0.296   20.5
#> [3,] 0.18     0.504    0.499   17.8     0.488    0.480   17.8
#> [4,] 0.22     0.679    0.673   15.0     0.663    0.652   15.1
#> [5,] 0.26     0.813    0.807   12.5     0.800    0.790   12.5
#> [6,] 0.30     0.902    0.897   10.3     0.893    0.885   10.3
```
