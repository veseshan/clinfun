# Simon's 2-stage Phase II design

Calculates Optimal and Minimax 2-stage Phase II designs given by Richard
Simon

## Usage

``` r
ph2simon(pu, pa, ep1, ep2, nmax=100)
# S3 method for class 'ph2simon'
print(x, ...)
# S3 method for class 'ph2simon'
plot(x, ...)
```

## Arguments

- pu:

  unacceptable response rate; baseline response rate that needs to be
  exceeded for treatment to be deemed promising

- pa:

  response rate that is desirable; should be larger than pu

- ep1:

  threshold for the probability of declaring drug desirable under pu
  (target type 1 error rate); between 0 and 1

- ep2:

  threshold for the probability of rejecting the drug under pa (target
  type 2 error rate); between 0 and 1

- nmax:

  maximum total sample size (default 100; can be at most 1000)

- x:

  object returned by ph2simon

- ...:

  arguments to be passed onto plot and print commands called within

## Value

ph2simon returns a list with pu, pa, alpha, beta and nmax as above and:

- out:

  matrix of best 2 stage designs for each value of total sample size n.
  The 6 columns in the matrix are:

  |         |                                                            |
  |---------|------------------------------------------------------------|
  | r1      | number of responses needed to exceeded in first stage      |
  | n1      | number of subjects treated in first stage                  |
  | r       | number of responses needed to exceeded at the end of trial |
  | n       | total number of subjects to be treated in the trial        |
  | EN(pu)  | expected number pf patients in the trial under pu          |
  | PET(pu) | probability of stopping after the first stage under pu     |

Trial is stopped early if \<= r1 responses are seen in the first stage
and treatment is considered desirable only when \>r responses seen.

## Methods (by generic)

- `print(ph2simon)`: formats and returns the minimax, optimal and any
  admissible designs.

- `plot(ph2simon)`: plots the expected sample size against the maximum
  sample size as in Jung et al., 2001

## References

Simon R. (1989). Optimal Two-Stage Designs for Phase II Clinical Trials.
*Controlled Clinical Trials* 10, 1-10.

Jung SH, Carey M and Kim KM. (2001). Graphical Search for Two-Stage
Designs for Phase II Clinical Trials. *Controlled Clinical Trials* 22,
367-372.

Jung SH, Lee T, Kim K, and George, SL. (2004). Admissible two-stage
designs for phase II cancer clinical trials. *Statistics in medicine*
23(4), 561-569.

## See also

[`twostage.inference`](twostage.inference.md),
[`oc.twostage.bdry`](oc.twostage.bdry.md)

## Examples

``` r
  ph2simon(0.2, 0.4, 0.1, 0.1)
#> 
#>  Simon 2-stage Phase II design 
#> 
#> Unacceptable response rate:  0.2 
#> Desirable response rate:  0.4 
#> Error rates: alpha =  0.1 ; beta =  0.1 
#> 
#>         r1 n1  r  n EN(p0) PET(p0)   qLo   qHi
#> Minimax  3 19 10 36  28.26  0.4551 0.691 1.000
#> Optimal  3 17 10 37  26.02  0.5489 0.000 0.691
#> 
  ph2simon(0.2, 0.35, 0.05, 0.05)
#> 
#>  Simon 2-stage Phase II design 
#> 
#> Unacceptable response rate:  0.2 
#> Desirable response rate:  0.35 
#> Error rates: alpha =  0.05 ; beta =  0.05 
#> 
#>         r1 n1  r  n EN(p0) PET(p0)   qLo   qHi
#> Minimax 15 68 25 95  75.44  0.7244 0.582 1.000
#> Optimal 11 53 26 99  69.88  0.6330 0.000 0.582
#> 
#> Warning:   Optimal sample size too close to nmax. 
#>   Try increasing nmax (current value = 100)
  ph2simon(0.2, 0.35, 0.05, 0.05, nmax=150)
#> 
#>  Simon 2-stage Phase II design 
#> 
#> Unacceptable response rate:  0.2 
#> Desirable response rate:  0.35 
#> Error rates: alpha =  0.05 ; beta =  0.05 
#> 
#>            r1 n1  r   n EN(p0) PET(p0)   qLo   qHi
#> Minimax    15 68 25  95  75.44  0.7244 0.582 1.000
#> Admissible 11 53 26  99  69.88  0.6330 0.306 0.582
#> Admissible 11 51 27 104  67.68  0.6852 0.109 0.306
#> Optimal    11 50 28 109  67.07  0.7107 0.000 0.109
#> 
```
