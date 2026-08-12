# Comparing the AUC from ROC curves from nested binary regression

Testing the differences in the area under the ROC curves (AUC) from
nested maximum rank correlation models.

## Usage

``` r
deltaAUC(y, x, z)
```

## Arguments

- y:

  binary response variable

- x:

  matrix of set of covariates that is the basis of the existing
  (reduced) model

- z:

  matrix of set of covariates that are added to to get the new (full)
  model

## Value

It returns a list with the following elements

- par.full:

  the MRC estimate of parameters for the full model

- par.red:

  the MRC estimate of parameters for the reduced model

- results:

  matrix od results which gives the full reduced model AUCs along with
  the test statistic and p-value

## Details

The models are fit using maximum rank correlation (MRC) method which is
an alternate approach to logistic regression. In MRC the area under the
ROC curve (AUC) is maximized as opposed to the likelihood in logistic
regression. Due to invariance of AUC to location and scale shifts one of
the parameters (anchor variable) is set to 1.

The first variable (column) in x is used as the anchor variable.

The IPMN data set used as an example in the paper below is included. The
columns are high risk lesion (V1), recent weight loss (V2), main duct
involvement (V4), presence of a solid component in imaging (V3), and
lesion size (V5).

## References

Heller G., Seshan V.E., Moskowitz C.S. and Gonen M. (2016) Inference for
the difference in the area under the ROC curve derived from nested
binary regression models. *Biostatistics* 18, 260-274.

## Examples

``` r
  data(ipmn)
  deltaAUC(ipmn$V1, cbind(ipmn$V4, ipmn$V3, ipmn$V5), ipmn$V2)
#> $par.full
#> [1] 1.0000000 0.1685387 0.7030273 0.4309026
#> 
#> $par.red
#> [1] 1.0000000 0.3973677 0.6374656
#> 
#> $results
#>             AUCfull AUCreduced statistic p-value
#> Empirical 0.8045688  0.7806801  9.842135   5e-05
#> Smooth    0.8007295  0.7777246  9.478028   8e-05
#> 
```
