
<!-- README.md is generated from README.Rmd. Please edit that file -->

# prolong

`prolong` is an implementation of PROLONG, a method designed for
regressing a longitudinal phenotype on longitudinal omics data. PROLONG
leverages the first differences of the data to address the piece-wise
linear structure and the observed time dependence, using a Laplacian
network constraint to incorporate the correlation structure of the
predictors, and a group lasso constraint induces sparsity while grouping
metabolites across their first differenced observations.

`prolong` currently supports continuous response variables and
continuous omics data. Next updates will expand to continuous multiomic
and demographic data, later updates will expand to discrete multiomics
and eventually discrete response variables.

## Installation

You can install the development version of prolong from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("stevebroll/prolong")
```
