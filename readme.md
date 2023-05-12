# Replication for the article 'Explicit non-asymptotic bounds for the distance to the first-order Edgeworth expansion'

Alexis Derumigny, Lucas Girard and Yannick Guyonvarch

---

This file describes the procedure to use in order to replicate the numerical results
of the article 'Explicit non-asymptotic bounds for the distance to the first-order Edgeworth expansion'.

# Summary of the files

These computations have been separated into different `.Rmd`
files. They can be run using the following commands:

``` r
rmarkdown::render("CoefficientsMainTheorems.Rmd")
rmarkdown::render("Comparison_Bounds_EE.Rmd")
rmarkdown::render("MinimumInformativeSize.Rmd")
rmarkdown::render("OptimizationIntegrals.Rmd")
rmarkdown::render("BoundingCDFs.Rmd")
rmarkdown::render("SufficientSampleSizes.Rmd")
rmarkdown::render("Computation_kappa.Rmd")
```

They respectively deal with:

  - the computations of the coefficients of $K_{4,n}$
  and $\lambda_{3,n}$
  
  - the production of Figures 1 and 2
  
  - the computation of Table 2
  
  - the optimization of the integrals for Appendices B.1 and B.2.
  
  - the production of Figures 3, 4 and 5
  
  - the computation of Table 3 and 4
  
  - the optimization of the values of $\kappa$ for several distributions


# Requirements

This study was done using R version 4.3.0 with the following packages:

| Package          | Version |
| :--------------- | :------ |
| `BoundEdgeworth` | 0.1.2   |
| `tidyverse`      | 2.0.0   |
| `ggplot2`        | 3.4.2   |
| `tibble`         | 3.2.1   |
| `tidyr`          | 1.3.0   |
| `readr`          | 2.1.4   |
| `purrr`          | 1.0.1   |
| `dplyr`          | 1.1.2   |
| `stringr`        | 1.5.0   |
| `forcats`        | 1.0.0   |
| `cubature`       | 2.0.4.6 |
| `expint`         | 0.1.8   |
| `rmarkdown`      | 2.21    |
| `gt`             | 0.9.0   |


They can be all installed using the following command.

``` r
install.packages(c("BoundEdgeworth", "tidyverse", "ggplot2",
  "tibble", "tidyr", "readr", "purrr" , "dplyr",
  "stringr", "forcats", "cubature", "expint", "rmarkdown", "gt"))
```
