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
```

They respectively deal with:

  - the computations of the coefficients of $K_{4,n}$
  and $\lambda_{3,n}$
  
  - the production of Figures 1 and 2
  
  - the computation of Table 2
  
  - the optimization of the integrals for Appendices B.1 and B.2.


# Requirements

This study was done using R version 4.2.1 with the following packages:

| Package          | Version |
| :--------------- | :------ |
| `BoundEdgeworth` | 0.1.0   |
| `tidyverse`      | 1.3.2   |
| `ggplot2`        | 3.3.6   |
| `tibble`         | 3.1.8   |
| `tidyr`          | 1.2.0   |
| `readr`          | 2.1.2   |
| `purrr`          | 0.3.4   |
| `dplyr`          | 1.0.9   |
| `stringr`        | 1.4.1   |
| `forcats`        | 0.5.2   |
| `cubature`       | 0.4.5   |
| `expint`         | 0.1-7   |
| `rmarkdown`      | 2.16.1  |
| `gt`             | 0.7.0   |


They can be all installed using the following command.

``` r
install.packages(c("BoundEdgeworth", "tidyverse", "ggplot2",
  "tibble", "tidyr", "readr", "purrr" , "dplyr",
  "stringr", "forcats", "cubature", "expint", "rmarkdown", "gt"))
```
