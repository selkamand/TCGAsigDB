
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TCGAsigDB

**This project is a work in progress and not yet ready for use!!**

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/TCGAsigDB)](https://CRAN.R-project.org/package=TCGAsigDB)
<!-- badges: end -->

The goal of TCGAsigDB is to provide access to the results of TCGA
signature analysis.

## Installation

``` r
# install.packages('remotes')
remotes::install_github('selkamand/TCGAsigDB')
```

## Quick Start

``` r
library(TCGAsigDB)

create_tcga_pancan_database(ids = c('TCGA-CA-6717-01', 'TCGA-A2-A0T5-01', 'TCGA-CF-A9FF-01'))
```
