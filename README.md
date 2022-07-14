<!-- README.md is generated from README.Rmd. Please edit that file -->

The `smms` package allows you to fit Semi-Markovian multi-state models
to panel datasets. The package constructs and optimises the likelihood
of arbitrary multi-state models where the possible state transitions can
be described by an acyclic graph with one or more initial states and one
or more absorbing states.

Designed for data where the exact state transitions are not necessarily
observed, so that we have interval censoring of the transition times.
Data like these are sometimes referred to as panel data.

The methodology is explained in the paper [A new framework for
semi-Markovian parametric multi-state models with interval
censoring](https://www.mn.uio.no/math/english/research/projects/focustat/publications_2/multistate_final_july2022.pdf)
by Aastveit, Cunen and Hjort (2022). The code used for the application
and simulations in the paper is available in the `scripts` folder.

# Installation

In order to directly install the package from github, you need the
package `devtools`.

``` r
install.packages("devtools")
devtools::install_github("NorskRegnesentral/smms")
```

<!--
# Whats new

### Version xx


See NEWS.md for changes for complete version history.
-->

# Overview

The main function in this package is the `smms` function, which is used
to fit a multi-state model. The user needs to provide a datasets of
states and time-points for multiple units of observation, typically
referred to as *patients*, a graph describing the states and the
possible transitions between them, and a set of parametric densities,
and survival functions.

In the next section, we will illustrate the use of the `smms` function
in a simple example.

# A simple example

We will use the CAV dataset from the `msm` package (Jackson 2011) as an
illustration. The dataset monitors a number of patients for a number of
years after heart transplantation. Coronary allograft vasculopathy (CAV)
is a condition potentially occurring after hear transplantation. At each
time-point the patients are assigned to one of four states: well, mild
CAV, severe CAV and death. The time of death is recorded precisely, but
the times of entrance into the CAV-states are interval censored

# References

-   Jackson CH (2011). [Multi-State Models for Panel Data: The msm
    Package for R.](https://www.jstatsoft.org/v38/i08/) Journal of
    Statistical Software, 38, 1–29.
