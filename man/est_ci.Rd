% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm_msm_optim_func.R
\name{est_ci}
\alias{est_ci}
\title{Estimates and confidence intervals}
\usage{
est_ci(param, hessian, level = 0.95, pos = FALSE)
}
\arguments{
\item{param}{The parameter estimates.}

\item{hessian}{The hessian matrix.}

\item{level}{The level of confidence. Default value is 0.95.}

\item{pos}{Boolean denoting whether the estimates and intervals should live on the real line (default), or on the postive half-line (pos=TRUE).}
}
\value{
A table with one row per estimate, and 3 columns (estimate, lower.ci, upper.ci).
}
\description{
Provide estimates and approximate confidence intervals for all parameters
(on the "original" scale - meaning that all parameters live on -Inf to +Inf).
}
