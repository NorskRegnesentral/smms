% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm_msm_optim_func.R
\name{overall_survival}
\alias{overall_survival}
\title{Compute the overall survival}
\usage{
overall_survival(time, param, graph, xval = NULL)
}
\arguments{
\item{time}{A vector of timepoints for which to compute the occupancy probability.}

\item{param}{The parameter values in which the probabilities should be evaluated. The dimension
and ordering is given by the user-specified densities.}

\item{graph}{A directed, acyclic graph giving the multistate structure, in the igraph format (igraph package).}

\item{xval}{A vector of covariate values.}
}
\value{
A vector of the same length as time.
}
\description{
...
}
