% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm_msm_optim_func.R
\name{occupancy_prob}
\alias{occupancy_prob}
\title{Compute the occupancy probability over time}
\usage{
occupancy_prob(state, time, param, graph, xval = NULL)
}
\arguments{
\item{state}{A string with the name of a state for which to compute the probability (user-defined names).}

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
The last state is often heavy to compute, but recall that is has to be equal to 1 minus the others
}
