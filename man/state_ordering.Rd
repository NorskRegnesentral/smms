% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm_msm_preprocessing_func.R
\name{state_ordering}
\alias{state_ordering}
\title{Produce an ordering of the states from a multi-state graph}
\usage{
state_ordering(graph)
}
\arguments{
\item{graph}{A directed, acyclic graph in the igraph format (igraph package).}
}
\value{
A data frame with 3 columns giving a mapping from the state names (first column) to a partial ordering of the states
(second column, from 0 to k, and where a higher number indicates that a state is further removed from the initial state).
The last column indicates whether each state is initial, absorbing or transient.
}
\description{
Produce an ordering of the states from a multi-state graph
}
