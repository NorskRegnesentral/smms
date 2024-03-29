% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm_msm_likelihood_func.R
\name{names_of_survival_density}
\alias{names_of_survival_density}
\title{Produce a table with each edge of the graph and the corresponding names of survival functions
and densities}
\usage{
names_of_survival_density(graph)
}
\arguments{
\item{graph}{A directed, acyclic graph in the igraph format (igraph package).}
}
\value{
A data frame with one row for each edge in the graph. First column gives the standard edge
name, second column indicates the names of the survival functions, third column gives the names of the
density functions. The two next columns give the names of the source node and the target node \emph{in the
original user-defined names}. The last column gives the type of edge: either "abs" for absorbing if the target node
is an absorbing state, or "trans" for transient if not.
}
\description{
Produce a table with each edge of the graph and the corresponding names of survival functions
and densities
}
