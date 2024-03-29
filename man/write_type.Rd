% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sm_msm_latex_func.R
\name{write_type}
\alias{write_type}
\title{Write likelihood contribution for one observation type to latex type format}
\usage{
write_type(obs_type, graph, abs_exact = TRUE)
}
\arguments{
\item{obs_type}{The observation type (as a string).}

\item{graph}{A directed, acyclic graph in the igraph format (igraph package).}

\item{abs_exact}{A boolean indicating whether the time of entrance into absorbing states is observed
exactly (TRUE) or not (FALSE). Default value is TRUE.}
}
\value{
A text string with the likelihood contribution formula written in a latex type format.
}
\description{
Write likelihood contribution for one observation type to latex type format
}
