% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pos.R
\name{pos}
\alias{pos}
\title{Pseudotime-temporal ordering score}
\usage{
pos(pseudotime, truetime)
}
\arguments{
\item{pseudotime}{An n-vector of pseudotime estimates}

\item{truetime}{An n-vector of true timepoint values}
}
\value{
The positional ordering score
}
\description{
Pseudotime-temporal ordering (POS) is a metric introduced in the TSCAN paper, used to compare the estimated ordering of cells to the true ordering.
}
