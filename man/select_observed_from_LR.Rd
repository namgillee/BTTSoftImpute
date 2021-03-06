% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_observed_from_LR.R
\name{select_observed_from_LR}
\alias{select_observed_from_LR}
\title{select_observed_from_LR}
\usage{
select_observed_from_LR(x, numobs = NULL, rate = 1)
}
\arguments{
\item{x}{A data matrix in LR (Low-rank (SVD)) format with attributes $u, $d, $v.
$u is a matrix, x$v is a blockTT format.}

\item{numobs}{Number of observed values to select.}

\item{rate}{If numobs==NULL, then numobs = rate * numels(full(x)).}
}
\value{
The returned value is a SparseMatrix format, with
        attributes $i, $j, $v, $nrow, $ncol, $dimnames.
}
\description{
Select observed values randomly from data x.
Here, 'randomly' means uniform random sampling.
}
