% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_minus_lr_btt.R
\name{sp_minus_lr_btt}
\alias{sp_minus_lr_btt}
\title{sp_minus_lr_btt}
\usage{
sp_minus_lr_btt(A, Z)
}
\arguments{
\item{A}{SParse (SP) format with elements i, j, v, nrow, ncol}

\item{Z}{LR format with elements u, d, and v.
Z$v is in BlockTT format.}
}
\description{
This computes a sparse (SP) matrix minus
low-rank (LR) matrix operation.
}
\details{
Last modified: 2017.08.25. by Namgil Lee
}
