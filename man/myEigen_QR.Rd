% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myRFunctions.R
\name{myEigen_QR}
\alias{myEigen_QR}
\title{Perform PCA on A}
\usage{
myEigen_QR(A, numIter = 1000)
}
\arguments{
\item{A}{is a square matrix, numIter (default 1000)}
}
\value{
a list with D and V
}
\description{
Perform PCA on A
}
\examples{
myEigen_QR(A,1000)
}
