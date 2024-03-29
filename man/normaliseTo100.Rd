% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{normaliseTo100}
\alias{normaliseTo100}
\title{Normalise read counts to 100, i.e. in percent relative abundance per sample}
\usage{
normaliseTo100(data)
}
\arguments{
\item{data}{(\emph{required}) Data list as loaded with \code{\link{amp_load}}.}
}
\value{
A modified ampvis2 object
}
\description{
Normalise read counts to 100, i.e. in percent relative abundance per sample
}
\examples{
data("AalborgWWTPs")
AalborgWWTPs
normalised <- normaliseTo100(AalborgWWTPs)
normalised
}
