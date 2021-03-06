% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmdvis.R
\name{plotstd}
\alias{plotstd}
\title{Plot the std mass from GlobalStd algorithm}
\usage{
plotstd(list)
}
\arguments{
\item{list}{a list from getstd function}
}
\description{
Plot the std mass from GlobalStd algorithm
}
\examples{
data(spmeinvivo)
pmd <- getpaired(spmeinvivo)
std <- getstd(pmd)
plotstd(std)
}
\seealso{
\code{\link{getstd}}, \code{\link{globalstd}}
}
