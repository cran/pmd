% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmdda.R
\name{gettarget}
\alias{gettarget}
\title{Get multiple injections index for selected retention time}
\usage{
gettarget(rt, drt = 10, n = 6)
}
\arguments{
\item{rt}{retention time vector for peaks in seconds}

\item{drt}{retention time drift for targeted analysis in seconds, default 10.}

\item{n}{max ions numbers within retention time drift windows}
}
\value{
index for each injection
}
\description{
Get multiple injections index for selected retention time
}
\examples{
data(spmeinvivo)
pmd <- getpaired(spmeinvivo)
std <- getstd(pmd)
index <- gettarget(std$rt[std$stdmassindex])
table(index)
}
