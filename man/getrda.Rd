% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sda.R
\name{getrda}
\alias{getrda}
\title{Perform structure/reaction directed analysis for mass only.}
\usage{
getrda(mz, freqcutoff = 10, digits = 3, top = 20, formula = NULL)
}
\arguments{
\item{mz}{numeric vector for independant mass or mass to charge ratio. Mass to charge ratio from GlobalStd algorithm is suggested. Isomers would be excluded automately}

\item{freqcutoff}{pmd freqency cutoff for structures or reactions, default 10}

\item{digits}{mass or mass to charge ratio accuracy for pmd, default 3}

\item{top}{top n pmd freqency cutoff when the freqcutoff is too small for large data set}

\item{formula}{vector for formula when you don't have mass or mass to charge ratio data}
}
\value{
logical matrix with row as the same order of mz or formula and column as high freqency pmd group
}
\description{
Perform structure/reaction directed analysis for mass only.
}
\examples{
data(spmeinvivo)
pmd <- getpaired(spmeinvivo)
std <- getstd(pmd)
sda <- getrda(spmeinvivo$mz[std$stdmassindex])
}
\seealso{
\code{\link{getsda}}
}
