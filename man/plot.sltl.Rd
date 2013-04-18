\name{plot.sltl}
\alias{plot.sltl}
\title{
Plot Method for "sltl" object
}
\description{Multiple figure plot with decomposed components as well as the original time series.}
\usage{
\method{plot}{sltl}(x, labels = colnames(X), set.pars = list(mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0), tck = -0.01, mfrow = c(nplot, 1)), main = NULL, ...)
}
\arguments{
  \item{x}{"sltl" object}
  \item{labels}{character with length 4}
  \item{set.pars}{setting for par(.)}
  \item{main}{main title}
  \item{\dots}{optional arguments}
}

\seealso{
\code{\link{sltl}}
}
