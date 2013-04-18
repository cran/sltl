\name{sltl}
\alias{sltl}
\title{
Time Series Decomposition using Loess and Harmonic Regression
}
\description{
Decompose a timeSeries object into seasonal, long-term and remainder components.
Loess is used for long-term components while harmonic regression is employed for seasonality.
The function can handle multiple long-term components cases. 
}
\usage{
sltl(x, c.window, c.degree = 1, Fm = 6, Fs = 0, searchQ = FALSE, ic = c("BIC", "AIC"), type = c("daily", "monthly"), filling = FALSE, npass = 1)
}
\arguments{
  \item{x}{daily or Montly timeSeries object}
  \item{c.window}{span of the loess for long-term extraction. Can be a vector.}
  \item{c.degree}{degree of local polynomial in long-term extraction. Must have the same demension as c.window.}
  \item{Fm}{number of Fourier components for seasonal mean. If Fm=0, then only the overall mean of series is used.}
  \item{Fs}{number of harmonics for seasonal standard deviations.If Fs=0, only overall standard deviation is used.}
  \item{searchQ}{IF TRUE, search for best BIC/AIC model using harmonic reqressions of maximum orders Fm and Fs
             for seasonal means and standard deviations.}
  \item{ic}{
"BIC" or "AIC" model selection
}
  \item{type}{type of timeSeries object. Must be 'daily' or 'monthly'}
  \item{filling}{IF TRUE, missing values are filled in by the seasonal adjustment algorithm prior to decomposition process.}
  \item{npass}{the number of passes in backfitting process. Should be less than 5.}
}
\details{See McLeod and Gweon(2013) for further details.}
\value{an object of class "sltl" with components}
\references{
McLeod, A. I. and  Hyukjun Gweon.  (2013). 
Optimal Deseasonalization for Geophysical Time Series. 
Journal of Environmental Statistics. Vol. 4, Issue 11
}
\author{Hyukjun Gweon and A.I. McLeod}

\examples{

#Example 1: With monthly co2 data
out <- sltl(as.timeSeries(co2),c.window=101,c.degree=1,type="monthly")
plot(out)

#Example 2: With additional long term oscillation 
out2 <- sltl(as.timeSeries(co2),c.window=c(101,35),c.degree=c(1,2),type="monthly")
plot(out2)
     }
