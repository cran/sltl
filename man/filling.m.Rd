\name{filling.m}
\alias{filling.m}
\title{Fill in missing values}
\description{
Filling in missing values of a monthly timeSeries object
}
\usage{
filling.m(x)
}
\arguments{
  \item{x}{monthly timeSeries object}
}
\details{See McLeod, Hipel and Camacho(1983) for further details.}
\value{monthly timeSeries object}
\references{
McLeod, A.I., Hipel, K.W. and Camacho, F. (1983), 
Trend assessment of water quality time series, Water Resources Bulletin, 
19 537-547.
}
\author{
Hyukjun Gweon and A.I. McLeod
}

\examples{
z <- as.timeSeries(co2)
set.seed(123)

#Removal several observations
index <- sample(1:length(co2),60,replace=FALSE)
z[index] <- NA

#Filling in the missing values
out <- filling.m(z)
}
\keyword{ timeSeries }
