\name{summary.sltl}
\alias{summary.sltl}
\title{ Summary Method for "sltl" Object }
\description{
summary for "sltl" object.
}
\usage{
\method{summary}{sltl}(object, ...)
}
\arguments{
  \item{object}{ "sltl" object }
  \item{...}{ optional arguments }
}
\value{
A printed summary is given
}
\author{
Hyukjun Gweon and A.I. McLeod
}
\seealso{
\code{\link{sltl}}
}
\examples{
out <- sltl(as.timeSeries(co2),c.window=101,c.degree=1,type="monthly")
summary(out)
}
\keyword{ timeSeries }
