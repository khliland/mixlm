\name{CIgrandMean}
\alias{CIgrandMean}
\alias{print.CIgm}
\title{
Confidence interval for the grand mean of a linear model
}
\description{
This function estimates the confidence interval for 
the grand mean of a balanced linear (mixed) model.
}
\usage{
CIgrandMean(object, alpha = 0.05)
\method{print}{CIgm}(x, ...)
}
\arguments{
  \item{object}{
	An \code{lm} object possibly containing random effects.
}
  \item{alpha}{
	A scalar significance level for the confidence interval.
}
  \item{x}{
	An object returned from \code{CIgrandMean}.
}
  \item{\dots}{
	Additional arguments (not used).
}
}
\details{
This implementation is only valid for models containing no
continuous effects and only balanced data.
}
\value{
\code{CIgrandMean} returns a vector of interval endpoints and center.
\code{print.CIgm} has no return.
}
\references{
Suggestions are welcome.
}
\author{
Kristian Hovde Liland
}

\examples{
set.seed(42)
dataset   <- data.frame(y=rnorm(8), x=factor(c(rep(1,4),rep(0,4))), z=factor(rep(c(1,0),4)))
mixlm <- lm(y~x*r(z), data = dataset)
CIgrandMean(mixlm)
}

\keyword{ Confidence interval }
\keyword{ Linear model }
