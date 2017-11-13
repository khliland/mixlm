\name{print.AnovaMix}
\alias{print.AnovaMix}
\alias{print}
\encoding{latin1}
\title{
	Print method for objects of class(AnovaMix)
}
\description{
	Prints relevant information like the ANOVA table, variance components and errors.
}
\usage{
\method{print}{AnovaMix}(x, ...)
}
\arguments{
  \item{x}{\code{AnovaMix} object to be printed.}
  \item{...}{Additional arguments (not supported yet).}
}
\author{
	Kristian Hovde Liland
}
\note{
	Only balanced models are fully supported.
}

\seealso{
\code{\link{AnovaMix}}, \code{\link{lm}}, \code{\link{Anova}}
}
\examples{
mixlm <- lm(y~x*r(z),
			data = data.frame(y = rnorm(8),
							  x = factor(c(rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),4))))
Anova(mixlm,type="III")
}
\keyword{ ANOVA }
