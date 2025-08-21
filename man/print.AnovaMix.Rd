\name{print.AnovaMix}
\alias{print.AnovaMix}
\alias{print.AnovaMixMLM}
\alias{print}
\encoding{latin1}
\title{
	Print method for objects of class(AnovaMix) and class(AnovaMixMLM)
}
\description{
	Prints relevant information like the ANOVA table, variance components and errors.
}
\usage{
\method{print}{AnovaMix}(x, ...)
\method{print}{AnovaMixMLM}(x, var, ...)
}
\arguments{
  \item{x}{\code{AnovaMix(MLM)} object to be printed.}
  \item{var}{which response variable to print.}
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
# Single response (AnovaMix)
mixlm <- lm(y~x*r(z),
			data = data.frame(y = rnorm(8),
							  x = factor(c(rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),4))))
Anova(mixlm,type="III")

# Double response (AnovaMixMLM)
mydata <- data.frame(Y = I(as.matrix(cbind(y1=rnorm(12), y2=rnorm(12)))),
							  x = factor(c(rep(2,4),rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),6)))
mixlm <- lm(Y~x*r(z),
			data = mydata)
Anova(mixlm,type="III")
}
\keyword{ ANOVA }
