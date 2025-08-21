\name{AnovaMix}
\alias{AnovaMix}
\alias{AnovaMixMLM}
\encoding{latin1}
\title{
Mixed model least squares analysis of variance (mixed ANOVA).
}
\description{
Uses output from lm() in mixlm package to compute ANOVA table,
variance components and errors.
}
\usage{
AnovaMix(object, SStype)
AnovaMixMLM(object, SStype)
}
\arguments{
  \item{object}{ object fitted by lm (mixlm package) containing at least one random effect.}
  \item{SStype}{ type of sums-of-squares (I/II/III) for Analysis of Variance.}
}
\details{
\code{AnovaMix} and \code{AnovaMixMLM} can either be invoked directly or through the Anova() function (with type I/II/III error).
}
\value{
  \item{lm }{linear model fitted by lm in package mixlm.}
  \item{anova }{ANOVA table.}
  \item{err.terms }{\code{list} of denominator information for F tests.}
  \item{denom.df }{\code{numeric} of denominator degrees of freedom for F tests.}
  \item{restricted }{\code{logical} indicating if ANOVA used restricted modelling.}
  \item{exp.mean.sq }{\code{character} containing expected mean squares.}
  \item{var.comps }{\code{numeric} containing variance components.}
  \item{random.effects }{\code{character} containing the random effects.}
  \item{ind.randoms }{\code{numeric} with indices of random effects in the model.}
  \item{formula.text }{\code{character} containing all effects of the model.}
}
\author{
	Kristian Hovde Liland
}
\note{
Only balanced models are fully supported.
}

\seealso{
\code{\link{print.AnovaMix}}, \code{\link{Anova}}, \code{\link{lm}}
}
\examples{
# Single response
mydata <- data.frame(y = rnorm(12),
							  x = factor(c(rep(2,4),rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),6)))
mixlm <- lm(y~x*r(z),
			data = mydata)
Anova(mixlm,type="III")

# Double response
mydata <- data.frame(Y = I(as.matrix(cbind(y1=rnorm(12), y2=rnorm(12)))),
							  x = factor(c(rep(2,4),rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),6)))
mixlm <- lm(Y~x*r(z),
			data = mydata)
Anova(mixlm,type="III")
}
\keyword{ ANOVA }
