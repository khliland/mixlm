\name{anova.lmm}
\alias{anova.lmm}
\alias{anova}
\encoding{latin1}
\title{
Analysis of variance (sequential SS)
}
\description{
Wrapper for \code{anova.lm} in package \code{stats} that halts execution
if unsupported input is detected.
}
\usage{
\method{anova}{lmm}(object, ...)
}
\arguments{
  \item{object}{object fitted by \code{lm}, \code{lmer} or similar.}
  \item{\dots}{possible additional argument to underlying functions.}
}
\value{
	Returns appropriate analysis of variance or halts if unsupported input is detected.
}
\author{
Kristian Hovde Liland
}

\seealso{
\code{\link{lm}}
}
\examples{
mixlm <- lm(y~x*r(z),
			data = data.frame(y = rnorm(8),
							  x = factor(c(rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),4))))
anova(mixlm)
}
\keyword{ ANOVA }
