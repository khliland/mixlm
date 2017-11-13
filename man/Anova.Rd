\name{Anova.lmm}
\alias{Anova}
\alias{Anova.lmm}
\encoding{latin1}
\title{
Analysis of variance with SS type II or III (including mixed models).
}
\description{
Replacement for \code{Anova.lm} in package \code{car}. This
version adds support for random effects when needed.
}
\usage{
\method{Anova}{lmm}(mod, ...)
}
\arguments{
  \item{mod}{\code{lm}, \code{aov}, \code{glm}, \code{multinom}, \code{polr}
    \code{mlm}, \code{coxph}, \code{lme}, \code{mer}, \code{svyglm} or other suitable model object.}
  \item{\dots}{do not use.}
}
\value{
	Returns appropriate analysis of variance or halts if unsupported input is detected.
}
\author{John Fox \email{jfox@mcmaster.ca}.
Extended by Kristian Hovde Liland.
}

\seealso{
\code{\link[car]{Anova}}, \code{\link{print.AnovaMix}}, \code{\link{AnovaMix}}, \code{\link{lm}}
}
\examples{
dataset <- data.frame(y = rnorm(8),
  x = factor(c(rep(1,4),rep(0,4))),
  z = factor(rep(c(1,0),4)))
mixlm <- lm(y~x*r(z),
  data = dataset)
Anova(mixlm, type="III")
}
\keyword{ ANOVA }
