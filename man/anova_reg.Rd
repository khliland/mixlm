\name{anova_reg}
\alias{anova_reg}
\title{
Analysis of variance for regression.
}
\description{
Summarizes all effects in one.
}
\usage{
anova_reg(lm.object)
}
\arguments{
  \item{lm.object}{
	an object of class \code{lm}.
}
}
\value{
Returns data.frame containing analysis of variance
}
\author{
Kristian Hovde Liland
}
\examples{
anova_reg(lm(y~x, data=data.frame(y=1:4,x=rnorm(4))))
}
