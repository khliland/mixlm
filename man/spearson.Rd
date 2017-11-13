\name{spearson}
\alias{spearson}
\title{
Standardized Pearson residuals
}
\description{
Standardized Pearson residuals.
}
\usage{
spearson(object)
}
\arguments{
  \item{object}{
	fitted model.
}
}
\details{
Takes ordinary Pearson residuals and standardizes them.
}
\value{
Returns the residuals.
}
\author{
Kristian Hovde Liland
}

\examples{
data <- data.frame(y = rnorm(8),
				   x = factor(c('a','a','a','a','b','b','b','b')),
				   z = factor(c('a','a','b','b','a','a','b','b')))
mod <- lm(y ~ x + z, data=data)
spearson(mod)
}

\keyword{ Regression }
