\name{best.subsets}
\alias{best.subsets}
\title{
F-test based best subset selection.
}
\description{
Adaptation of existing methods based on AIC/BIC.
}
\usage{
best.subsets(model, nbest = 5, nvmax, digits, force.in = "NULL")
}
\arguments{
  \item{model}{
	object class \code{lm} to select effects from.
}
  \item{nbest}{
	\code{numeric} indicating number of models to report of each size.
}
  \item{nvmax}{
	\code{numeric} maximum size of subsets to examine.
}
  \item{digits}{
	\code{numeric} giving number of digits in format of output.
}
  \item{force.in}{
	\code{character} vector indicating effects to keep in all models.
}
}
\details{
F-based versions of built in subset method.
}
\value{
No return, only print.
}
\author{
Kristian Hovde Liland
}

\examples{
data <- data.frame(y = rnorm(8),
				   x = factor(c('a','a','a','a','b','b','b','b')),
				   z = factor(c('a','a','b','b','a','a','b','b')))
mod <- lm(y ~ x + z, data=data)
best.subsets(mod)
}
\keyword{ regression }
