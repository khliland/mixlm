\name{is.balanced}
\alias{is.balanced}
\encoding{latin1}
\title{
Balance cheking of models.
}
\description{
Checks if models have balanced data.
}
\usage{
is.balanced(object)
}
\arguments{
  \item{object}{fitted model that includes \code{variables} attribute and \code{model} slot.}
}
\value{
	Returns TRUE if balanced, FALSE if not.
}
\author{
Kristian Hovde Liland
}

\examples{
mixlm <- lm(y~x*r(z),
		    data = data.frame(y = rnorm(8),
							  x = factor(c(rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),4))))
is.balanced(mixlm)
}
