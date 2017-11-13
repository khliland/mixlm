\name{fparse}
\alias{fparse}
\encoding{latin1}
\title{
Effects of formulas.
}
\description{
Extracts all effects from a formula, even though inside functions or interactions.
}
\usage{
fparse(f)
}
\arguments{
  \item{f}{\code{formula} to be parsed.}
}
\value{
	Returns a \code{character} vector containing all effects.
}
\author{
Bjørn-Helge Mevik
}

\seealso{
\code{\link{rparse}}
}
\examples{
f <- formula(y~x*r(z))
fparse(f)
}
