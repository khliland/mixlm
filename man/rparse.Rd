\name{rparse}
\alias{rparse}
\encoding{latin1}
\title{
	Removes function r() from formulas.
}
\description{
	Removes function r() from formulas. Can also convert to lmer formula.
}
\usage{
rparse(f, REML = FALSE)
}
\arguments{
  \item{f}{formula to be stripped of r().}
  \item{REML}{logical indicating if lmer conversion should be done.}
}
\value{
	Formula without r(), possibly converted to lmer mixed model format.
}
\author{
	Kristian Hovde Liland
}

\seealso{
\code{\link{fparse}}
}
\examples{
f <- formula(y~x*r(z))
rparse(f)
}