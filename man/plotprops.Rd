\name{plotprops}
\alias{plotprops}
\title{
Property plots for relevant component analysis
}
\description{
Plot summary of relevant component analysis.
}
\usage{
plotprops(Y, X, doscaleX = FALSE, docenterX = TRUE, ncomp, subset)
}
\arguments{
  \item{Y}{
	Response matrix.
}
  \item{X}{
	Predictor matrix.
}
  \item{doscaleX}{
	Standardize predictors.
}
  \item{docenterX}{
	Center predictors.
}
  \item{ncomp}{
	Number of singular values and eigenvalues to extract.
}
  \item{subset}{
	Subset of predictor and response.
}
}
\value{
Only plotting.
}
\references{
Helland, I.S. & T. Almøy (1994) Comparison of prediction methods when only a few components are relevant. JASA 89, 583-591.
}
\author{
Solve Sæbø
}

\examples{
X <- matrix(rnorm(100),20,5)
Y <- matrix(rnorm(20),20,1)
plotprops(Y, X, doscaleX = FALSE, docenterX = TRUE, 5)
}
\keyword{ Regression }
\keyword{ Classification }
