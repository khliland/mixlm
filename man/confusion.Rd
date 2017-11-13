\name{confusion}
\alias{confusion}
\title{
Confusion matrix.
}
\description{
Computes the confusion matrix of a classification
result.
}
\usage{
confusion(true, predicted)
}
\arguments{
  \item{true}{
	true classes.
}
  \item{predicted}{
	predicted classes.
}
}
\details{
This is a pure print function.
}

\examples{
true <- c('a','a','b','b','c','c')
predicted <- c('a','c','b','b','a','c')
confusion(true, predicted)
}
\keyword{ Classification }
