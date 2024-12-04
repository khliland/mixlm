\name{effect.labels}
\alias{effect.labels}
\title{
Create new effect labels for \code{lm}
}
\description{
Alternative notation of effect labels including levels.
}
\usage{
effect.labels(t, data, contrasts)
}
\arguments{
  \item{t}{
	Terms object.
}
  \item{data}{
	Corresponding model.matrix.
}
  \item{contrasts}{
	List of contrasts, similar to the \code{contrast} argument of \code{lm}.
}
}
\value{
  \item{names }{Character vector of effect labels.}
}
\author{
Kristian Hovde Liland
}

\keyword{ effect }
\keyword{ label }
