\name{contr.weighted}
\alias{contr.weighted}
\title{
Contrast matrix for weighted effect coding
}
\description{
Weighted contrast coding for linear models.
}
\usage{
contr.weighted(x, base)
}
\arguments{
  \item{x}{ \code{factor} for which a contrast matrix should be made.}
  \item{base}{ \code{factor} level used as basis for contrast coding.
  Default is the (first) level with maximum frequency.}
}
\details{
Different from the contrasts made throught the \code{stats} package functions
this contrast requires a full \code{factor} vector as input rather than its
respective \code{levels} as weights are computed from the frequencies of the
\code{factor} \code{levels}.
}
\value{
A \code{matrix} with \code{n} rows and \code{n-1} values.
}
\references{
Nieuwenhuis, R.; Grotenhuis, M.; Pelzer, B. Weighted Effect Coding for Observational Data with wec. R. J. 2017, 9, 477–485.
}
\author{
Kristian Hovde Liland
}
\note{
\code{contr.weighted} cannot be used directly as a replacement for other
contrasts by name, but must be used via contrasts matrix computations.
}

\seealso{
\code{\link{lm}}
}
\examples{
balanced   <- factor(c(rep("A", 3), rep("B", 3), rep("C", 3)))
unbalanced <- factor(c(rep("A", 3), rep("B", 3), rep("C", 2)))
# Weighted coding when applied to balanced data
contr.weighted(balanced)
# Weighted coding when applied to unbalanced data (default base level)
contr.weighted(unbalanced)
# Weighted coding when applied to unbalanced data (base level = "C")
contr.weighted(unbalanced, "C")
}
