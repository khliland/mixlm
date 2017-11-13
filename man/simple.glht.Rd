\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{simple.glht}
\alias{simple.glht}
\alias{print.simple.glht}
\alias{cld}
\title{
Pairwise comparison with multiple testing compensation.
}
\description{
  Extension of glht from the multcomp package to handle
  Fisher family-wise error and Bonferroni testing.
  Create a set of confidence intervals on the differences between the
  means of the levels of a factor with the specified family-wise
  probability of coverage.  The intervals are based on the Studentized
  range statistic, Tukey's \sQuote{Honest Significant Difference}
  method, Fisher's family-wise error, or Bonferroni testing.

}
\usage{
simple.glht(mod, effect, corr = c("Tukey","Bonferroni","Fisher"),
	level = 0.95, df = NULL, \dots)
}
\arguments{
 \item{mod}{A fitted model object, usually an \code{\link{lm}} or \code{\link{glm}} fit.}
 \item{effect}{A character vector giving the term of the fitted model for
   which the intervals should be calculated. This can also be an interaction.}
  \item{corr}{A character vector giving the multiple testing correction
  method. Defaults to \code{Tukey}.}
 \item{level}{A numeric value between zero and one giving the
   family-wise confidence level to use.}
 \item{df}{User supplied number of degrees of freedom. If not supplied or
   NULL, the default is to extract these from the model.}
 \item{\dots}{Optional additional arguments.  None are used at present.}
}
\details{
  When comparing the means for the levels of a factor in an analysis of
  variance, a simple comparison using t-tests will inflate the
  probability of declaring a significant difference when it is not in
  fact present.  This because the intervals are calculated with a
  given coverage probability for each interval but the interpretation of
  the coverage is usually with respect to the entire family of
  intervals.

  John Tukey introduced intervals based on the range of the
  sample means rather than the individual differences.  The intervals
  returned by this function are based on this Studentized range
  statistics.

  The intervals constructed in this way would only apply exactly to
  balanced designs where there are the same number of observations made
  at each level of the factor.  This function incorporates an adjustment
  for sample size that produces sensible intervals for mildly unbalanced
  designs.

  If \code{which} specifies non-factor terms these will be dropped with
  a warning: if no terms are left this is a an error.
}
\value{
  An object of classes \code{"simple.glht"}, \code{"summary.glht"}
  and \code{"glht"} containing information to produce confidence intervals,
  tests and plotting.

  There are \code{print}, \code{plot} and \code{cld} methods for class
  \code{"simple.glht"}.
  The \code{plot} method does not accept
  \code{xlab}, \code{ylab} or \code{main} arguments and creates its own
  values for each plot.
}
\references{
  Miller, R. G. (1981)
  \emph{Simultaneous Statistical Inference}. Springer.

  Yandell, B. S. (1997)
  \emph{Practical Data Analysis for Designed Experiments}.
  Chapman & Hall.
}
\author{
  Douglas Bates, extended to mixed effect models by Kristian Hovde Liland.
}
\seealso{
  \code{\link{aov}}, \code{\link{qtukey}}, \code{\link{model.tables}},
  \code{\link[multcomp]{glht}} in package \CRANpkg{multcomp}.
}
\examples{
require(graphics)

summary(fm1 <- lm(breaks ~ wool + tension, data = warpbreaks))
simple.glht(fm1, "tension")
plot(simple.glht(fm1, "tension"))
cld(simple.glht(fm1, "tension"))
}
\keyword{models}
\keyword{design}
