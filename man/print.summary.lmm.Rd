\name{print.summary.lmm}
\alias{print.summary.lmm}
\title{Summarizing Linear Model Fits}
\usage{
\method{print}{summary.lmm}(x, digits = max(3, getOption("digits") - 3),
      symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{x}{an object of class \code{"summary.lmm"}, usually, a result of a
    call to \code{summary.lmm}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in
    a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{signif.stars}{logical. If \code{TRUE}, \sQuote{significance stars}
    are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
\code{summary} method for class \code{"lmm"}.
}
\details{
  This adaptation of \code{print.summary.lm} from package \code{stats}
  slightly alters the output to better conform with text-book notation.
  
  \code{print.summary.lm} tries to be smart about formatting the
  coefficients, standard errors, etc. and additionally gives
  \sQuote{significance stars} if \code{signif.stars} is \code{TRUE}.

  Correlations are printed to two decimal places (or symbolically): to
  see the actual correlations print \code{summary(object)$correlation}
  directly.
}
\value{
  The function \code{summary.lm} computes and returns a list of summary
  statistics of the fitted linear model given in \code{object}, using
  the components (list elements) \code{"call"} and \code{"terms"}
  from its argument, plus
  \item{residuals}{the \emph{weighted} residuals, the usual residuals
    rescaled by the square root of the weights specified in the call to
    \code{lm}.}
  \item{coefficients}{a \eqn{p \times 4}{p x 4} matrix with columns for
    the estimated coefficient, its standard error, t-statistic and
    corresponding (two-sided) p-value.  Aliased coefficients are omitted.}
  \item{aliased}{named logical vector showing if the original
    coefficients are aliased.} 
  \item{sigma}{the square root of the estimated variance of the random
    error
    \deqn{\hat\sigma^2 = \frac{1}{n-p}\sum_i{w_i R_i^2},}{\sigma^2 = 1/(n-p) Sum(w[i] R[i]^2),}
    where \eqn{R_i}{R[i]} is the \eqn{i}-th residual, \code{residuals[i]}.}
  \item{df}{degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the last
    being the number of non-aliased coefficients.}
  \item{fstatistic}{(for models including non-intercept terms)
    a 3-vector with the value of the F-statistic with
    its numerator and denominator degrees of freedom.}
  \item{r.squared}{\eqn{R^2}, the \sQuote{fraction of variance explained by
    the model},
    \deqn{R^2 = 1 - \frac{\sum_i{R_i^2}}{\sum_i(y_i- y^*)^2},}{R^2 = 1 - Sum(R[i]^2) / Sum((y[i]- y*)^2),}
    where \eqn{y^*}{y*} is the mean of \eqn{y_i}{y[i]} if there is an
    intercept and zero otherwise.}
  \item{adj.r.squared}{the above \eqn{R^2} statistic
    \sQuote{\emph{adjusted}}, penalizing for higher \eqn{p}.}
  \item{cov.unscaled}{a \eqn{p \times p}{p x p} matrix of (unscaled)
    covariances of the \eqn{\hat\beta_j}{coef[j]}, \eqn{j=1, \dots, p}.}
  \item{correlation}{the correlation matrix corresponding to the above
    \code{cov.unscaled}, if \code{correlation = TRUE} is specified.}
  \item{symbolic.cor}{(only if \code{correlation} is true.)  The value
    of the argument \code{symbolic.cor}.}
  \item{na.action}{from \code{object}, if present there.}
}
\seealso{
  The model fitting function \code{\link{lm}}, \code{\link{summary}}.

  Function \code{\link{coef}} will extract the matrix of coefficients
  with standard errors, t-statistics and p-values.
}
\examples{
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
sld90 <- summary(lm.D90 <- lm(weight ~ group -1))# omitting intercept
sld90
}
\keyword{regression}
\keyword{models}
