\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{lm}
\alias{lm}
\alias{summary.lmm}
\alias{lmer}
\alias{qr.lmm}
\alias{mixlm}
\alias{print.lmm}
\alias{random.worker}
\concept{regression}
\title{Fitting Linear Models}
\description{
  \code{lm} is used to fit linear models.
  It can be used to carry out regression,
  single stratum analysis of variance and
  analysis of covariance (although \code{\link{aov}} may provide a more
  convenient interface for these).
  The version distributed through the package mixlm
  extends the capabilities with balanced mixture
  models and \code{lmer} interfacing. Random effects 
  are indicated by wrapping their formula entries in \code{r()}.
  Also, effect level names are kept in printing.
}
\usage{
lm(formula, data, subset, weights, na.action,
   method = "qr", model = TRUE, x = TRUE, y = TRUE, qr = TRUE,
   singular.ok = TRUE, contrasts = "contr.sum", offset, 
   unrestricted = TRUE, REML = NULL, equal_baseline = FALSE, \dots)
}
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{lm} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Should be \code{NULL} or a numeric vector.
    If non-NULL, weighted least squares is used with weights
    \code{weights} (that is, minimizing \code{sum(w*e^2)}); otherwise
    ordinary least squares is used.  See also \sQuote{Details},}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \sQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{method}{the method to be used; for fitting, currently only
    \code{method = "qr"} is supported; \code{method = "model.frame"} returns
    the model frame (the same as with \code{model = TRUE}, see below).}

  \item{model, x, y, qr}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned.
  }

  \item{singular.ok}{logical. If \code{FALSE} (the default in S but
    not in \R) a singular fit is an error.}

  \item{contrasts}{character indicating which coding should be applied
    to all \code{factor}s or an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}. Defaults to \code{"contr.sum"}. 
    See Details for more information.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during fitting.
    This should be \code{NULL} or a numeric vector of length equal to
    the number of cases.  One or more \code{\link{offset}} terms can be
    included in the formula instead or as well, and if more than one are
    specified their sum is used.  See \code{\link{model.offset}}.}

  \item{unrestricted}{additional argument for switching between unrestricted and restricted models if including random variables.}

  \item{REML}{is used to invoke restricted maximum likelihood (TRUE)
  or maximum likelihood (FALSE) estimation instead of least squares.}
  
  \item{equal_baseline}{Experimental: Set to \code{TRUE} to let interactions, where a main effect is missing,
    e.g., a nested model, 
    be handled with the same baseline as a cross effect model. If \code{TRUE} the corresponding interactions
    will be put in quotation marks and included in the \code{model.frame}.}

  \item{\dots}{additional arguments to be passed to the low level
    regression fitting functions (see below).}
}
\details{
  Models for \code{lm} are specified symbolically.  A typical model has
  the form \code{response ~ terms} where \code{response} is the (numeric)
  response vector and \code{terms} is a series of terms which specifies a
  linear predictor for \code{response}.  A terms specification of the form
  \code{first + second} indicates all the terms in \code{first} together
  with all the terms in \code{second} with duplicates removed.  A
  specification of the form \code{first:second} indicates the set of
  terms obtained by taking the interactions of all terms in \code{first}
  with all terms in \code{second}.  The specification \code{first*second}
  indicates the \emph{cross} of \code{first} and \code{second}.  This is
  the same as \code{first + second + first:second}.

  If the formula includes an \code{\link{offset}}, this is evaluated and
  subtracted from the response.

  If \code{response} is a matrix a linear model is fitted separately by
  least-squares to each column of the matrix.

  See \code{\link{model.matrix}} for some further details.  The terms in
  the formula will be re-ordered so that main effects come first,
  followed by the interactions, all second-order, all third-order and so
  on: to avoid this pass a \code{terms} object as the formula (see
  \code{\link{aov}} and \code{demo(glm.vr)} for an example).

  A formula has an implied intercept term.  To remove this use either
  \code{y ~ x - 1} or \code{y ~ 0 + x}.  See \code{\link{formula}} for
  more details of allowed formulae.

  Non-\code{NULL} \code{weights} can be used to indicate that different
  observations have different variances (with the values in
  \code{weights} being inversely proportional to the variances); or
  equivalently, when the elements of \code{weights} are positive
  integers \eqn{w_i}, that each response \eqn{y_i} is the mean of
  \eqn{w_i} unit-weight observations (including the case that there are
  \eqn{w_i} observations equal to \eqn{y_i} and the data have been
  summarized).

  \code{lm} calls the lower level functions \code{\link{lm.fit}}, etc,
  see below, for the actual numerical computations.  For programming
  only, you may consider doing likewise.

  All of \code{weights}, \code{subset} and \code{offset} are evaluated
  in the same way as variables in \code{formula}, that is first in
  \code{data} and then in the environment of \code{formula}.
  
  The \code{contrasts} argument is applied when the \code{lm} model is
  fitted. In additional to the standard contrasts from the \code{stats}
  package, one can choose weighted coding: \code{contr.weighted} which
  balances unbalanced data through factor weighting. Different from the
  original version of \code{lm}, a single contrast can be indicated
  to automatically be applied to all \code{factor}s.
}
\value{
  \code{lm} returns an object of \code{\link{class}} \code{c("lmm","lm")} or for
  multiple responses of class \code{c("mlm", "lm")}.

  The functions \code{summary} and \code{\link{anova}} are used to
  obtain and print a summary and analysis of variance table of the
  results.  The generic accessor functions \code{coefficients},
  \code{effects}, \code{fitted.values} and \code{residuals} extract
  various useful features of the value returned by \code{lm}.

  An object of class \code{"lm"} is a list containing at least the
  following components:

  \item{coefficients}{a named vector of coefficients}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted.values}{the fitted mean values.}
  \item{rank}{the numeric rank of the fitted linear model.}
  \item{weights}{(only for weighted fits) the specified weights.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{call}{the matched call.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{contrasts}{(only where relevant) the contrasts used.}
  \item{xlevels}{(only where relevant) a record of the levels of the
    factors used in fitting.}
  \item{offset}{the offset used (missing if none were used).}
  \item{y}{if requested, the response used.}
  \item{x}{if requested, the model matrix used.}
  \item{model}{if requested (the default), the model frame used.}
  \item{na.action}{(where relevant) information returned by
    \code{\link{model.frame}} on the special handling of \code{NA}s.}

  In addition, non-null fits will have components \code{assign},
  \code{effects} and (unless not requested) \code{qr} relating to the linear
  fit, for use by extractor functions such as \code{summary} and
  \code{\link{effects}}.
  
  And models containing random effect will contain
  \code{random} having additional information about
  the model.
}
\section{Using time series}{
  Considerable care is needed when using \code{lm} with time series.

  Unless \code{na.action = NULL}, the time series attributes are
  stripped from the variables before the regression is done.  (This is
  necessary as omitting \code{NA}s would invalidate the time series
  attributes, and if \code{NA}s are omitted in the middle of the series
  the result would no longer be a regular time series.)

  Even if the time series attributes are retained, they are not used to
  line up series, so that the time shift of a lagged or differenced
  regressor would be ignored.  It is good practice to prepare a
  \code{data} argument by \code{\link{ts.intersect}(\dots, dframe = TRUE)},
  then apply a suitable \code{na.action} to that data frame and call
  \code{lm} with \code{na.action = NULL} so that residuals and fitted
  values are time series.
}
\seealso{
  \code{\link{summary.lmm}} for summaries and \code{\link{anova.lmm}} for
  the ANOVA table; \code{\link{aov}} for a different interface.

  The generic functions \code{\link{coef}}, \code{\link{effects}},
  \code{\link{residuals}}, \code{\link{fitted}}, \code{\link{vcov}}.

  \code{\link{predict.lm}} (via \code{\link{predict}}) for prediction,
  including confidence and prediction intervals;
  \code{\link{confint}} for confidence intervals of \emph{parameters}.

  \code{\link{lm.influence}} for regression diagnostics, and
  \code{\link{glm}} for \bold{generalized} linear models.

  The underlying low level functions,
  \code{\link{lm.fit}} for plain, and \code{\link{lm.wfit}} for weighted
  regression fitting.

  More \code{lm()} examples are available e.g., in
  \code{\link{anscombe}}, \code{\link{attitude}}, \code{\link{freeny}},
  \code{\link{LifeCycleSavings}}, \code{\link{longley}},
  \code{\link{stackloss}}, \code{\link{swiss}}.

  \code{biglm} in package \CRANpkg{biglm} for an alternative
  way to fit linear models to large datasets (especially those with many
  cases).
  
  \code{\link{print.AnovaMix}}, \code{\link{AnovaMix}}, \code{\link{Anova}}
}
\references{
  Chambers, J. M. (1992)
  \emph{Linear models.}
  Chapter 4 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

  Wilkinson, G. N. and Rogers, C. E. (1973)
  Symbolic descriptions of factorial models for analysis of variance.
  \emph{Applied Statistics}, \bold{22}, 392--9.
}
\author{
  The design was inspired by the S function of the same name described
  in Chambers (1992).  The implementation of model formula by Ross Ihaka
  was based on Wilkinson & Rogers (1973). Mixed model extensions
  by Kristian Hovde Liland.
}

\note{
  Offsets specified by \code{offset} will not be included in predictions
  by \code{\link{predict.lm}}, whereas those specified by an offset term
  in the formula will be.
}
\examples{
require(graphics)

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
lm.D90 <- lm(weight ~ group - 1) # omitting intercept
anova(lm.D9)
summary(lm.D90)

opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.D9, las = 1)      # Residuals, Fitted, ...
par(opar)

# Linear mixed model
dataset   <- data.frame(y=rnorm(8), x=factor(c(rep(1,4),rep(0,4))), z=factor(rep(c(1,0),4)))
mixlm <- lm(y ~ x*r(z), data = dataset)
Anova(mixlm,type="III")

# Nested model in default version and with equal_baseline = TRUE
anova(lm(y ~ x + x:z, data = dataset))
anova(lm(y ~ x + x:z, data = dataset, equal_baseline=TRUE))

### less simple examples in "See Also" above
}
\keyword{regression}
