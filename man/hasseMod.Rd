\name{hasseMod}
\alias{hasseMod}
\title{
Hasse Diagram from Linear Model
}
\description{
This function extracts terms from a linear model object and creates a Hasse 
diagram of the terms. The function is useful for visualizing the structure of 
a linear model. If the model contains random effects, these are placed in
parentheses in the diagram. Manually placed terms are supported to some
extent, see example for usage.
}
\usage{
hasseMod(object, manualTerms=NULL, manualParents=NULL, 
         meanName="M", errorName="(E)")
}
\arguments{
  \item{object}{A linear model object, e.g., \code{\link{lm}}.}
  \item{manualTerms}{A \code{vector} of terms that should be placed manually in the diagram.}
  \item{manualParents}{A \code{list} of vectors with the parents of the terms in \code{manualTerms}.}
  \item{meanName}{The name of the mean term (default = "M").}
  \item{errorName}{The name of the error term (default = "(E)").}
}
\value{
	A list with the levels of the diagram and the adjacency matrix.
}
\author{
Kristian Hovde Liland
}

\seealso{
\code{\link{lm}}
}
\examples{
# Random data
dat <- data.frame(A = factor(rep(c(1,2),32)), 
                  B = factor(rep(c(1,1,2,2),16)), 
                  C = factor(rep(c(1,1,1,1,2,2,2,2),8)), 
                  Ind = factor(rep(c(1,1,1,1,2,2,2,2),8)),
                  D = factor(rep(c(rep(1,8),rep(2,8)),4)), 
                  E = factor(rep(c(rep(1,16),rep(2,16)),2)))
dat$y = rnorm(64)

# Linear model with interactions and nested factors
mod <- lm(y~A*B*C + D + E\%in\%D, data=dat)
(an <- Anova(mod, type="II"))
H <- hasseMod(mod)
\dontrun{ # Requires installation of Rgraphviz
library(Rgraphviz)
hasse(H$hasse, parameters=list(cluster = FALSE, arrows = "none", edgeColor = "darkred"))
}

# Linear model with repeated measures where Ind is nested in A
modv <- lm(y~A*r(B) + r(Ind), data=dat)
(anv <- Anova(mod, type="II"))
Hv <- hasseMod(modv, manualTerms=c("Ind"), manualParents=list(c("A")))
\dontrun{ # Requires installation og Rgraphviz
hasse(Hv$hasse, parameters=list(cluster = FALSE, arrows = "none", edgeColor = "darkred"))
}
}
