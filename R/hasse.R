#' Hasse Diagram from Linear Model
#'
#' @param object A linear model object, e.g., \code{\link{lm}}.
#' @param manualTerms A vector of terms that should be placed manually in the diagram.
#' @param manualParents A list of vectors with the parents of the terms in \code{manualTerms}.
#' @param meanName The name of the mean term (default = "M").
#' @param errorName The name of the error term (default = "(E)").
#' 
#' @description
#' This function extracts terms from a linear model object and creates a Hasse 
#' diagram of the terms. The function is useful for visualizing the structure of 
#' a linear model. If the model contains random effects, these are placed in
#' parentheses in the diagram. Manually placed terms are supported to some
#' extent, see example for usage.
#' 
#' @details
#' Plotting is handled by the \code{hasse} function from the \code{hasseDiagram}
#' package, requiring also the \code{Rgraphviz} package.
#'
#' @return A list with the levels of the diagram and the adjacency matrix.
#' @export
#'
#' @examples
#' # Random data
#' dat <- data.frame(A = factor(rep(c(1,2),32)), 
#'                   B = factor(rep(c(1,1,2,2),16)), 
#'                   C = factor(rep(c(1,1,1,1,2,2,2,2),8)), 
#'                   Ind = factor(rep(c(1,1,1,1,2,2,2,2),8)),
#'                   D = factor(rep(c(rep(1,8),rep(2,8)),4)), 
#'                   E = factor(rep(c(rep(1,16),rep(2,16)),2)))
#' dat$y = rnorm(64)
#' 
#' # Linear model with interactions and nested factors
#' mod <- lm(y~A*B*C + D + E%in%D, data=dat)
#' (an <- Anova(mod, type="II"))
#' H <- hasseMod(mod)
#' \notrun{ # Requires installation of Rgraphviz and hasseDiagram
#' library(Rgraphviz)
#' hasse(H$hasse, parameters=list(cluster = FALSE, arrows = "none", edgeColor = "darkred"))
#' }
#' 
#' # Linear model with repeated measures where Ind is nested in A
#' modv <- lm(y~A*r(B) + r(Ind), data=dat)
#' (anv <- Anova(mod, type="II"))
#' Hv <- hasseMod(modv, manualTerms=c("Ind"), manualParents=list(c("A")))
#' \notrun{ # Requires installation of Rgraphviz and hasseDiagram
#' hasse(Hv$hasse, parameters=list(cluster = FALSE, arrows = "none", edgeColor = "darkred"))
#' }
hasseMod <- function(object, manualTerms=NULL, manualParents=NULL, 
                     meanName="M", errorName="(E)"){
  tt <- attr(terms(object),"factors")
  terms <- colnames(tt)
  tn <- c(meanName, terms, errorName)
  M <- matrix(FALSE, length(tn), length(tn))
  dimnames(M) <- list(tn, tn)
  
  ## Add standard contents
  levs <- list()
  levs[[1]] <- meanName
  levs[[2]] <- character(0)
  depth <- 2
  
  ## Loop over terms to set right level
  for(i in 1:length(terms)){
    # Check if special term handling is needed
    if(!(terms[i] %in% manualTerms)){
      splat <- strsplit(terms[i], ":")[[1]]
      # Add level if needed
      if(length(splat) >= depth){
        depth <- length(splat)+1
        levs[[depth]] <- character(0)
      }
      # Add term in right level
      levs[[depth]] <- c(levs[[depth]], terms[i])
    }
  }
  ## Special terms
  for(i in 1:length(terms)){
    # Check if special term handling is needed
    if(terms[i] %in% manualTerms){
      curLev <- 1
      for(j in 2:depth){
        for(k in 1:length(manualParents)){
          if(manualParents[k] %in% levs[[j]]){
            curLev <- j+1
          }
        }
      }
      if(curLev > depth){
        depth <- curLev
        levs[[depth]] <- character(0)
      }
      # Add term in right level
      levs[[curLev]] <- c(levs[[curLev]], terms[i])
    }
  }
  levs[[depth+1]] <- errorName
  
  ## Connect terms
  # Mean to main effects
  for(i in levs[[2]]){
    M[1,i] <- TRUE
  }
  # Loop over remaining levels
  for(i in 3:depth){
    for(j in 1:length(levs[[i]])){
      splatj <- strsplit(levs[[i]][j], ":")[[1]]
      for(k in 1:length(levs[[i-1]])){
        splatk <- strsplit(levs[[i-1]][k], ":")[[1]]
        if(all(splatk %in% splatj))
          M[levs[[i-1]][k], levs[[i]][j]] <- TRUE
      }
    }
  }
  # Add special connections
  for(i in 1:length(manualTerms)){
    for(j in 1:length(manualParents[[i]])){
      M[manualParents[[i]][j], manualTerms[i]] <- TRUE
    }
  }
  # Connect loose ends
  for(i in 2:(length(terms)+1)){
    if(sum(M[i,])==0)
      M[i,errorName] <- TRUE
  }
  # Add parentheses around random effects
  if(!is.null(object$random)){
    for(r in object$random$random){
      tn[tn==r] <- paste("(", r, ")", sep="")
    }
    dimnames(M) <- list(tn, tn)
  }
  ret <- list(levels=levs, hasse=M)
  class(ret) <- "hasseMod"
  return(ret)
}

# datt <- data.frame(A = factor(rep(c(1,2),32)), 
#                    B = factor(rep(c(1,1,2,2),16)), 
#                    C = factor(rep(c(1,1,1,1,2,2,2,2),8)), 
#                    Ind = factor(rep(c(1,1,1,1,2,2,2,2),8)),
#                    D = factor(rep(c(rep(1,8),rep(2,8)),4)), 
#                    E = factor(rep(c(rep(1,16),rep(2,16)),2)))
# datt$y = rnorm(64)
# mod <- lm(y~A*B*C + D + E%in%D, data=datt)
# (an <- Anova(mod, type="II"))
# H <- hasseMod(mod)
# hasse(H$hasse, parameters=list(cluster = FALSE, arrows = "none", edgeColor = "darkred"))
# 
# datv <- data.frame(A = factor(rep(c(1,2),32)), 
#                    B = factor(rep(c(1,1,2,2),16)), 
#                    Ind = factor(rep(c(1,1,1,1,2,2,2,2),8)))
# datv$y = rnorm(64)
# modv <- lm(y~A*r(B) + r(Ind), data=datv)
# (anv <- Anova(mod, type="II"))
# Hv <- hasseMod(modv, manualTerms=c("Ind"), manualParents=list(c("A")))
# hasse(Hv$hasse, parameters=list(cluster = FALSE, arrows = "none", edgeColor = "darkred"))

