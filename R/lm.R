##
## This file contains MODIFIED COPIES of the lm, summary.lm and confint.lm methods from the base package stats (2014-10-31)##

# The lm function replaces stats::lm, organizes fixed and random effects, removes r() from formula and parses to lm or lmer.

# lm from stats, edited to use treatment names in sum contrasts
# and enable limited classical least squares mixed models
if(requireNamespace("lme4", quietly = TRUE)){
  lmer <- lme4::lmer
} else {
  lmer <- function()warning("Install package lme4 to enable REML/ML modelling.")
}
lm <- function (formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = TRUE, y = TRUE,
                qr = TRUE, singular.ok = TRUE, contrasts = "contr.sum",
                offset, unrestricted = TRUE, REML = NULL, equal_baseline=FALSE, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  ## Edited by KHL
  mfd <- match(c("formula","data"), names(mf), 0L)
  if(length(mfd)==2){ # Has formula and data
    is.random <- TRUE
    if( any(grepl("r(",formula,fixed=TRUE)) ){
      rw <- random.worker(formula, data, REML)
    } else {
      rw <- list(0)
    }
    if(length(rw) == 1){
      is.random <- FALSE
    } else { # Removed r() from formula
      formula <- rw$formula
      mf$formula <- rw$formula
      rw$unrestricted <- unrestricted
    }
  } else {
    is.random <- FALSE
  }
  ## End of edit
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if (method != "model.frame")
    if (method != "qr")
      warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
            domain = NA)
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  # Create contrast list if single character argument is supplied to contrasts (edit by KHL)
  contrasts.orig <- contrasts
  if(!is.null(contrasts)){
    if(is.character(contrasts) && length(contrasts)==1){
      # Handle contrasts given as a single string
      facs <- which(unlist(lapply(mf, inherits, what = "factor")))
      contrasts <- as.list(rep(contrasts, length(facs)))
      names(contrasts) <- names(facs)
    }
    contrasts.names <- contrasts
    # Handle contrasts given as lists
    if(is.list(contrasts) && !is.logical(REML)){# && length(contrasts) > 1){
      facs <- which(unlist(lapply(mf, inherits, what = "factor")))
      if(length(facs) != length(contrasts))
        stop("Number of contrasts must match number of factors when specified separately")
      if(length(facs)>0)
        for(i in 1:length(contrasts)){
          if(!is.matrix(contrasts[[i]])){
            nl <- nlevels(mf[[names(contrasts)[i]]])
            if(contrasts[[i]] == "contr.treatment.last"){ # Force last level of factor to base level
              contrasts[[i]] <- contr.treatment(nl,nl)
            } else
              if(contrasts[[i]] == "contr.weighted")
                contrasts[[i]] <- contr.weighted(mf[[names(contrasts)[i]]])
              else {
                if(contrasts[[i]] == "contr.treatment")
                  contrasts[[i]] <- contr.treatment(nl)
                else {
                  if(contrasts[[i]] == "contr.sum")
                    contrasts[[i]] <- contr.sum(nl)
                  else {
                    if(contrasts[[i]] == "contr.poly")
                      contrasts[[i]] <- contr.poly(nl)
                    else {
                      if(contrasts[[i]] == "contr.SAS")
                        contrasts[[i]] <- contr.SAS(nl)
                      else {
                        if(contrasts[[i]] == "contr.helmert")
                          contrasts[[i]] <- contr.helmert(nl)
                      }
                    }
                  }
                }
              }
          }
        }
      else
        contrasts <- NULL
    }
  }
  # lme4 handling
  if(length(rw) != 1){
    if(is.logical(REML)){ # Perform 
      if(requireNamespace("lme4", quietly = TRUE)){
        cl[[1]] <- as.name("lmer")
        cl[["contrasts"]] <- contrasts[rw$fixed[!grepl(":",rw$fixed)]]
        cl[["formula"]] <- rw$reml.formula
        object <- eval(cl,parent.frame())
        object@call <- cl
        return(object)
      } else {
        warning('Package lme4 required for random REML/ML models.')
      }
    }
  }
  ## Expand model.frame and adapt formula if necessary for missing main effects
  if(equal_baseline){
    eformula <- .extend_formula(formula)
    # Check that there are interactions and all included interaction variables are factors
    if(!is.null(eformula$interactions) && !length(eformula$interactions)==0 && all(unlist(lapply(mf[eformula$variables], inherits, what="factor")))){
      # Expand mf
      for(i in 1:length(eformula$interactions)){
        mf[ncol(mf)+1] <- interaction(mf[strsplit(eformula$interactions[i],":")[[1]]])
      }
      colnames(mf)[(ncol(mf)-length(eformula$interactions)+1):ncol(mf)] <- eformula$interactions
      # Adapt formula
      for(i in 1:length(eformula$interactions)){
        formula <- formula(paste(gsub(eformula$interactions[[i]], paste0("`",eformula$interactions[[i]],"`"),as.character(formula))[c(2,1,3)], collapse=" "))
      }
      mt <- terms(formula)
      for(i in 1:length(eformula$interactions)){
        contrasts[[length(contrasts)+1]] <- .krons(contrasts[strsplit(eformula$interactions[[i]],":")[[1]]])
        names(contrasts)[length(contrasts)] <- eformula$interactions[[i]]
      }
      contrasts <- contrasts[setdiff(names(contrasts), eformula$missing)]
    } else {
      stop("'equal_baseline' should only be used with interactions where one of its factors are missing from the formula, and all interaction-variables should be factors.")
    }
  }
  if (method == "model.frame")
    return(mf)
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  
  # ccs <- FALSE
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y))
      matrix(,0,3) else numeric(), residuals = y,
      fitted.values = 0 * y, weights = w, rank = 0L,
      df.residual = if(!is.null(w)) sum(w != 0) else
        if (is.matrix(y)) nrow(y) else length(y))
    if(!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    # Alternative handling of missing main effects in model.matrix
    x <- model.matrix(object=mt, data=mf, contrasts.arg=contrasts)
    effect.sources <- effect.source(mt,mf)
    ## Edited by KHL (CCS = Cell Count Scaling)
    col.names   <- effect.labels(mt,mf,contrasts.names) # mt is "terms" from formula, x is model.matrix
    if(length(col.names)==length(colnames(x))){
      colnames(x) <- col.names
    }
    if((is.null(contrasts.orig) && options("contrasts")[[1]][1] %in% c("contr.sum", "contr.weighted")) ||
       (is.list(contrasts.orig) && all(unlist(contrasts.orig) %in% c("contr.sum", "contr.weighted"))) ||
       (is.character(contrasts.orig) && contrasts.orig %in% c("contr.sum", "contr.weighted"))){
      # Special handling of interactions for weighted coding
      if((is.null(contrasts.orig) && options("contrasts")[[1]][1] %in% c("contr.weighted")) ||
         (is.list(contrasts.orig) && all(unlist(contrasts.orig) %in% c("contr.weighted"))) ||
         (is.character(contrasts.orig) && contrasts.orig %in% c("contr.weighted"))){      
        mt_factors <- attr(mt, "factors")
        main_interactions <- colSums(mt_factors)
        if(any(main_interactions>1)){
          # Use contr.sum as basis for weighted interactions
          contsum <- as.list(rep("contr.sum", length(facs)))
          names(contsum) <- names(facs)
          x_sum <- model.matrix(mt, mf, contsum)
          ass <- attr(x_sum, "assign")
          
          for(i in which(main_interactions>1)){
            # Convert columns of model.matrix to factor and use to find weights
            int_fac <- interaction(mf[rownames(mt_factors)[mt_factors[,i]==1]])
            n_each  <- table(int_fac)
            if(any(n_each==0))
              warning(paste0("Contrast error due to empty cell"))
            x_col <- which(ass==i)
            for(lev in levels(int_fac)){
              x[int_fac == lev, x_col] = x_sum[int_fac == lev, x_col] * min(n_each)/n_each[lev]
            }
          }
        }
      }
    }
    ## End edit
    z <- if(is.null(w)) lm.fit(x, y, offset = offset,
                               singular.ok=singular.ok, ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok=singular.ok, ...)
  }
  if(is.matrix(y)){
    class(z) <- c( "lmm","mlm", "lm")
  } else {
    class(z) <- c( "lmm", "lm")
  }
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$contrasts.names <- contrasts.names
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr) z$qr <- NULL
  ## Edited by KHL
  if( is.random ){
    z$random <- rw
    if(!all(grepl("factor",attr(mt,"dataClasses")[-1])|grepl("ordered",attr(mt,"dataClasses")[-1]))){
      stop("Mixed models containing continuous effects not supported")
    }
  }
  if(exists("effect.sources") && !is.null(effect.sources))
    z$effect.sources <- effect.sources
  #  if(ccs) # Save contr.sum_ccs weights
  #    z$ccs <- wgt
  ## End edit
  z
}

## Collect and extract randomness
random.worker <- function(formula, data, REML = NULL){
  formula <- formula(formula)
  terms <- terms(formula)
  effsr <- attr(terms,"term.labels")
  effs  <- attr(terms(rparse(formula)),"term.labels")
  if(length(effs)==0){
    return( list(0) )
  }
  
  has.intercept <- attr(terms,"intercept")==1
  rands <- sort(unique(c(grep("[:]r[(]",effsr),   # Match random in interaction
                         grep("^r[(]",  effsr),   # Match random in the beginning
                         grep("[(]r[(]",effsr)))) # Match random inside function
  
  # which.rands <- match(rands,effsr)
  eff.splits <- list()
  for(i in 1:length(effs)){ # Split effect to look for hidden random interactions
    eff.splits[[i]] <- fparse(formula(paste("1~", effs[i],sep="")))
  }
  eff.lengths <- lapply(eff.splits,length)
  main.effs   <- effs[eff.lengths==1]
  main.rands  <- main.effs[main.effs%in%effs[rands]]
  main.rands.only.inter <- character(0)
  for(i in rands){
    main.rands.only.inter <- c(main.rands.only.inter, setdiff(eff.splits[[i]],main.effs)) # Random main effects only present in interactions
  }
  inter.rands <- which(unlist(lapply(eff.splits,function(i) any(main.rands%in%i))))
  # Check if any interactions containing random effects are not labeled as random
  if(any(is.na(match(inter.rands,rands)))){
    extra.randoms <- inter.rands[which(is.na(match(inter.rands,rands)))]
    warning(paste(paste(effs[extra.randoms],sep="",collapse=", "), " included as random interaction",ifelse(length(extra.randoms)==1,"","s"),sep=""))
    rands <- cbind(rands,extra.randoms)
    effs  <- effs[!(extra.randoms%in%effs)]
  }
  if(length(rands)==0){
    return( list(0) ) 
  } else {
    if(is.logical(REML)){
      remleffs     <- c(effs[setdiff(1:length(effs),rands)],paste("(1|",effs[rands],")",sep=""))
      reml.formula <- formula(paste(formula[[2]],"~",paste(remleffs,collapse="+"),ifelse(has.intercept,"","-1"),sep=""))
      
      return( list(rformula = formula, formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")), random = effs[rands], main.rands.only.inter = main.rands.only.inter, fixed = effs[setdiff(1:length(effsr),rands)], all = effs, allr = effsr, has.intercept = has.intercept, remleffs = remleffs, reml.formula = reml.formula))
    } else {
      return( list(rformula = formula, formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")), random = effs[rands], main.rands.only.inter = main.rands.only.inter, fixed = effs[setdiff(1:length(effsr),rands)], all = effs, allr = effsr, has.intercept = has.intercept))
    }
  }
}


###########################################
# summary.lm from stats, edited to enable limited classical least squares mixed models
summary.lmm <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
                                               "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- qr.lmm(object)
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  p1 <- 1L:p # Moved up for the sake of "random"
  if(is.null(object$random)){
    resvar <- rss/rdf
  } else {
    An <- Anova(object,type=3)
    effect.names <- rownames(An$anova)
    if(is.null(object$effect.sources))
      stop("Effect sources not found")
    errors <- An$errors
    err.df <- An$denom.df
    inds <- match(object$effect.sources[-1],effect.names)
    #  resvar <- c(rss/rdf,errors[inds]/err.df[inds])
    resvar <- c(rss/rdf,errors[inds])[Qr$pivot[p1]]
  }
  if (any(is.finite(resvar) & resvar < (mean(f)^2 + var(f)) * 
          1e-30)) 
    warning("essentially perfect fit: summary may be unreliable")
  R    <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se   <- sqrt(diag(R) * resvar)
  est  <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans  <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
                                                  rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$aliased <- is.na(coef(object))
  ans$sigma <- sqrt(resvar[1])
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
                                                       df.int)/rdf[1])
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar[1], 
                        numdf = p - df.int, dendf = rdf[1])
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
                                                             1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  if(!is.null(object$random) && !is.balanced(object)){
    cat("\nWARNING: Unbalanced data may lead to poor estimates\n")
  }
  class(ans) <- "summary.lmm"
  ans
}


###########################################
# confint.lm from stats, edited to enable limited classical least squares mixed models
confint.lmm <- function (object, parm, level = 0.95, ...) 
{
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qt(a, object$df.residual)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  if(is.null(object$random)){
    ses <- sqrt(diag(vcov(object)))[parm]
  } else {
    p <- object$rank
    p1 <- 1L:p # Moved up for the sake of "random"
    Qr <- qr.lmm(object)
    rdf <- object$df.residual
    r <- object$residuals
    w <- object$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    } else {
      rss <- sum(w * r^2)
    }
    An <- Anova(object,type=3)
    effect.names <- rownames(An$anova)
    if(is.null(object$effect.sources))
      stop("Effect sources not found")
    errors <- An$errors
    err.df <- An$denom.df
    inds <- match(object$effect.sources[-1],effect.names)
    resvar <- c(rss/rdf,errors[inds])[Qr$pivot[p1]]
    R    <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    ses   <- sqrt(diag(R) * resvar)[match(parm,colnames(vcov(object)))]
  }
  
  ci[] <- cf[parm] + ses %o% fac
  ci
}
format.perc <- function (x, digits, ...) 
  paste(format(100 * x, trim = TRUE, scientific = FALSE, digits = digits, ...), 
        "%")
