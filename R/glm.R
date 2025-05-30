##
## This file contains a MODIFIED COPY (2014) of the glm method from the base package stats (2014-10-31).
##

# The function replaces stats::glm, organizes fixed and random effects, removes r() from formula and parses to glmer.
# 
glm <- function(formula, family = gaussian, data, weights,
                subset, na.action, start = NULL,
                etastart, mustart, offset,
                control = list(...),
                model = TRUE, method = "glm.fit",
                x = FALSE, y = TRUE,
                contrasts = "contr.sum", REML = TRUE, ...)
{
  call <- match.call()
  ## family
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ## Edited by KHL
  mf <- match.call(expand.dots = FALSE)
  mfd <- match(c("formula","data"), names(mf), 0L)
  if(length(mfd)==2){ # Has formula and data
    is.random <- TRUE
    if( any(grepl("r(",formula,fixed=TRUE)) ){
      rw <- random.worker(formula, data, REML)
    } else {
      rw <- list(0)
    }
    if(length(rw) != 1){  # Removed r() from formula
      if(requireNamespace("lme4", quietly = TRUE)){
        formula <- rw$formula
        mf$formula <- rw$formula
        if(!is.logical(REML)){ # Perform 
          REML <- TRUE
          warning("REML must be logical")
        }
        if(family$family == "gaussian" && family$link =="identity"){
          object <- lme4::lmer(rw$reml.formula, data, REML = REML, contrasts = contrasts, na.action = na.action, ...)
        } else {
          object <- lme4::glmer(rw$reml.formula, data, family = family, contrasts = contrasts, na.action = na.action, ...)
        }
        object@call <- call
        return(object)
      } else {
        warning('Package lme4 required for random glm')
      }
    }
  }
  ## End of edit
  
  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame")) return(mf)
  
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  ## for back-compatibility in return result
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)
  
  mt <- attr(mf, "terms") # allow model.frame to have updated it
  
  Y <- model.response(mf, "any") # e.g. factors are allowed
  # Create contrast list if single character argument is supplied to contrasts (edit by KHL)
  contrasts.orig <- contrasts
  if(!is.null(contrasts)){
    if(is.character(contrasts) && length(contrasts)==1){
      facs <- which(unlist(lapply(mf, inherits, what = "factor")))
      if(contrasts == "contr.treatment.last"){ # Force last level of factor to base level
        contrasts <- lapply(mf[names(facs)], function(f){nl <- nlevels(f); contr.treatment(nl,nl)})
      } else
        if(contrasts == "contr.weighted")
          contrasts <- lapply(mf[facs], "contr.weighted")
        else {
          contrasts <- as.list(rep(contrasts, length(facs)))
          names(contrasts) <- names(facs)
        }
    }
  }
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  ## null model support
  if (!is.empty.model(mt)){
    X <- model.matrix(mt, mf, contrasts)
    # ## Edited by KHL 
    # if(is.null(contrasts) && (options("contrasts")[[1]][1]!="contr.treatment" || options("contrasts")[[1]][1]!="contr.poly") && !missing(data)){
    #   col.names   <- effect.labels(mt,data)
    #   if(length(col.names)==length(colnames(X))){
    #     colnames(X) <- effect.labels(mt,data)
    #   }
    # }
    # ## End edit
    ## Edited by KHL CCS = Cell Count Scaling
    col.names   <- effect.labels(mt,mf,contrasts) # mt is "terms" from formula, x is model.matrix
    if(length(col.names)==length(colnames(X))){
      colnames(X) <- col.names
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
              X[int_fac == lev, x_col] = x_sum[int_fac == lev, x_col] * min(n_each)/n_each[lev]
            }
          }
        }
      }
    }
    ## End edit
  } else {
    X <- matrix(,NROW(Y), 0L)
  }
  
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) )
    stop("negative weights not allowed")
  
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  
  ## We want to set the name on this call and the one below for the
  ## sake of messages from the fitter function
  fit <- eval(call(if(is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = family, control = control,
                   intercept = attr(mt, "intercept") > 0L))
  
  ## This calculated the null deviance from the intercept-only model
  ## if there is one, otherwise from the offset-only model.
  ## We need to recalculate by a proper fit if there is intercept and
  ## offset.
  ##
  ## The glm.fit calculation could be wrong if the link depends on the
  ## observations, so we allow the null deviance to be forced to be
  ## re-calculated by setting an offset (provided there is an intercept).
  ## Prior to 2.4.0 this was only done for non-zero offsets.
  if(length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <-
      eval(call(if(is.function(method)) "method" else method,
                x = X[, "(Intercept)", drop=FALSE], y = Y,
                weights = weights, offset = offset, family = family,
                control = control, intercept = TRUE))
    ## That fit might not have converged ....
    if(!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase maxit?")
    fit$null.deviance <- fit2$deviance
  }
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula,
                     terms = mt, data = data,
                     offset = offset, control = control, method = method,
                     contrasts = attr(X, "contrasts"),
                     xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, c("glm", "lm"))
  fit
}
