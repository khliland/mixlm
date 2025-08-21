##
## This file contains MODIFIED COPIES of the anova.lm method from the base package stats (2014-10-31),
## and Anova.lm from the package car (2014-10-10).
##

anova.lmm <- function(object, ...){
  if(!is.null(object$random)){
    if(inherits(object, "mlm"))
      return(AnovaMixMLM(object, 1))
    else
      return(AnovaMix(object, 1))
  } else {
    class(object) <- setdiff(class(object), "lmm")
    return(anova(object,...))
  }
}

Anova.lmm <- function(mod, ...){
  mf <- match.call()
  if(!is.null(mod$random)){
    if(mf$type=="I" || mf$type=="1"){
      if(inherits(mod, "mlm"))
        return(AnovaMixMLM(mod, 1))
      else
        return(AnovaMix(mod, 1))
    } else { 
      if(mf$type=="II" || mf$type=="2"){
        if(inherits(mod, "mlm"))
          return(AnovaMixMLM(mod, 2))
        else
          return(AnovaMix(mod, 2))
      } else {
        if(inherits(mod, "mlm"))
          return(AnovaMixMLM(mod, 3))
        else
          return(AnovaMix(mod, 3))
      }
    }
  } else {
    if(inherits(mod, "mlm"))
      class(mod) <- c("mlm", "lm")
    else
      class(mod) <- "lm"
    return(Anova(mod, ...))
  }
}

# Mixed model ANOVA
AnovaMix <- function(object, SStype){
  formula         <- formula(object)
  formula.text    <- as.character(formula)
  all.effects     <- object$random$all							  # All model effects and interactions
  fixed.effects   <- object$random$fixed							# All fixed effects
  random.effects  <- object$random$random						  # All random effects
  main.rands.only.inter <- object$random$main.rands.only.inter     # Random effects only present in interactions
  restrictedModel <- !object$random$unrestricted
  data    <- object$model
  n.effects    <- length(all.effects)
  main.effects <- fparse(formula)							  # All main effects (even though only included in interactions)
  n.levels     <- numeric(length(main.effects))
  for(i in 1:length(main.effects)){
    n.levels[i] <- length(levels(data[,main.effects[i]])) # Number of levels per main effect
  }
  names(n.levels) <- main.effects
  N <- dim(data)[1]
  
  ind.randoms <- numeric()
  ind.randoms <- match(random.effects,all.effects) # Placement of random effects in "all.effects"
  ind.fixed   <- match(fixed.effects,all.effects)  # Placement of fixed effects in "all.effects"
  ind.fixed   <- setdiff(1:n.effects,ind.randoms)										
  n.randoms   <- length(ind.randoms)
  
  # Estimate fixed effect Anova
  noRandom <- object
  noRandom$random <- NULL
  if(inherits(noRandom, "mlm"))
    class(noRandom) <- c("mlm", "lm")
  else
    class(noRandom) <- "lm"
  if(SStype == 1 || SStype == "I")
    fixed.model <- as.data.frame(stats::anova(noRandom))
  if(SStype == 2 || SStype == "II")
    fixed.model <- as.data.frame(car::Anova(noRandom, type='II', singular.ok=TRUE))
  if(SStype == 3 || SStype == "III")
    fixed.model <- as.data.frame(car::Anova(noRandom, type='III', singular.ok=TRUE))
  fixed.model <- fixed.model[c(all.effects,"Residuals"),] # Sort according to all.effects
  if(!any("Mean Sq"%in%colnames(fixed.model))){
    fixed.model <- cbind(fixed.model[,"Sum Sq"]/fixed.model[,"Df"], fixed.model)
    colnames(fixed.model)[1] <- "Mean Sq"
  }
  
  # Check which effects should use interactions as denominators instead of error
  approved.interactions <- list()
  approved.interactions.fixed <- list()
  for(i in 1:n.effects){
    this.effect <- strsplit(all.effects[i],":")[[1]]
    which.contains <- numeric()
    for(j in 1:n.effects){ # Find all other effects containing this.effect
      effect.names <- is.element(strsplit(all.effects[j],":")[[1]],this.effect)
      # Check if current effect is contained in another effect of higher interaction level
      if(i!=j && sum(effect.names)==length(this.effect) && length(effect.names)>length(this.effect)){
        which.contains <- union(which.contains,j)}
    }
    which.contains <- sort(which.contains)
    if(length(which.contains)>0){
      approved.interaction <- numeric(length(which.contains))
      approved.interaction.fixed <- numeric(length(which.contains))
      for(j in 1:length(which.contains)){
        if(restrictedModel){
          # Check if any of the other main effect contained in the higher order interaction is random
          approved.interaction[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),c(random.effects,main.rands.only.inter)))
        } else {
          if(any(is.element(ind.fixed,i))){
            # Check if any of the other main effects contained in the higher order interaction is fixed
            approved.interaction.fixed[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),fixed.effects))
          }
          # Check if all of the main effects contained in the higher order interaction are random
          approved.interaction[j] <- 1-prod(!is.element(strsplit(all.effects[which.contains],":")[[j]],c(random.effects,main.rands.only.inter)))
        }
      }
      if(length(which(approved.interaction==1))>0){
        approved.interactions[[i]] <- which.contains[which(approved.interaction==1)]}
      else{
        approved.interactions[[i]] <- FALSE}
      if(length(which(approved.interaction.fixed==1))>0){
        approved.interactions.fixed[[i]] <- which.contains[which(approved.interaction.fixed==1)]}
      else{
        approved.interactions.fixed[[i]] <- FALSE}
    }
    else{
      approved.interactions[[i]] <- FALSE
      approved.interactions.fixed[[i]] <- FALSE}
  }
  
  # Find variance components (except MSerror), 
  # and find linear combinations needed to produce denominators of F-statistics
  mix.model.attr <- list()
  denom.df <- numeric(n.effects+1)
  exp.mean.sq <- rep(paste("(",n.effects+1,")", sep=""), n.effects+1)
  var.comps <- numeric(n.effects+1)*NA
  var.comps[n.effects+1] <- fixed.model[n.effects+1,"Mean Sq"]
  errors <- numeric(n.effects)
  for(i in 1:n.effects) {
    if(!is.logical(approved.interactions[[i]])){
      # Set up matrix A and vector b to find linear combinations of effects to use as denominators in F statistics
      ## This is probably where unbalancedness should be included !!!!!!!!
      lap <- length(approved.interactions[[i]])
      A <- matrix(0,lap+1,n.effects+1)
      b <- rep(1,lap+1)
      for(j in 1:lap){
        A[j,approved.interactions[[approved.interactions[[i]][j]]]] <- 1
        A[j,approved.interactions[[i]][j]] <- 1
        k <- length(approved.interactions[[i]])+1-j
        exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[approved.interactions[[i]][k]],":")[[1]]]), " (",which(all.effects==all.effects[approved.interactions[[i]][k]]),")", sep="")
      }
      A[, n.effects+1] <- 1
      A <- A[,apply(A,2,sum)>0]
      denominator <- solve(t(A),b)
      denominator.id <- c(approved.interactions[[i]],n.effects+1)
      denominator.id <- denominator.id[denominator!=0]
      mix.model.attr[[i]] <- denominator <- denominator[denominator!=0]
      names(mix.model.attr[[i]]) <- denominator.id
      if(length(denominator)==1){ # Original df
        denom.df[i] <- fixed.model[denominator.id,"Df"]}
      else{ # Satterthwaite's df correction
        denom.df[i] <- sum(fixed.model[denominator.id,"Mean Sq"]*denominator)^2/sum((fixed.model[denominator.id,"Mean Sq"]*denominator)^2/fixed.model[denominator.id,"Df"])} 
    } else{
      denominator.id <- n.effects+1
      mix.model.attr[[i]] <- 1
      names(mix.model.attr[[i]]) <- denominator.id
      denom.df[i] <- fixed.model[denominator.id,"Df"]
      denominator <- 1
    }
    if(sum(ind.randoms==i)>0){
      exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " (",i,")", sep="")
      var.comps[i] <- (fixed.model[i,"Mean Sq"]-fixed.model[denominator.id,"Mean Sq"]%*%denominator)/(N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]))}
    else{
      if(!is.logical(approved.interactions.fixed[[i]])){
        ex.ind <- paste(",", paste(approved.interactions.fixed[[i]], sep="", collapse=","),sep="")}
      else{
        ex.ind <- ""}
      exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " Q[",i,ex.ind,"]", sep="")
    }
    errors[i] <- fixed.model[denominator.id,"Mean Sq"]%*%denominator
    fixed.model[i,"F value"] <- fixed.model[i,"Mean Sq"]/(fixed.model[denominator.id,"Mean Sq"]%*%denominator)
    if(is.na(fixed.model[i,"F value"]) || fixed.model[i,"F value"]<0){
      fixed.model[i,"F value"] <- NA
    }
    fixed.model[i,"Pr(>F)"] <- 1-pf(fixed.model[i,"F value"],fixed.model[i,"Df"],denom.df[i])
  }
  names(denom.df) <- rownames(fixed.model)
  object <- list(lm=object, anova=fixed.model, err.terms=c(mix.model.attr,NA), denom.df=denom.df, restricted=restrictedModel,
                 exp.mean.sq=exp.mean.sq, var.comps=var.comps, random.effects=random.effects, ind.randoms=ind.randoms, formula.text=formula.text, errors=errors)
  class(object) <- "AnovaMix"
  object
}


## Print method for object from AnovaMix
print.AnovaMix <- function(x,...){
  object <- x
  N <- length(object$err.terms)
  output1 <- object$anova
  Fs <- PrF <- character(N)
  PrF[!is.na(output1$"Pr(>F)")] <- format(round(output1$"Pr(>F)"[!is.na(output1$"Pr(>F)")],4), digits=1, scientific=FALSE, nsmall=4)
  PrF[is.na(output1$"Pr(>F)")] <- "-"
  output1$"Pr(>F)" <- PrF
  Fs[!is.na(output1$"F value")] <- format(output1$"F value"[!is.na(output1$"F value")], digits=1, scientific=FALSE, nsmall=2)
  Fs[is.na(output1$"F value")] <- "-"
  output1$"F value" <- Fs
  output1$"Sum Sq" <- format(output1$"Sum Sq", digits=1, scientific=FALSE, nsmall=2)
  output1$"Mean Sq" <- format(output1$"Mean Sq", digits=1, scientific=FALSE, nsmall=2)
  
  err.terms <- character(length(object$err.terms))
  for(i in 1:N){
    if(length(object$err.terms[[i]])==1 && is.na(object$err.terms[[i]])){
      err.terms[i] <- "-"
    }
    else{
      err.terms[i] <- paste(ifelse(object$err.terms[[i]][1]>1,paste(object$err.terms[[i]][1],"*",sep=""),""),"(",names(object$err.terms[[i]][1]),")",sep="")
      if(length(object$err.terms[[i]])>1){
        for(j in 2:length(object$err.terms[[i]])){
          if(object$err.terms[[i]][j]<0){
            err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]<(-1),paste(abs(object$err.terms[[i]][j]),"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" - ")
          } else {
            err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]>1,paste(object$err.terms[[i]][j],"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" + ")}
        }
      }
    }
  }
  var.comps <- format(object$var.comps, digits=3)
  var.comps[setdiff(1:(N-1), object$ind.randoms)] <- "fixed"
  
  denom.df <- character(N)
  denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df], digits=3)
  denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df], digits=3)
  denom.df[object$denom.df==0] <- "-"
  output2 <- data.frame("Err.terms"=err.terms, "Denom.df"=denom.df, "VC(SS)"=var.comps)
  colnames(output2) <- c("Err.term(s)", "Err.df", "VC(SS)")
  output3 <- data.frame("E(MS)"=format(object$exp.mean.sq))
  colnames(output3) <- "Expected mean squares"
  rownames(output2) <- paste(1:N," ",rownames(object$anova), sep="")
  rownames(output3) <- rownames(object$anova)
  if(!object$restricted){
    un <- "un"}
  else{
    un <- ""}
  cat("Analysis of variance (", un, "restricted model)\n", sep="")
  cat("Response: ", object$formula.text[2], "\n", sep="")
  print(format(output1, digits=3))
  cat("\n")
  print(output2)
  cat("(VC = variance component)\n\n")
  print(output3)
}

AnovaMixMLM <- function(object, SStype){
  formula         <- formula(object)
  formula.text    <- as.character(formula)
  all.effects     <- object$random$all							  # All model effects and interactions
  fixed.effects   <- object$random$fixed							# All fixed effects
  random.effects  <- object$random$random						  # All random effects
  main.rands.only.inter <- object$random$main.rands.only.inter     # Random effects only present in interactions
  restrictedModel <- !object$random$unrestricted
  data    <- object$model
  n.effects    <- length(all.effects)
  main.effects <- fparse(formula)							  # All main effects (even though only included in interactions)
  n.levels     <- numeric(length(main.effects))
  for(i in 1:length(main.effects)){
    n.levels[i] <- length(levels(data[,main.effects[i]])) # Number of levels per main effect
  }
  names(n.levels) <- main.effects
  N <- dim(data)[1]
  
  ind.randoms <- numeric()
  ind.randoms <- match(random.effects,all.effects) # Placement of random effects in "all.effects"
  ind.fixed   <- match(fixed.effects,all.effects)  # Placement of fixed effects in "all.effects"
  ind.fixed   <- setdiff(1:n.effects,ind.randoms)										
  n.randoms   <- length(ind.randoms)
  
  # Estimate fixed effect Anova
  noRandom <- object
  noRandom$random <- NULL
  if(inherits(noRandom, "mlm"))
    class(noRandom) <- c("mlm", "lm")
  else
    class(noRandom) <- "lm"
  univar <- noRandom
  class(univar) <- c("lm")
  fixed.models <- list()
  if(SStype == 1 || SStype == "I")
    for(i in 1:ncol(object$coefficients)){
      univar$coefficients <- object$coefficients[,i]
      univar$residuals <- object$residuals[,i]
      fixed.models[[i]] <- as.data.frame(stats::anova(univar))[c(all.effects,"Residuals"),]
    }
  if(SStype == 2 || SStype == "II"){
    # Loop over all responses, convert to univariate and perform ANOVA
    for(i in 1:ncol(object$coefficients)){
      univar$coefficients <- object$coefficients[,i]
      univar$residuals <- object$residuals[,i]
      fixed.models[[i]] <- as.data.frame(car::Anova(univar, type='II', singular.ok=TRUE))[c(all.effects,"Residuals"),]
    }
  }
  if(SStype == 3 || SStype == "III")
    for(i in 1:ncol(object$coefficients)){
      univar$coefficients <- object$coefficients[,i]
      univar$residuals <- object$residuals[,i]
      fixed.models[[i]] <- as.data.frame(car::Anova(univar, type='III', singular.ok=TRUE))[c(all.effects,"Residuals"),]
    }
  if(!any("Mean Sq"%in%colnames(fixed.models[[1]]))){
    for(i in 1:length(fixed.models)){
      fixed.models[[i]] <- cbind("Mean Sq"=fixed.models[[i]][,"Sum Sq"]/fixed.models[[i]][,"Df"], fixed.models[[i]])
      colnames(fixed.models[[i]])[1] <- "Mean Sq"
    }
  }
  
  # Check which effects should use interactions as denominators instead of error
  approved.interactions <- list()
  approved.interactions.fixed <- list()
  for(i in 1:n.effects){
    this.effect <- strsplit(all.effects[i],":")[[1]]
    which.contains <- numeric()
    for(j in 1:n.effects){ # Find all other effects containing this.effect
      effect.names <- is.element(strsplit(all.effects[j],":")[[1]],this.effect)
      # Check if current effect is contained in another effect of higher interaction level
      if(i!=j && sum(effect.names)==length(this.effect) && length(effect.names)>length(this.effect)){
        which.contains <- union(which.contains,j)}
    }
    which.contains <- sort(which.contains)
    if(length(which.contains)>0){
      approved.interaction <- numeric(length(which.contains))
      approved.interaction.fixed <- numeric(length(which.contains))
      for(j in 1:length(which.contains)){
        if(restrictedModel){
          # Check if any of the other main effect contained in the higher order interaction is random
          approved.interaction[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),c(random.effects,main.rands.only.inter)))
        } else {
          if(any(is.element(ind.fixed,i))){
            # Check if any of the other main effects contained in the higher order interaction is fixed
            approved.interaction.fixed[j] <- prod(is.element(setdiff(strsplit(all.effects[which.contains],":")[[j]],strsplit(all.effects[i],":")[[1]]),fixed.effects))
          }
          # Check if all of the main effects contained in the higher order interaction are random
          approved.interaction[j] <- 1-prod(!is.element(strsplit(all.effects[which.contains],":")[[j]],c(random.effects,main.rands.only.inter)))
        }
      }
      if(length(which(approved.interaction==1))>0){
        approved.interactions[[i]] <- which.contains[which(approved.interaction==1)]}
      else{
        approved.interactions[[i]] <- FALSE}
      if(length(which(approved.interaction.fixed==1))>0){
        approved.interactions.fixed[[i]] <- which.contains[which(approved.interaction.fixed==1)]}
      else{
        approved.interactions.fixed[[i]] <- FALSE}
    }
    else{
      approved.interactions[[i]] <- FALSE
      approved.interactions.fixed[[i]] <- FALSE}
  }
  
  # Find variance components (except MSerror), 
  # and find linear combinations needed to produce denominators of F-statistics
  mix.model.attr <- errorss <- var.compss <- list()
  for(j1 in 1:ncol(object$coefficients)){
    fixed.model <- fixed.models[[j1]]
    denom.df <- numeric(n.effects+1)
    exp.mean.sq <- rep(paste("(",n.effects+1,")", sep=""), n.effects+1)
    var.comps <- numeric(n.effects+1)*NA
    var.comps[n.effects+1] <- fixed.model[n.effects+1,"Mean Sq"]
    errors <- numeric(n.effects)
    for(i in 1:n.effects) {
      if(!is.logical(approved.interactions[[i]])){
        # Set up matrix A and vector b to find linear combinations of effects to use as denominators in F statistics
        ## This is probably where unbalancedness should be included !!!!!!!!
        lap <- length(approved.interactions[[i]])
        A <- matrix(0,lap+1,n.effects+1)
        b <- rep(1,lap+1)
        for(j in 1:lap){
          A[j,approved.interactions[[approved.interactions[[i]][j]]]] <- 1
          A[j,approved.interactions[[i]][j]] <- 1
          k <- length(approved.interactions[[i]])+1-j
          exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[approved.interactions[[i]][k]],":")[[1]]]), " (",which(all.effects==all.effects[approved.interactions[[i]][k]]),")", sep="")
        }
        A[, n.effects+1] <- 1
        A <- A[,apply(A,2,sum)>0]
        denominator <- solve(t(A),b)
        denominator.id <- c(approved.interactions[[i]],n.effects+1)
        denominator.id <- denominator.id[denominator!=0]
        mix.model.attr[[i]] <- denominator <- denominator[denominator!=0]
        names(mix.model.attr[[i]]) <- denominator.id
        if(length(denominator)==1){ # Original df
          denom.df[i] <- fixed.model[denominator.id,"Df"]}
        else{ # Satterthwaite's df correction
          denom.df[i] <- sum(fixed.model[denominator.id,"Mean Sq"]*denominator)^2/sum((fixed.model[denominator.id,"Mean Sq"]*denominator)^2/fixed.model[denominator.id,"Df"])} 
      } else{
        denominator.id <- n.effects+1
        mix.model.attr[[i]] <- 1
        names(mix.model.attr[[i]]) <- denominator.id
        denom.df[i] <- fixed.model[denominator.id,"Df"]
        denominator <- 1
      }
      if(sum(ind.randoms==i)>0){
        exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " (",i,")", sep="")
        var.comps[i] <- (fixed.model[i,"Mean Sq"]-fixed.model[denominator.id,"Mean Sq"]%*%denominator)/(N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]))
      } else{
        if(!is.logical(approved.interactions.fixed[[i]])){
          ex.ind <- paste(",", paste(approved.interactions.fixed[[i]], sep="", collapse=","),sep="")
        } else{
          ex.ind <- ""}
        exp.mean.sq[i] <- paste(exp.mean.sq[i], " + ", N/prod(n.levels[strsplit(all.effects[i],":")[[1]]]), " Q[",i,ex.ind,"]", sep="")
      }
      errors[i] <- fixed.model[denominator.id,"Mean Sq"]%*%denominator
      fixed.model[i,"F value"] <- fixed.model[i,"Mean Sq"]/(fixed.model[denominator.id,"Mean Sq"]%*%denominator)
      if(is.na(fixed.model[i,"F value"]) || fixed.model[i,"F value"]<0){
        fixed.model[i,"F value"] <- NA
      }
      fixed.model[i,"Pr(>F)"] <- 1-pf(fixed.model[i,"F value"],fixed.model[i,"Df"],denom.df[i])
    }
    fixed.models[[j1]] <- fixed.model
    errorss[[j1]] <- errors
    var.compss[[j1]] <- var.comps
  }
  if(!is.null(colnames(object$coefficients))){
    names(fixed.models) <- colnames(object$coefficients)
    names(errorss) <- colnames(object$coefficients)
    names(var.compss) <- colnames(object$coefficients)
  }
  names(denom.df) <- rownames(fixed.model)
  object <- list(lm=object, anova=fixed.models, err.terms=c(mix.model.attr,NA), denom.df=denom.df, restricted=restrictedModel,
                 exp.mean.sq=exp.mean.sq, var.comps=var.compss, random.effects=random.effects, ind.randoms=ind.randoms, formula.text=formula.text, errors=errorss)
  class(object) <- "AnovaMixMLM"
  object
}

## Print method for object from AnovaMixMLM
print.AnovaMixMLM <- function(x, var = 1, ...){
  object <- x
  if(is.numeric(var) && !is.null(names(object$anova)))
    var <- names(object$anova)[var]
  N <- length(object$err.terms)
  output1 <- object$anova[[var]]
  Fs <- PrF <- character(N)
  PrF[!is.na(output1$"Pr(>F)")] <- format(round(output1$"Pr(>F)"[!is.na(output1$"Pr(>F)")],4), digits=1, scientific=FALSE, nsmall=4)
  PrF[is.na(output1$"Pr(>F)")] <- "-"
  output1$"Pr(>F)" <- PrF
  Fs[!is.na(output1$"F value")] <- format(output1$"F value"[!is.na(output1$"F value")], digits=1, scientific=FALSE, nsmall=2)
  Fs[is.na(output1$"F value")] <- "-"
  output1$"F value" <- Fs
  output1$"Sum Sq" <- format(output1$"Sum Sq", digits=1, scientific=FALSE, nsmall=2)
  output1$"Mean Sq" <- format(output1$"Mean Sq", digits=1, scientific=FALSE, nsmall=2)
  
  err.terms <- character(length(object$err.terms))
  for(i in 1:N){
    if(length(object$err.terms[[i]])==1 && is.na(object$err.terms[[i]])){
      err.terms[i] <- "-"
    }
    else{
      err.terms[i] <- paste(ifelse(object$err.terms[[i]][1]>1,paste(object$err.terms[[i]][1],"*",sep=""),""),"(",names(object$err.terms[[i]][1]),")",sep="")
      if(length(object$err.terms[[i]])>1){
        for(j in 2:length(object$err.terms[[i]])){
          if(object$err.terms[[i]][j]<0){
            err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]<(-1),paste(abs(object$err.terms[[i]][j]),"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" - ")
          } else {
            err.terms[i] <- paste(err.terms[i], paste(ifelse(object$err.terms[[i]][j]>1,paste(object$err.terms[[i]][j],"*",sep=""),""), "(", names(object$err.terms[[i]][j]), ")", sep=""), sep=" + ")}
        }
      }
    }
  }
  var.comps <- format(object$var.comps[[var]], digits=3)
  var.comps[setdiff(1:(N-1), object$ind.randoms)] <- "fixed"
  
  denom.df <- character(N)
  denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)==object$denom.df], digits=3)
  denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df] <- format(object$denom.df[!is.na(object$denom.df)&round(object$denom.df)!=object$denom.df], digits=3)
  denom.df[object$denom.df==0] <- "-"
  output2 <- data.frame("Err.terms"=err.terms, "Denom.df"=denom.df, "VC(SS)"=var.comps)
  colnames(output2) <- c("Err.term(s)", "Err.df", "VC(SS)")
  output3 <- data.frame("E(MS)"=format(object$exp.mean.sq))
  colnames(output3) <- "Expected mean squares"
  rownames(output2) <- paste(1:N," ",rownames(object$anova[[var]]), sep="")
  rownames(output3) <- rownames(object$anova[[var]])
  if(!object$restricted){
    un <- "un"}
  else{
    un <- ""}
  cat("Analysis of variance (", un, "restricted model)\n", sep="")
  cat("Response: ", object$formula.text[2], ", variable: ", var, "\n", sep="")
  print(format(output1, digits=3))
  cat("\n")
  print(output2)
  cat("(VC = variance component)\n\n")
  print(output3)
}
