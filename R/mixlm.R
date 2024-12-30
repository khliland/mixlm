# Startup
# Set contrasts and report
.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
#  options(contrasts=c('contr.sum','contr.poly'))
  if (interactive()) {
    packageStartupMessage('mixlm ',as.character(utils::packageVersion("mixlm")))
  }
}

# Effect labels
effect.labels <- function(t, data, contrasts){
	effects   <- attr(t, "term.labels")
	factors   <- attr(t, "factors")
	intercept <- attr(t, "intercept")
	n.eff     <- length(effects)
	if(n.eff==0){
		return(NULL)
	}
	split.effects <- strsplit(effects,":")
	names     <- "(Intercept)"

	for(i in 1:n.eff){ 
		cur <- split.effects[[i]]
		if(i == 1 && intercept == 0){
			levs <- levels(data[[cur]])
			names <- paste(cur,"(",levs,")",sep="")
		} else {
			inter <- list()
			for(j in 1:length(cur)){
			  if(inherits(data[[cur[j]]],'factor')){ # Handle factor main effect
					levs <- levels(data[[cur[j]]])
					# Individual factor contrast handling
					csum <- ifelse(contrasts[[cur[j]]][1] %in% c("contr.sum","contr.treatment.last"),	TRUE, FALSE)
					if(csum){
						n.lev <- length(levs)
						# Nested factor contrast handling
						if(factors[cur[j],i]==2){
							inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
						} else {
							inter[[j]] <- paste(cur[j],"(",levs[-n.lev],")",sep="")
						}
					} else {
					  if(contrasts[[cur[j]]][1]=="contr.weighted"){
					    omit <- which(!(levs %in% colnames(contr.weighted(data[[cur[j]]]))))
					  } else {
					    omit <- 1
					  }
						if(factors[cur[j],i]==2){
							inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
						} else {
							inter[[j]] <- paste(cur[j],"(",levs[-omit],")",sep="")
						}
					}			
				} else {
					inter[[j]] <- cur[j]
				}
			}
			names <- c(names, apply(expand.grid(inter),1,paste,sep="",collapse=":"))
		}
	}
	names
}


effect.source <- function(t,data){
  # csum is a silly check, but most of this code is solely for the purpose of finding the number of effect levels.
  csum <- ifelse(options("contrasts")[[1]][1] %in% c("contr.weighted","contr.sum"),	TRUE, FALSE)
  effects   <- attr(t, "term.labels")
  factors   <- attr(t, "factors")
  intercept <- attr(t, "intercept")
  n.eff     <- length(effects)
  if(n.eff==0){
    return(NULL)
  }
  split.effects <- strsplit(effects,":")
  names     <- "(Intercept)"
  
  for(i in 1:n.eff){ 
    cur <- split.effects[[i]]
    if(i == 1 && intercept == 0){
      levs <- levels(data[[cur]])
      names <- paste(cur,"(",levs,")",sep="")
    } else {
      inter <- list()
      for(j in 1:length(cur)){
        if(inherits(data[[cur[j]]], "factor")){ # Handle factor main effect
          levs <- levels(data[[cur[j]]])
          if(csum){
            n.lev <- length(levs)
            if(factors[cur[j],i]==2){
              inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
            } else {
              inter[[j]] <- paste(cur[j],"(",levs[-n.lev],")",sep="")
            }
          } else {
            if(factors[cur[j],i]==2){
              inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
            } else {
              inter[[j]] <- paste(cur[j],"(",levs[-1],")",sep="")
            }
          }			
        } else {
          inter[[j]] <- cur[j]
        }
      }
      tmp <- apply(expand.grid(inter),1,paste,sep="",collapse=":")
      names <- c(names, rep(paste(cur,sep="",collapse=":"),length(tmp)))
    }
  }
  names
}

