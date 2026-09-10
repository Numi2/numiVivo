args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% c(2,3))
tight <- length(args)==3
if(tight) stopifnot(args[[3]]=="--tighter-full-fit")
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
stopifnot(as.character(packageVersion("edgeR"))=="4.10.5",as.character(packageVersion("jsonlite"))=="2.0.0")
input <- fromJSON(gzfile(args[[1]]),simplifyVector=FALSE)
results <- list(); messages <- character(); warnings <- character()
start <- proc.time()[[3]]
withCallingHandlers({
 for(group in input$groups) {
  x <- do.call(rbind,lapply(group$design,unlist));y <- do.call(rbind,lapply(group$counts,unlist))
  dispersion <- unlist(group$dispersions);offset <- matrix(unlist(group$offsets),nrow=nrow(y),ncol=ncol(y),byrow=TRUE)
  contrast <- unlist(group$contrast)
  fit <- glmFit(y,design=x,dispersion=dispersion,offset=offset,prior.count=0)
  refined <- tight && fit$method=="levenberg"
  if(refined) {
   # glmFit.default does not forward tol/maxit from ... to its likelihood fit.
   # Call that same exported owning solver explicitly; retain the model metadata.
   update <- mglmLevenberg(y,design=x,dispersion=dispersion,offset=offset,
    coef.start=fit$coefficients,tol=1e-10,maxit=1000)
   for(name in names(update)) fit[[name]] <- update[[name]]
  }
  test <- glmLRT(fit,contrast=contrast)
  # glmLRT does not return its constrained fit. Refit that same null space to
  # retain optimizer diagnostics instead of inferring convergence from a p-value.
  basis <- qr.Q(qr(matrix(contrast,ncol=1)),complete=TRUE)[,-1,drop=FALSE]
  x0 <- x %*% basis
  null <- glmFit(y,design=x0,dispersion=dispersion,offset=offset,prior.count=0)
  scaled <- function(f,design) {
   mu <- f$fitted.values
   apply(abs(((y-mu)/(1+dispersion*mu)) %*% design)/sqrt((mu/(1+dispersion*mu)) %*% (design^2)),1,max)
  }
  results[[length(results)+1]] <- data.frame(featureIndex=unlist(group$featureIndices),statistic=test$table$LR,pValue=test$table$PValue,
   fullSolver=fit$method,fullToleranceRefined=refined,fullFailureFlagAvailable=!is.null(fit$failed),nullFailureFlagAvailable=!is.null(null$failed),
   failed=(if(is.null(fit$failed)) FALSE else fit$failed) | (if(is.null(null$failed)) FALSE else null$failed),
   fullScaledScore=scaled(fit,x),nullScaledScore=scaled(null,x0),nullRefitLRDifference=null$deviance-fit$deviance-test$table$LR)
 }
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")})
write_json(list(results=do.call(rbind,results),warnings=warnings,messages=messages,seconds=proc.time()[[3]]-start,session=capture.output(sessionInfo()),tighterFullFit=tight,qualification="edgeR conditional full/null NB fits; fixed native dispersion and offsets; no extra count or coefficient prior"),args[[2]],auto_unbox=TRUE,pretty=FALSE,digits=17)
