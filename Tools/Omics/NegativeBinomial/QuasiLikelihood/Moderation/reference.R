#!/usr/bin/env Rscript
# Independent qualification against installed, pinned reference functions.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
library(jsonlite)
stopifnot(as.character(packageVersion('limma'))=='3.68.5',
          as.character(packageVersion('edgeR'))=='4.10.5')
read_gzip <- function(path) {
  con <- gzfile(path,open='rt'); on.exit(close(con))
  fromJSON(paste(readLines(con,warn=FALSE),collapse='\n'),simplifyVector=FALSE)
}
input <- read_gzip(args[1]); native <- read_gzip(args[2])
stopifnot(is.null(native$error), !is.null(native$fit))
v <- function(x) as.numeric(unlist(x))
x <- v(input$variances); df <- v(input$degreesOfFreedom); covariate <- v(input$abundance)
if (!length(covariate)) covariate <- NULL
robust <- if(is.null(input$robust)) TRUE else input$robust
n <- length(x); fit <- native$fit
stopifnot(identical(v(input$featureIndices),v(native$featureIndices)))
h <- get('logmdigamma',asNamespace('limma'))
tight_optimize <- function(f,interval,...) {
  out <- stats::optimize(f,interval,tol=1e-12)
  locations <- c(interval,out$minimum)
  values <- vapply(locations,f,numeric(1)); i <- which.min(values)
  list(minimum=locations[i],objective=values[i])
}
# Change only the optimizer binding in a runtime copy; installed source is
# neither modified nor copied into the native owner or this repository.
tight <- limma::fitFDistUnequalDF1
tight_environment <- new.env(parent=environment(tight))
optimizer_calls <- list()
tight_environment$optimize <- function(f,interval,...) {
  out <- tight_optimize(f,interval,...)
  optimizer_calls[[length(optimizer_calls)+1L]] <<- out
  out
}
environment(tight) <- tight_environment
default <- limma::fitFDistUnequalDF1(x,df,covariate=covariate,robust=robust)
reference <- tight(x,df,covariate=covariate,robust=robust)
posterior <- function(prior) {
  prior_df <- if(is.null(prior$df2.shrunk)) prior$df2 else prior$df2.shrunk
  (df*x+prior_df*prior$scale)/(df+prior_df)
}
relative <- function(a,b) max(abs(a-b)/pmax(1,abs(b)))
checks <- list()
check <- function(name,error,tolerance) {
  checks[[length(checks)+1L]] <<- list(name=name,error=error,tolerance=tolerance,
    passed=is.finite(error) && error<=tolerance)
}
eligible <- df>=.01 & x>0
floor <- 1e-12*median(x[eligible]); xpos <- pmax(x,floor)
effective_df <- pmax(df,0); effective_df[df<.01] <- 1
a <- effective_df/2; corrected <- log(xpos)+h(a)
span <- limma::chooseLowessSpan(n,small.n=500)
profile_checks <- list()
for (p in fit$profiles) {
  weights <- v(p$priorWeights); w <- weights/trigamma(a)
  trend <- v(p$logVarianceTrend)
  expected_trend <- if(is.null(covariate)) rep(sum(w*corrected)/sum(w),n) else
    limma::loessFit(corrected,covariate,weights=w/quantile(w,.75),min.weight=1e-8,
      max.weight=100,span=span,iterations=1)$fitted
  check(paste0('trend-',length(profile_checks)+1L),relative(trend,expected_trend),2e-7)
  objective <- function(parameter) {
    d <- parameter/(1-parameter); scaled <- d*exp(trend-h(d))
    -2*sum(weights*(lgamma(a+d)-lgamma(d)-a*log(scaled)-(a+d)*log1p(a*xpos/scaled)))
  }
  optimum <- tight_optimize(objective,c(.5,.9998))
  native_value <- objective(p$shape/(1+p$shape))
  check(paste0('objective-value-',length(profile_checks)+1L),abs(native_value-p$objective)/sum(weights),2e-8)
  check(paste0('objective-minimum-',length(profile_checks)+1L),max(0,native_value-optimum$objective)/sum(weights),2e-8)
  endpoints <- vapply(p$evaluations,function(e)e$parameter,numeric(1))
  stopifnot(.5 %in% endpoints,.9998 %in% endpoints)
  profile_checks[[length(profile_checks)+1L]] <- list(referenceOptimum=optimum,
    nativeObjective=native_value,nativeShape=p$shape,referenceTrend=expected_trend)
}
initial <- fit$profiles[[1]]
initial_scale <- exp(v(initial$logVarianceTrend)-h(initial$shape))
F <- x/initial_scale
right_log <- pf(F,effective_df,2*initial$shape,lower.tail=FALSE,log.p=TRUE)
left_log <- pf(F,effective_df,2*initial$shape,lower.tail=TRUE,log.p=TRUE)
if(robust) {
  right <- v(fit$screeningRightProbabilities)
  represented <- right>0 & is.finite(right_log)
  check('right-log-tail',max(c(0,abs(log(right[represented])-right_log[represented]))),2e-7)
  check('right-underflow',sum((right==0)!=(exp(right_log)==0)),0)
  weights <- p.adjust(2*exp(pmin(right_log,left_log)),'BH'); weights[weights>.3] <- 1
  check('screening-weights',relative(v(fit$screeningFDRWeights),weights),2e-7)
  if(length(fit$profiles)>1L) {
    check('refit-weights',relative(v(fit$profiles[[2]]$priorWeights),weights*as.numeric(df>=.01)),2e-7)
    expected_probability <- pmin(exp(right_log)/((n-rank(F)+.5)/n),1)
    check('outlier-probability',relative(v(fit$notOutlierProbabilities),expected_probability),2e-7)
  }
}
prior_df <- v(fit$priorDegreesOfFreedom); scales <- v(fit$priorScales)
check('posterior-arithmetic',relative(v(fit$posteriorVariances),(df*x+prior_df*scales)/(df+prior_df)),2e-12)
check('tight-reference-posterior',relative(v(fit$posteriorVariances),posterior(reference)),2e-5)
result <- list(status=if(all(vapply(checks,function(c)c$passed,logical(1)))) 'passed' else 'failed',
  versions=list(R=R.version.string,limma=as.character(packageVersion('limma')),edgeR=as.character(packageVersion('edgeR'))),
  checks=checks,profiles=profile_checks,tightOptimizerCalls=optimizer_calls,
  defaultReference=default,tightReference=reference,
  defaultPosteriorRelativeError=relative(v(fit$posteriorVariances),posterior(default)),
  tightCommonDFRelativeError=relative(fit$commonPriorDegreesOfFreedom,reference$df2),
  referencePosterior=posterior(reference))
write_json(result,args[3],auto_unbox=TRUE,digits=NA,null='null')
cat(result$status,'features',n,'profiles',length(fit$profiles),'tight posterior',tail(checks,1)[[1]]$error,'\n')
for(c in checks) if(!c$passed) cat('FAILED',c$name,c$error,'limit',c$tolerance,'\n')
quit(status=if(result$status=='passed') 0 else 1)
