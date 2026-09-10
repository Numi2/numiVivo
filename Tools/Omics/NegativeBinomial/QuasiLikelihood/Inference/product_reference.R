args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3L)
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
stopifnot(as.character(packageVersion('edgeR'))=='4.10.5',as.character(packageVersion('limma'))=='3.68.5')
read <- function(p) {con<-gzfile(p,open='rt');on.exit(close(con));fromJSON(paste(readLines(con,warn=FALSE),collapse='\n'),simplifyVector=FALSE)}
v <- function(x)as.numeric(unlist(x));m<-function(x)do.call(rbind,lapply(x,v))
out<-read(args[1]);frozen<-read(args[2]);r<-out$expression;q<-r$negativeBinomial$quasiLikelihood
indices<-v(q$featureIndices)+1L
fits<-lapply(indices,function(g)r$negativeBinomial$features[[g]]$finalFit)
x<-m(r$design$rows);y<-m(frozen$counts);mu<-m(lapply(fits,function(f)f$means))
phi<-v(lapply(fits,function(f)f$dispersion));c<-v(r$design$contrast)
offset<-matrix(log(v(r$design$sizeFactorValues)),nrow(y),ncol(y),byrow=TRUE)
pivot<-which.max(abs(c));free<-setdiff(seq_along(c),pivot)
nx<-x[,free,drop=FALSE]-outer(x[,pivot],c[free]/c[pivot])
null<-mglmLevenberg(y,nx,dispersion=phi,offset=offset,maxit=500,tol=1e-12);stopifnot(!any(null$failed))
tests<-q$inference$tests;actual_null<-m(lapply(tests,function(t)t$likelihoodRatio$nullMeans))
lr<-v(lapply(tests,function(t)t$likelihoodRatio$statistic))
reference_lr<-nbinomDeviance(y,null$fitted.values,phi)-nbinomDeviance(y,mu,phi)
post<-v(q$moderation$posteriorVariances);prior_df<-v(q$moderation$priorDegreesOfFreedom)
# Original native adjusted DF belongs to the same count/design/dispersion/mean
# problem. Validate those inputs below; this comparison is numerical, not exact
# equality after the public owner uses a different common offset constant.
df<-pmin(v(frozen$expectedResidualDF)+prior_df,nrow(y)*(ncol(y)-ncol(x)))
logp<-pf(lr/post,1,df,lower.tail=FALSE,log.p=TRUE)
rel<-function(a,b)max(abs(a-b)/pmax(1,abs(b)))
checks<-list();check<-function(name,error,tolerance){checks[[length(checks)+1L]]<<-list(name=name,error=error,tolerance=tolerance,passed=is.finite(error)&&error<=tolerance)}
old_mu<-m(lapply(frozen$expectedRefits,function(f)f$means));old_phi<-v(lapply(frozen$expectedRefits,function(f)f$dispersion))
check('published-refit-means',rel(mu,old_mu),2e-6);check('published-refit-dispersion',rel(phi,old_phi),2e-12)
score<-function(mm,xx)max(abs(((y-mm)/(1+phi*mm))%*%xx)/sqrt((mm/(1+phi*mm))%*%(xx^2)))
check('full-score',score(mu,x),1e-7+1e-10);check('null-score',score(actual_null,nx),1e-7+1e-10)
check('constrained-means',rel(actual_null,null$fitted.values),2e-6)
check('LR',rel(lr,reference_lr),2e-6)
check('denominator-DF',rel(v(lapply(tests,function(t)t$denominatorDegreesOfFreedom)),df),2e-12)
check('log-F-tail',max(abs(v(lapply(tests,function(t)t$logPValue))-logp)),2e-7)
check('BH',rel(v(lapply(tests,function(t)t$adjustedPValue)),p.adjust(exp(logp),'BH')),2e-12)
result<-list(status=if(all(vapply(checks,function(c)c$passed,logical(1))))'passed' else 'failed',checks=checks,
 qualification='Existing pseudobulk owner on unchanged public Kang aggregation; independent constrained fits and F/BH arithmetic, with offset-invariance comparison to published adjusted-DF inputs; not a new raw-H5AD import or calibration proof.')
write_json(result,args[3],auto_unbox=TRUE,digits=NA)
cat(result$status,'\n');for(c in checks)if(!c$passed)cat('FAILED',c$name,c$error,'\n')
quit(status=if(result$status=='passed')0 else 1)
