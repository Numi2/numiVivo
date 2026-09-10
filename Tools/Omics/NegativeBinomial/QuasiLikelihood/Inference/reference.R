args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==3L)
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
stopifnot(as.character(packageVersion('edgeR'))=='4.10.5',as.character(packageVersion('limma'))=='3.68.5')
read <- function(path) {
  con <- gzfile(path,open='rt');on.exit(close(con))
  fromJSON(paste(readLines(con,warn=FALSE),collapse='\n'),simplifyVector=FALSE)
}
vec <- function(x) as.numeric(unlist(x)); mat <- function(x) do.call(rbind,lapply(x,vec))
i <- read(args[1]); native <- read(args[2]); result <- native$inference
stopifnot(native$completed,native$refitsExactlyEqual,native$residualDFExactlyEqual,
          native$abundanceExactlyEqual,native$moderationExactlyEqual,result$completed,
          identical(vec(i$featureIndices),vec(native$featureIndices)),length(result$failures)==0)
y <- mat(i$counts); x <- mat(i$design); offset <- matrix(vec(i$offsets),nrow(y),ncol(y),byrow=TRUE)
c <- vec(i$contrast); pivot <- which.max(abs(c)); free <- setdiff(seq_along(c),pivot)
stopifnot(length(free)>0,qr(x)$rank==ncol(x))
null_x <- x[,free,drop=FALSE]-outer(x[,pivot],c[free]/c[pivot])
phi <- vapply(i$expectedRefits,function(f)f$dispersion,numeric(1))
full_mu <- mat(lapply(i$expectedRefits,function(f)f$means))
null <- mglmLevenberg(y,design=null_x,dispersion=phi,offset=offset,maxit=500,tol=1e-12)
stopifnot(!any(null$failed))
reference_mu <- null$fitted.values
native_mu <- mat(lapply(result$tests,function(t)t$likelihoodRatio$nullMeans))
native_beta <- mat(lapply(result$tests,function(t)t$likelihoodRatio$nullCoefficients))
raw_lr <- vec(lapply(result$tests,function(t)t$likelihoodRatio$rawStatistic))
lr <- vec(lapply(result$tests,function(t)t$likelihoodRatio$statistic))
prior <- i$expectedModeration
post <- vec(prior$posteriorVariances); residual <- vec(i$expectedResidualDF)
denominator_df <- pmin(residual+vec(prior$priorDegreesOfFreedom),nrow(y)*(ncol(y)-ncol(x)))
reference_lr <- nbinomDeviance(y,reference_mu,phi)-nbinomDeviance(y,full_mu,phi)
same_native_lr <- nbinomDeviance(y,native_mu,phi)-nbinomDeviance(y,full_mu,phi)
relative <- function(a,b) max(abs(a-b)/pmax(1,abs(b)))
row_relative <- function(a,b) apply(abs(a-b)/pmax(1,abs(b)),1,max)
checks <- list()
check <- function(name,error,tolerance) {
 checks[[length(checks)+1L]] <<- list(name=name,error=error,tolerance=tolerance,passed=is.finite(error)&&error<=tolerance)
}
scaled_score <- apply(abs(((y-native_mu)/(1+phi*native_mu))%*%null_x)/sqrt((native_mu/(1+phi*native_mu))%*%(null_x^2)),1,max)
constraint <- as.vector(native_beta%*%c)
check('constrained-means',relative(native_mu,reference_mu),2e-6)
check('null-score',max(scaled_score),1e-7+1e-10)
check('contrast-zero',max(abs(constraint)/pmax(1,apply(abs(native_beta),1,max)*sum(abs(c)))),2e-12)
check('LR-tight-null',relative(raw_lr,reference_lr),2e-6)
check('LR-native-means',relative(raw_lr,same_native_lr),2e-6)
check('denominator-DF',relative(vec(lapply(result$tests,function(t)t$denominatorDegreesOfFreedom)),denominator_df),2e-12)
check('F-arithmetic',relative(vec(lapply(result$tests,function(t)t$fStatistic)),lr/post),2e-12)
log_p <- pf(lr/post,df1=1,df2=denominator_df,lower.tail=FALSE,log.p=TRUE)
p <- exp(log_p); q <- p.adjust(p,'BH')
check('log-F-tail',max(abs(vec(lapply(result$tests,function(t)t$logPValue))-log_p)),2e-7)
check('probability-arithmetic',relative(vec(lapply(result$tests,function(t)t$pValue)),p),2e-12)
check('BH',relative(vec(lapply(result$tests,function(t)t$adjustedPValue)),q),2e-12)
check('ordinary-DF-cap',abs(result$ordinaryResidualDFCap-nrow(y)*(ncol(y)-ncol(x))),0)
stopifnot(result$poissonBound=='not-applicable-to-modern-adjusted-QL')
out <- list(status=if(all(vapply(checks,function(c)c$passed,logical(1)))) 'passed' else 'failed',
 checks=checks,features=nrow(y),nativeCalls=sum(q<=.05),
 perGene=list(featureIndices=vec(i$featureIndices),nullMeanRelativeError=row_relative(native_mu,reference_mu),
   nullScaledScore=scaled_score,constraintResidual=constraint,referenceLR=reference_lr,sameNativeLR=same_native_lr,
   referenceF=lr/post,referenceDF=denominator_df,referenceLogP=log_p,referenceBH=q),
 versions=list(edgeR=as.character(packageVersion('edgeR')),limma=as.character(packageVersion('limma')),R=R.version.string))
write_json(out,args[3],auto_unbox=TRUE,digits=NA,null='null')
cat(out$status,'features',nrow(y),'calls',out$nativeCalls,'\n')
for(c in checks) if(!c$passed)cat('FAILED',c$name,c$error,'limit',c$tolerance,'\n')
quit(status=if(out$status=='passed') 0 else 1)
