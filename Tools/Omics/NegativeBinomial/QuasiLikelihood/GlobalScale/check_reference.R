args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3)
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
stopifnot(packageVersion('edgeR')=='4.10.5',packageVersion('limma')=='3.68.5')
i<-fromJSON(gzfile(args[1]));n<-fromJSON(gzfile(args[2]),simplifyDataFrame=FALSE);y<-as.matrix(i$counts);x<-as.matrix(i$design)
offset<-matrix(i$offsets,nrow=nrow(y),ncol=ncol(y),byrow=TRUE);phi<-n$trendDispersions
smootherChecks<-list()
for(k in seq_along(n$updates)) {
 u<-n$updates[[k]];e<-u$eligibleIndices+1L
 sm<-lowess(n$abundanceCovariates[e],u$quasiDispersions[e]^0.25,f=.5,iter=3,delta=.01*diff(range(n$abundanceCovariates[e])));f<-numeric(length(e));f[order(n$abundanceCovariates[e])]<-sm$y
 err<-abs(f-u$smoother)/pmax(1,abs(f));smootherChecks[[k]]<-list(maximumRelativeError=max(err),relativeErrors=err,reference=f,referenceScale=max(1,unname(quantile(sm$y,.9,type=7)))^4)
}
initial<-mglmLevenberg(y,design=x,offset=offset,dispersion=phi,tol=1e-12,maxit=10000)
final<-mglmLevenberg(y,design=x,offset=offset,dispersion=phi/n$averageQuasiDispersion,tol=1e-12,maxit=10000)
relative<-function(a,b) apply(abs(a-b)/pmax(1,abs(b)),1,max)
effect<-function(f) as.vector(f$coefficients%*%i$contrast)
score<-function(f,disp) apply(abs(((y-f$fitted.values)/(1+disp*f$fitted.values))%*%x)/sqrt((f$fitted.values/(1+disp*f$fitted.values))%*%(x^2)),1,max)
table<-data.frame(featureIndex=i$featureIndices,initialMeanError=relative(n$initialMeans,initial$fitted.values),finalMeanError=relative(n$finalMeans,final$fitted.values),initialEffectError=abs(n$initialEffects-effect(initial))/pmax(1,abs(effect(initial))),finalEffectError=abs(n$finalEffects-effect(final))/pmax(1,abs(effect(final))),initialReferenceScore=score(initial,phi),finalReferenceScore=score(final,phi/n$averageQuasiDispersion),initialReferenceFailed=initial$failed,finalReferenceFailed=final$failed)
errors<-unlist(table[c('initialMeanError','finalMeanError','initialEffectError','finalEffectError')]);smoothMax<-max(vapply(smootherChecks,function(s)s$maximumRelativeError,numeric(1)))
write_json(list(status=if(all(errors<=2e-5)&&smoothMax<=2e-7&&!any(initial$failed|final$failed)) 'passed' else 'failed',solverFunction='edgeR::mglmLevenberg',solverTolerance=1e-12,solverMaximumIterations=10000,maximumFitRelativeError=max(errors),maximumSmootherRelativeError=smoothMax,table=table,smootherChecks=smootherChecks,initialReferenceMeans=initial$fitted.values,finalReferenceMeans=final$fitted.values,session=capture.output(sessionInfo())),args[3],auto_unbox=TRUE,digits=17,na='string',null='null')
