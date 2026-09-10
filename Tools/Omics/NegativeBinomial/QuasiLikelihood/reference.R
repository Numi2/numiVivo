args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
expected<-c(edgeR="4.10.5",limma="3.68.5",statmod="1.5.2",jsonlite="2.0.0")
stopifnot(identical(vapply(names(expected),function(p)as.character(packageVersion(p)),character(1)),expected))
i<-fromJSON(gzfile(args[1]));y<-as.matrix(i$counts);x<-as.matrix(i$design);ng<-nrow(y)
stopifnot(ncol(y)==nrow(x),ng==length(i$featureIndices),qr(x)$rank==ncol(x),all(y>=0),all(y==floor(y)))
offset<-matrix(i$offsets,nrow=ng,ncol=ncol(y),byrow=TRUE)
results<-list();stages<-list();packageTrend<-NULL
for(owner in c("nativeTrend","edgeRTrend")) for(legacy in c(FALSE,TRUE)) {
 name<-paste(owner,if(legacy) "legacy" else "adjusted",sep="-");warnings<-messages<-character();start<-proc.time()[[3]]
 out<-tryCatch(withCallingHandlers({
  if(owner=="edgeRTrend" && is.null(packageTrend)) {
   d<-DGEList(y,lib.size=exp(i$offsets));d<-estimateDisp(d,x,robust=TRUE)
   packageTrend<-d$trended.dispersion
  }
  requested<-if(owner=="nativeTrend") i$nativeTrend else packageTrend
  fit<-glmQLFit(y,design=x,dispersion=requested,offset=offset,robust=TRUE,
    abundance.trend=TRUE,winsor.tail.p=c(.05,.1),legacy=legacy,prior.count=0,keep.unit.mat=!legacy)
  test<-glmQLFTest(fit,contrast=i$contrast);lr<-glmLRT(fit,contrast=i$contrast)
  poissonBound<-rep(0,ng)
  if(legacy) {
   need<-get('.isBelowPoissonBound',asNamespace('edgeR'))(fit)
   if(any(need)) {
    pf<-glmFit(y[need,,drop=FALSE],design=x,offset=offset[need,,drop=FALSE],dispersion=0,prior.count=0)
    poissonBound[need]<-glmLRT(pf,contrast=i$contrast)$table$PValue
   }
  }
  average<-if(legacy) 1 else fit$average.ql.dispersion
  fittedDispersion<-fit$dispersion/average
  score<-apply(abs(((y-fit$fitted.values)/(1+fittedDispersion*fit$fitted.values))%*%x)/sqrt((fit$fitted.values/(1+fittedDispersion*fit$fitted.values))%*%(x^2)),1,max)
  residual<-if(legacy) fit$df.residual.zeros else fit$df.residual.adj
  deviance<-if(legacy) fit$deviance else fit$deviance.adj
  table<-data.frame(featureIndex=i$featureIndices,requestedDispersion=requested,actualDispersion=fit$dispersion,
   fittedDispersion=fittedDispersion,averageQLScale=average,rawDeviance=fit$deviance,
   residualDeviance=deviance,residualDF=residual,priorDF=fit$df.prior,priorScale=fit$s2.prior,
   posteriorScale=fit$s2.post,LR=lr$table$LR,F=test$table$F,denominatorDF=test$df.total,
   pValue=test$table$PValue,BH=p.adjust(test$table$PValue,"BH"),poissonBound=poissonBound,
   scaledScore=score,failed=if(is.null(fit$failed)) FALSE else fit$failed)
  list(table=table,means=fit$fitted.values,unitDF=fit$unit.df.adj,unitDeviance=fit$unit.deviance.adj,
   leverage=fit$leverage,coefficients=fit$coefficients,ordinaryResidualDF=fit$df.residual,
   failureFlagAvailable=!is.null(fit$failed),solver=fit$method,
   methodConstraint=if(all(fit$dispersion==requested)) "passed" else "failed-dispersion-cap")
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart('muffleMessage')}),error=function(e)e)
 ok<-!inherits(out,'error');stages[[name]]<-list(status=if(ok) 'completed' else 'failed',warnings=warnings,messages=messages,seconds=proc.time()[[3]]-start,error=if(ok) NULL else conditionMessage(out))
 if(ok) results[[name]]<-out
 cat(name,stages[[name]]$status,'\n')
}
write_json(list(case=i$case,stages=stages,results=results,session=capture.output(sessionInfo())),args[2],auto_unbox=TRUE,digits=17,na='string',null='null')
if(length(results)!=4) quit(status=1)
