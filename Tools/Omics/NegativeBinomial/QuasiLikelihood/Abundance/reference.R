args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
stopifnot(packageVersion('edgeR')=='4.10.5',packageVersion('limma')=='3.68.5',packageVersion('jsonlite')=='2.0.0')
input<-fromJSON(args[1],simplifyDataFrame=FALSE)
results<-lapply(input$rows,function(r) {
 warnings<-character()
 result<-tryCatch(withCallingHandlers({
  y<-matrix(r$counts,nrow=1)
  prior<-if(is.null(r$priorCount)) 2 else r$priorCount
  value<-aveLogCPM(y,offset=r$offsets,dispersion=r$dispersion,prior.count=prior)
  augmented<-addPriorCount(y,offset=r$offsets,prior.count=prior)
  tighter<-(mglmOneGroup(augmented$y,offset=augmented$offset,dispersion=r$dispersion,maxit=10000,tol=1e-12)+log(1e6))/log(2)
  list(value=unname(value),tighterValue=unname(tighter))
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e)))
 result$warnings<-warnings;result
})
write_json(list(results=results,session=capture.output(sessionInfo())),args[2],auto_unbox=TRUE,digits=17,na='string',null='null')
