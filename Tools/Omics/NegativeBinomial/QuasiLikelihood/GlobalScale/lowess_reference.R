args<-commandArgs(trailingOnly=TRUE);suppressPackageStartupMessages(library(jsonlite));stopifnot(packageVersion('jsonlite')=='2.0.0')
i<-fromJSON(args[1],simplifyVector=FALSE)
results<-lapply(i,function(v){
 x<-unlist(v$x);y<-unlist(v$y);delta<-if(is.null(v$deltaFraction)) .01 else v$deltaFraction
 sm<-lowess(x,y,f=if(is.null(v$span)) .5 else v$span,iter=if(is.null(v$robustnessIterations)) 3L else v$robustnessIterations,delta=delta*diff(range(x)))
 out<-numeric(length(x));out[order(x)]<-sm$y
 list(id=v$id,fitted=out)
})
write_json(list(results=results,session=capture.output(sessionInfo())),args[2],auto_unbox=TRUE,digits=17)
