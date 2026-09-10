args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
root <- normalizePath(args[[1]],mustWork=TRUE);out <- args[[2]]
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
stopifnot(as.character(packageVersion("DESeq2"))=="1.52.0")
records <- list()
for(study in c("kang","hagai"))for(seed in 1:10){
  path <- file.path(root,study,seed)
  n <- fromJSON(gzfile(file.path(path,"native-stages.json.gz")))
  r <- read.delim(file.path(path,"stages.tsv.gz"),check.names=FALSE)
  meta <- fromJSON(file.path(path,"stages.json"))
  nt <- n$features;at <- match(r$featureID,nt$featureID);stopifnot(!anyNA(at))
  nativeRefs <- n$trend$referenceFeatureIndices+1L
  rm <- match(nt$featureID[nativeRefs],r$featureID);stopifnot(!anyNA(rm))
  cases <- list(
    reference_original=list(means=r$meanNormalizedCount[r$geneWiseDispersion>=1e-6],alpha=r$geneWiseDispersion[r$geneWiseDispersion>=1e-6]),
    native_original=list(means=nt$meanNormalizedCount[nativeRefs],alpha=nt$geneWiseDispersion[nativeRefs]),
    reference_on_native_set=list(means=r$meanNormalizedCount[rm],alpha=r$geneWiseDispersion[rm]),
    native_above_reference_floor=list(means=nt$meanNormalizedCount[nativeRefs][nt$geneWiseDispersion[nativeRefs]>=1e-6],alpha=nt$geneWiseDispersion[nativeRefs][nt$geneWiseDispersion[nativeRefs]>=1e-6]))
  for(name in names(cases)){
    warnings <- messages <- character();c <- cases[[name]]
    value <- tryCatch(withCallingHandlers({
      f <- DESeq2:::parametricDispersionFit(c$means,c$alpha)
      baseline <- if(name=="reference_original")r$trendDispersion else n$trend$intercept+n$trend$inverseMeanCoefficient/r$meanNormalizedCount
      error <- max(abs(f(r$meanNormalizedCount)/baseline-1))
      if(name=="reference_original")stopifnot(error<1e-8)
      list(status="completed",coefficients=as.list(attr(f,"coefficients")),maximumCurveRelativeDifference=error)
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")} ),error=function(e)list(status="failed",error=conditionMessage(e)))
    value$study <- study;value$seed <- seed;value$case <- name;value$genes <- length(c$means)
    value$warnings <- warnings;value$messages <- messages
    records[[length(records)+1L]] <- value
    write_json(records,file.path(out,"trends.json"),pretty=TRUE,auto_unbox=TRUE,digits=16,null="null")
  }
}
write_json(list(status=if(all(vapply(records,function(r)r$status=="completed",logical(1))))"completed-all-eighty-curves" else "completed-with-failures",curves=length(records)),file.path(out,"complete.json"),pretty=TRUE,auto_unbox=TRUE)
