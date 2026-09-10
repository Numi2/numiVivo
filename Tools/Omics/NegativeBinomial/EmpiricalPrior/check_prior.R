# Independent installed-DESeq2 weighted-quantile reference on exact native inputs.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
root <- args[[1]];protocol <- fromJSON(file.path(root,"protocol.json"),simplifyVector=FALSE)
records <- lapply(protocol$records,function(rec) {
  d<-file.path(root,rec$id);x<-read.delim(file.path(d,"prior-input.tsv"),check.names=FALSE)
  keep<-abs(x$effectLog2)<10;stopifnot(sum(keep)>=20)
  weights<-1/(1/x$mean[keep]+x$trendDispersion[keep])
  raw<-sqrt(DESeq2:::matchWeightedUpperQuantileForVariance(x$effectLog2[keep],weights,upperQuantile=.05))
  list(case=rec$id,referenceFeatureIndices=unname(x$featureIndex[keep]),rawSDLog2=raw,
    priorSDLog2=max(.01,raw),weightSum=sum(weights),effectiveReferenceFeatures=sum(weights)^2/sum(weights^2))
})
write_json(list(method="Installed DESeq2 matchWeightedUpperQuantileForVariance; native full-design tested effects, means and trend dispersions; explicit abs(effect)<10 and SD floor .01",DESeq2=as.character(packageVersion("DESeq2")),records=records),file.path(root,"prior-reference.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(root,"prior-reference-session.txt"))
