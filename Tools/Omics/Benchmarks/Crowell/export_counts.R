# Author SCE -> count-assay wire; original Rda remains the complete source object.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
suppressPackageStartupMessages({library(SingleCellExperiment);library(jsonlite)})
load(args[[1]]);out<-args[[2]];dir.create(out,recursive=TRUE,showWarnings=FALSE)
y<-assay(sce,"counts");stopifnot(inherits(y,"dgCMatrix"),all(is.finite(y@x)),all(y@x>0),all(y@x==floor(y@x)),max(y@x)<=.Machine$integer.max)
write_bin<-function(name,value) { c<-gzfile(file.path(out,name),"wb");on.exit(close(c));writeBin(as.integer(value),c,size=4,endian="little") }
write_bin("values.i32.gz",y@x);write_bin("indices.i32.gz",y@i);write_bin("indptr.i32.gz",y@p)
column<-function(x) list(type=if(is.factor(x)) "factor" else typeof(x),values=if(is.factor(x)) as.character(x) else unname(x),levels=if(is.factor(x)) levels(x) else NULL,ordered=is.ordered(x))
meta<-list(dim=dim(y),nonzeros=length(y@x),umiTotal=sum(y@x),obsNames=colnames(sce),varNames=rownames(sce),obs=lapply(as.data.frame(colData(sce)),column),var=lapply(as.data.frame(rowData(sce)),column),assaySelection="counts only; original Rda retains logcounts and source reductions",sourceAssays=assayNames(sce),sourceReducedDimensions=reducedDimNames(sce),experimentInfo=metadata(sce)$experiment_info)
c<-gzfile(file.path(out,"metadata.json.gz"),"wt");writeLines(toJSON(meta,auto_unbox=TRUE,na="null",null="null",digits=17),c);close(c)
write_json(list(rows=nrow(y),columns=ncol(y),nonzeros=length(y@x),umiTotal=sum(y@x),packages=list(SingleCellExperiment=as.character(packageVersion("SingleCellExperiment")),Matrix=as.character(packageVersion("Matrix"))),scope=meta$assaySelection),file.path(out,"export.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
writeLines(capture.output(sessionInfo()),file.path(out,"session-info.txt"))
