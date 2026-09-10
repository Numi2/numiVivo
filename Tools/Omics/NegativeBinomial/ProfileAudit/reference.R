args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
input <- normalizePath(args[[1]],mustWork=TRUE);output <- args[[2]]
stopifnot(!dir.exists(output));dir.create(output,recursive=TRUE)
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
stopifnot(as.character(packageVersion("DESeq2"))=="1.52.0")
y <- as.matrix(read.delim(file.path(input,"counts.tsv"),row.names=1,check.names=FALSE))
x <- as.matrix(read.delim(file.path(input,"design.tsv"),row.names=1,check.names=FALSE))
s <- read.delim(file.path(input,"samples.tsv"),check.names=FALSE);meta <- fromJSON(file.path(input,"input.json"))
stopifnot(identical(colnames(y),s$sampleID),identical(rownames(x),s$sampleID))
storage.mode(y) <- "integer";keep <- rowSums(y)>=meta$minimumFeatureCounts & rowSums(y>0)>=meta$minimumExpressingPseudobulks
stopifnot(sum(keep)==meta$eligibleFeatures);y <- y[keep,,drop=FALSE]
warnings <- messages <- character()
result <- tryCatch(withCallingHandlers({
 d <- DESeqDataSetFromMatrix(y,data.frame(row.names=s$sampleID,donor=s$donor,condition=s$condition),design=x)
 sizeFactors(d) <- s$sizeFactor;d <- DESeq2:::getBaseMeansAndVariances(d)
 initial <- pmin(pmax(1e-8,pmin(DESeq2:::roughDispEstimate(counts(d,normalized=TRUE),x),DESeq2:::momentsDispEstimate(d))),max(10,ncol(y)))
 d <- estimateDispersionsGeneEst(d,quiet=TRUE)
 mu <- assays(d)[["mu"]]
 grid <- DESeq2:::fitDispGridWrapper(y=y,x=x,mu=mu,logAlphaPriorMean=rep(0,nrow(y)),logAlphaPriorSigmaSq=1,usePrior=FALSE,
   weightsSEXP=matrix(1,nrow(y),ncol(y)),useWeightsSEXP=FALSE,weightThresholdSEXP=1e-2,useCRSEXP=TRUE)
 table <- data.frame(featureID=rownames(y),initialDispersion=initial,originalDispersion=mcols(d)$dispGeneEst,
   originalIterations=mcols(d)$dispGeneIter,gridDispersion=grid)
 f <- gzfile(file.path(output,"reference.tsv.gz"),"wt");write.table(table,f,sep="\t",quote=FALSE,row.names=FALSE);close(f)
 f <- gzfile(file.path(output,"means.tsv.gz"),"wt");write.table(data.frame(featureID=rownames(y),mu,check.names=FALSE),f,sep="\t",quote=FALSE,row.names=FALSE);close(f)
 list(status="completed",genes=nrow(y),residualDF=nrow(x)-ncol(x),gridUpper=max(10,ncol(y)))
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")}),error=function(e)list(status="failed",error=conditionMessage(e)))
result$warnings <- warnings;result$messages <- messages
write_json(result,file.path(output,"status.json"),pretty=TRUE,auto_unbox=TRUE,digits=16,null="null")
writeLines(capture.output(sessionInfo()),file.path(output,"session-info.txt"))
if(result$status!="completed")quit(status=1)
