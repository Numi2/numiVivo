# Stage diagnostics for the already scored fixed-normalization null benchmark.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[[1]],mustWork=TRUE)
output <- args[[2]]
stopifnot(!dir.exists(output)); dir.create(output,recursive=TRUE)
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
stopifnot(as.character(packageVersion("DESeq2"))=="1.52.0")
writeLines(capture.output(sessionInfo()),file.path(output,"session-info.txt"))
meta <- fromJSON(file.path(input,"input.json"))
y <- as.matrix(read.delim(file.path(input,"counts.tsv"),row.names=1,check.names=FALSE))
x <- as.matrix(read.delim(file.path(input,"design.tsv"),row.names=1,check.names=FALSE))
s <- read.delim(file.path(input,"samples.tsv"),check.names=FALSE)
stopifnot(identical(colnames(y),s$sampleID),identical(rownames(x),s$sampleID))
storage.mode(y) <- "integer"
keep <- rowSums(y)>=meta$minimumFeatureCounts & rowSums(y>0)>=meta$minimumExpressingPseudobulks
stopifnot(sum(keep)==meta$eligibleFeatures)
y <- y[keep,,drop=FALSE]
warnings <- messages <- character()
result <- tryCatch(withCallingHandlers({
  d <- DESeqDataSetFromMatrix(y,data.frame(row.names=s$sampleID,donor=s$donor,condition=s$condition),design=x)
  sizeFactors(d) <- s$sizeFactor
  d <- DESeq(d,test="Wald",fitType="parametric",betaPrior=FALSE,minReplicatesForReplace=Inf,quiet=TRUE,parallel=FALSE)
  r <- results(d,contrast=as.numeric(meta$contrast),independentFiltering=FALSE,cooksCutoff=FALSE)
  m <- mcols(d); f <- dispersionFunction(d)
  table <- data.frame(featureID=rownames(d),meanNormalizedCount=m$baseMean,
    geneWiseDispersion=m$dispGeneEst,trendDispersion=m$dispFit,
    finalDispersion=dispersions(d),dispersionOutlier=m$dispOutlier,
    log2FoldChange=r$log2FoldChange,standardError=r$lfcSE,pValue=r$pvalue,
    adjustedPValue=p.adjust(r$pvalue,"BH"),betaConverged=m$betaConv)
  stream <- gzfile(file.path(output,"stages.tsv.gz"),"wt")
  write.table(table,stream,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(stream)
  list(status="completed",residualDF=nrow(x)-ncol(x),priorLogVariance=attr(f,"dispPriorVar"),
    residualLogVariance=attr(f,"varLogDispEsts"),trendCoefficients=as.list(attr(f,"coefficients")),
    trendMethod=attr(f,"fitType"),eligibleGenes=nrow(y))
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},
message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")}),error=function(e)list(status="failed",error=conditionMessage(e)))
result$warnings <- warnings;result$messages <- messages
write_json(result,file.path(output,"stages.json"),pretty=TRUE,auto_unbox=TRUE,digits=16,null="null")
if(result$status!="completed")quit(status=1)
