args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==4)
input <- normalizePath(args[[1]],mustWork=TRUE);audit <- normalizePath(args[[2]],mustWork=TRUE)
original <- normalizePath(args[[3]],mustWork=TRUE);output <- args[[4]]
stopifnot(!dir.exists(output));dir.create(output,recursive=TRUE)
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
stopifnot(as.character(packageVersion("DESeq2"))=="1.52.0")
y <- as.matrix(read.delim(file.path(input,"counts.tsv"),row.names=1,check.names=FALSE));storage.mode(y) <- "integer"
x <- as.matrix(read.delim(file.path(input,"design.tsv"),row.names=1,check.names=FALSE))
s <- read.delim(file.path(input,"samples.tsv"),check.names=FALSE);meta <- fromJSON(file.path(input,"input.json"))
stopifnot(identical(colnames(y),s$sampleID),identical(rownames(x),s$sampleID))
keep <- rowSums(y)>=meta$minimumFeatureCounts & rowSums(y>0)>=meta$minimumExpressingPseudobulks;y <- y[keep,,drop=FALSE]
checked <- read.delim(gzfile(file.path(audit,"checked.tsv.gz")),check.names=FALSE)
stopifnot(identical(rownames(y),checked$featureID))
records <- list()
finish <- function(d){
 d <- estimateDispersionsFit(d,fitType="parametric",quiet=TRUE)
 d <- estimateDispersionsMAP(d,quiet=TRUE)
 d <- nbinomWaldTest(d,betaPrior=FALSE,modelMatrix=x,quiet=TRUE)
 r <- results(d,contrast=as.numeric(meta$contrast),independentFiltering=FALSE,cooksCutoff=FALSE)
 list(data=data.frame(featureID=rownames(d),log2FoldChange=r$log2FoldChange,standardError=r$lfcSE,
  pValue=r$pvalue,adjustedPValue=p.adjust(r$pvalue,"BH"),dispersion=dispersions(d),dispersionOutlier=mcols(d)$dispOutlier),
  nonconvergedCoefficients=sum(!mcols(d)$betaConv),
  priorLogVariance=attr(dispersionFunction(d),"dispPriorVar"),residualLogVariance=attr(dispersionFunction(d),"varLogDispEsts"),
  trendCoefficients=as.list(attr(dispersionFunction(d),"coefficients")),trendMethod=attr(dispersionFunction(d),"fitType"))
}
warnings <- messages <- character()
result <- tryCatch(withCallingHandlers({
 d <- DESeqDataSetFromMatrix(y,data.frame(row.names=s$sampleID,donor=s$donor,condition=s$condition),design=x)
 sizeFactors(d) <- s$sizeFactor;d <- estimateDispersionsGeneEst(d,quiet=TRUE)
 stopifnot(isTRUE(all.equal(mcols(d)$dispGeneEst,checked$originalDispersion,tolerance=1e-12)))
 baseline <- finish(d)
 old <- read.delim(gzfile(original),check.names=FALSE);stopifnot(identical(old$featureID,baseline$data$featureID))
 columns <- c("log2FoldChange","standardError","pValue","adjustedPValue","dispersion")
 errors <- vapply(columns,function(col){stopifnot(isTRUE(all.equal(old[[col]],baseline$data[[col]],tolerance=1e-10)));max(abs(old[[col]]-baseline$data[[col]]),na.rm=TRUE)},numeric(1))
 accepted <- checked$gridObjectiveGain>1e-4 & checked$gridDispersion>=1e-8 & checked$gridDispersion<=max(10,ncol(y))
 revised <- d;mcols(revised)$dispGeneEst[accepted] <- checked$gridDispersion[accepted]
 candidate <- finish(revised)
 for(name in c("baseline","candidate")){
  value <- get(name);stream <- gzfile(file.path(output,paste0(name,".tsv.gz")),"wt")
  write.table(value$data,stream,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(stream)
  value$testedGenes <- sum(is.finite(value$data$pValue));value$bhBelow005 <- sum(value$data$adjustedPValue<.05,na.rm=TRUE)
  value$data <- NULL;records[[name]] <- value
 }
 records$changedGenes <- sum(accepted);records$outsideBoundsGridGenes <- sum(checked$gridDispersion<1e-8 | checked$gridDispersion>max(10,ncol(y)))
 records$baselineMaximumErrors <- as.list(errors)
 list(status="completed",results=records)
},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")}),error=function(e)list(status="failed",error=conditionMessage(e)))
result$warnings <- warnings;result$messages <- messages
write_json(result,file.path(output,"status.json"),pretty=TRUE,auto_unbox=TRUE,digits=16,null="null")
if(result$status!="completed")quit(status=1)
