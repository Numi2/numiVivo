# Explicitly paired donor pseudobulk comparisons: native-normalized and package-normalized.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[[1]], mustWork=TRUE)
output <- args[[2]]
stopifnot(!dir.exists(output));dir.create(output, recursive=TRUE)
suppressPackageStartupMessages({library(edgeR);library(limma);library(DESeq2);library(jsonlite)})
expected_versions <- c(edgeR="4.10.5",limma="3.68.5",DESeq2="1.52.0",statmod="1.5.2",jsonlite="2.0.0")
stopifnot(identical(vapply(names(expected_versions),function(p)as.character(packageVersion(p)),character(1)),expected_versions))
writeLines(capture.output(sessionInfo()), file.path(output,"session-info.txt"))
meta <- fromJSON(file.path(input,"input.json"))
y <- as.matrix(read.delim(file.path(input,"counts.tsv"),row.names=1,check.names=FALSE))
x <- as.matrix(read.delim(file.path(input,"design.tsv"),row.names=1,check.names=FALSE))
samples <- read.delim(file.path(input,"samples.tsv"),check.names=FALSE)
stopifnot(identical(colnames(y), samples$sampleID), identical(rownames(x), samples$sampleID),qr(x)$rank==ncol(x),nrow(x)>ncol(x), all(is.finite(y)),all(y>=0),all(y==floor(y)),max(y)<=.Machine$integer.max)
storage.mode(y) <- "integer"
keep <- rowSums(y)>=meta$minimumFeatureCounts & rowSums(y>0)>=meta$minimumExpressingPseudobulks
stopifnot(sum(keep)==meta$eligibleFeatures)
y <- y[keep,,drop=FALSE]
contrast <- as.numeric(meta$contrast)
stopifnot(length(contrast)==ncol(x))
write.table(data.frame(featureID=rownames(y)),file.path(output,"tested-features.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
records <- list()
run_method <- function(name,fun) {
  warnings <- character();messages <- character();start <- proc.time()[[3]]
  result <- tryCatch(withCallingHandlers(fun(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")},message=function(m){messages<<-c(messages,conditionMessage(m));invokeRestart("muffleMessage")}),error=function(e)e)
  elapsed <- proc.time()[[3]]-start
  ok <- !inherits(result,"error")
  records[[name]] <<- list(status=if(ok) "completed" else "failed",seconds=elapsed,warnings=warnings,messages=messages,error=if(ok) NULL else conditionMessage(result))
  write_json(records,file.path(output,"runs.json"),pretty=TRUE,auto_unbox=TRUE,null="null",digits=16)
  if(ok) write.table(result,file.path(output,paste0(name,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
for(mode in c("native_size_factors","package_normalization")) {
  # A common multiplicative offset changes only the intercept. Keep effective
  # library sizes in count units for voom logCPM; native factors are relative.
  effective <- samples$sizeFactor * exp(mean(log(samples$libraryCounts)))
  make_y <- function() {
    d <- DGEList(counts=y,lib.size=if(mode=="native_size_factors") effective else samples$libraryCounts)
    if(mode=="package_normalization") d <- calcNormFactors(d,method="TMM")
    d
  }
  d <- make_y()
  write.table(data.frame(sampleID=samples$sampleID,librarySize=d$samples$lib.size,normFactor=d$samples$norm.factors,effectiveLibrary=d$samples$lib.size*d$samples$norm.factors),file.path(output,paste0(mode,"-edgeR-normalization.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  run_method(paste0(mode,"-edgeR-QL"),function() {
    fit <- glmQLFit(estimateDisp(d,x,robust=TRUE),x,robust=TRUE)
    test <- glmQLFTest(fit,contrast=contrast)
    t <- topTags(test,n=Inf,sort.by="none")$table
    data.frame(featureID=rownames(t),log2FoldChange=t$logFC,pValue=t$PValue,adjustedPValue=p.adjust(t$PValue,"BH"),statistic=t$F,averageLogCPM=t$logCPM,dispersion=fit$dispersion,qlPriorDF=fit$df.prior)
  })
  run_method(paste0(mode,"-limma-voom"),function() {
    v <- voom(d,x,plot=FALSE)
    fit <- eBayes(contrasts.fit(lmFit(v,x),matrix(contrast,ncol=1)),robust=TRUE)
    t <- topTable(fit,coef=1,number=Inf,sort.by="none")
    data.frame(featureID=rownames(t),log2FoldChange=t$logFC,pValue=t$P.Value,adjustedPValue=p.adjust(t$P.Value,"BH"),statistic=t$t,averageLogCPM=t$AveExpr)
  })
  run_method(paste0(mode,"-DESeq2"),function() {
    coldata <- data.frame(row.names=samples$sampleID,donor=samples$donor,condition=samples$condition)
    dds <- DESeqDataSetFromMatrix(y,coldata,design=x)
    if(mode=="native_size_factors") sizeFactors(dds) <- samples$sizeFactor
    # Common-universe comparison: no count replacement, independent filtering
    # or Cook's suppression. Also export the default result policy separately.
    dds <- DESeq(dds,test="Wald",fitType="parametric",betaPrior=FALSE,minReplicatesForReplace=Inf,quiet=TRUE,parallel=FALSE)
    rr <- results(dds,contrast=contrast,independentFiltering=FALSE,cooksCutoff=FALSE)
    default <- results(dds,contrast=contrast)
    write.table(data.frame(featureID=rownames(default),as.data.frame(default)),file.path(output,paste0(mode,"-DESeq2-default-results.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    write.table(data.frame(sampleID=samples$sampleID,sizeFactor=sizeFactors(dds)),file.path(output,paste0(mode,"-DESeq2-normalization.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    data.frame(featureID=rownames(rr),log2FoldChange=rr$log2FoldChange,pValue=rr$pvalue,adjustedPValue=p.adjust(rr$pvalue,"BH"),statistic=rr$stat,standardError=rr$lfcSE,dispersion=dispersions(dds),betaConverged=mcols(dds)$betaConv,dispersionOutlier=mcols(dds)$dispOutlier,maximumCooks=apply(assays(dds)[["cooks"]],1,max))
  })
}
write_json(list(status=if(all(vapply(records,function(r)r$status=="completed",logical(1)))) "completed" else "incomplete",packages=lapply(c("edgeR","limma","DESeq2","statmod","jsonlite"),function(p)list(package=p,version=as.character(packageVersion(p)))),input=meta,policy="Same predeclared count filter, observations and design; common-universe BH. QL F, voom moderated t and NB Wald tests remain distinct. DESeq2 default result filtering retained separately."),file.path(output,"summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
if(any(vapply(records,function(r)r$status!="completed",logical(1)))) quit(status=1)
