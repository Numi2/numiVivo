args <- commandArgs(TRUE); stopifnot(length(args)>=1, length(args)<=3)
root <- args[1]; out <- file.path(root,if(length(args)>=2) args[2] else 'reference');stopifnot(!dir.exists(out));dir.create(out)
suppressPackageStartupMessages({library(edgeR);library(limma);library(DESeq2);library(jsonlite)})
stopifnot(as.character(packageVersion('edgeR'))=='4.10.5',as.character(packageVersion('limma'))=='3.68.5',as.character(packageVersion('DESeq2'))=='1.52.0')
sites <- if(length(args)==3) strsplit(args[3],',',fixed=TRUE)[[1]] else c('AGR','BGI','CNL','COH','MAY','NVS')
stopifnot(all(sites %in% c('AGR','BGI','CNL','COH','MAY','NVS')), !anyDuplicated(sites))
write_gzip <- function(table, path) {
  con <- gzfile(path, 'wt'); on.exit(close(con))
  write.table(table, con, sep='\t', quote=FALSE, row.names=FALSE)
}
for(site in sites) {
  directory <- file.path(out,site);dir.create(directory)
  native <- fromJSON(gzfile(file.path(root,'native',paste0(site,'.json.gz'))))
  stopifnot(is.null(native$error),!is.null(native$result))
  d <- native$result$design
  y <- as.matrix(read.delim(gzfile(file.path(root,'counts',site,'counts.tsv.gz')),check.names=FALSE))
  y <- y[,d$sourcePseudobulkIndices+1,drop=FALSE]
  rownames(y) <- paste0('SEQC:RefSeq:row:',seq_len(nrow(y))-1)
  condition <- factor(substr(colnames(y),1,1),levels=c('A','B'))
  design <- model.matrix(~condition);sf <- d$sizeFactorValues
  stopifnot(all(design==d$rows),all(colSums(y)==d$libraryCounts),all(y>=0),all(y==floor(y)),max(y)<.Machine$integer.max)
  effective <- sf*exp(mean(log(colSums(y))))
  data <- DGEList(counts=y,lib.size=effective,norm.factors=rep(1,ncol(y)))
  data <- estimateDisp(data,design,robust=TRUE)
  fit <- glmFit(data,design,prior.count=0)
  edge <- glmLRT(fit,coef=2)$table
  write_gzip(data.frame(featureID=rownames(edge),edge),file.path(directory,'edgeR.tsv.gz'))
  v <- voom(DGEList(counts=y,lib.size=effective,norm.factors=rep(1,ncol(y))),design,plot=FALSE)
  lm <- eBayes(lmFit(v,design),robust=TRUE)
  tab <- topTable(lm,coef=2,number=Inf,sort.by='none')
  write_gzip(data.frame(featureID=rownames(tab),tab),file.path(directory,'limma.tsv.gz'))
  dds <- DESeqDataSetFromMatrix(round(y),data.frame(condition=condition,row.names=colnames(y)),~condition)
  sizeFactors(dds) <- sf
  dds <- DESeq(dds,fitType='parametric',betaPrior=FALSE,quiet=TRUE,minReplicatesForReplace=Inf)
  ds <- as.data.frame(results(dds,contrast=c('condition','B','A'),independentFiltering=FALSE,cooksCutoff=FALSE))
  write_gzip(data.frame(featureID=rownames(ds),ds),file.path(directory,'DESeq2.tsv.gz'))
  write_json(list(site=site,design=unname(design),sizeFactors=sf,features=nrow(y),libraries=ncol(y),
                  DESeq2FitType=attr(dispersionFunction(dds),'fitType'),normalization='matched-native-median-ratio',
                  edgeR=as.character(packageVersion('edgeR')),limma=as.character(packageVersion('limma')),DESeq2=as.character(packageVersion('DESeq2'))),
             file.path(directory,'receipt.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
  cat(site,'complete\n');flush.console()
}
