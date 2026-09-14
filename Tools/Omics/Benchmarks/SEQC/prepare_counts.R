args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
source <- args[1]; out <- args[2]
stopifnot(!dir.exists(out)); dir.create(out, recursive=TRUE)
for (site in c('AGR','BGI','CNL','COH','MAY','NVS')) {
  env <- new.env()
  file <- file.path(source, paste0('ILM_refseq_gene_',site,'.rda'))
  load(file, envir=env); x <- get(paste0('ILM_refseq_gene_',site), envir=env)
  stopifnot(identical(names(x)[1:4], c('EntrezID','Symbol','GeneLength','IsERCC')))
  columns <- grep('^[AB]_[0-9]+_', names(x), value=TRUE)
  stopifnot(length(columns)>0, !anyDuplicated(columns))
  library <- sub('^([AB]_[0-9]+)_.*$', '\\1', columns)
  ids <- sort(unique(library))
  values <- as.matrix(x[,columns,drop=FALSE])
  stopifnot(all(is.finite(values)), all(values>=0), all(values==floor(values)))
  counts <- sapply(ids, function(id) rowSums(values[,library==id,drop=FALSE]))
  stopifnot(all(counts<2^53), sum(counts)==sum(values))
  directory <- file.path(out,site); dir.create(directory)
  con <- gzfile(file.path(directory,'counts.tsv.gz'),'wt')
  write.table(counts,con,sep='\t',quote=FALSE,row.names=FALSE);close(con)
  write.table(x[,1:4],file.path(directory,'features.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  write.table(data.frame(sourceColumn=columns,library=library),file.path(directory,'lanes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(library=ids,condition=substr(ids,1,1),replicate=sub('^[AB]_','',ids),
                         laneColumns=as.integer(table(factor(library,levels=ids)))),
              file.path(directory,'libraries.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  cat(site,nrow(x),length(columns),length(ids),format(sum(counts),scientific=FALSE),'\n')
}
