#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
source <- args[1]; output <- args[2]
env <- new.env(parent = emptyenv())
load(source, envir = env)
stopifnot(identical(ls(env), "taqman"))
x <- get("taqman", envir = env)
expected <- c("EntrezID", "Symbol",
              as.vector(rbind(paste0("A", 1:4, "_value"), paste0("A", 1:4, "_detection"))),
              as.vector(rbind(paste0("B", 1:4, "_value"), paste0("B", 1:4, "_detection"))),
              as.vector(rbind(paste0("C", 1:4, "_value"), paste0("C", 1:4, "_detection"))),
              as.vector(rbind(paste0("D", 1:4, "_value"), paste0("D", 1:4, "_detection"))))
stopifnot(identical(names(x), expected), nrow(x) == 1044)
for (name in names(x)) if (is.factor(x[[name]])) x[[name]] <- as.character(x[[name]])
write.table(x, output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("rows", nrow(x), "columns", ncol(x), "\n")
