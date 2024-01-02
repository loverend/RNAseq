#!/usr/bin/env Rscript

# Mapping QC
library(ggplot2)

log.paths <- list.files(path=".", pattern="*Log.final.out", full.names=TRUE, recursive=TRUE)
sample.names <- list.files(path=".", pattern="*Log.final.out", full.names=FALSE, recursive=TRUE)
sample.names <- unlist(strsplit(sample.names, "/"))[seq(1, 2*length(sample.names), 2)]
logs <- lapply(log.paths, read.delim, row.names=1, header=FALSE)
logs <- lapply(logs, as.data.frame)

cbindlist <- function(list) {
  n <- length(list)
  res <- list[[1]]
  for (i in 2:n) res <- cbind(res, list[[i]])
  return(res)
}

for (i in 1:length(logs)){
	colnames(logs[[i]]) <- sample.names[i]
}

mapping.stats <- cbindlist(logs)
mapping.stats <- t(mapping.stats)
colnames(mapping.stats) <- gsub(" ", "", colnames(mapping.stats))
colnames(mapping.stats) <- gsub("\\|", "", colnames(mapping.stats))
mapping.stats <- data.frame(mapping.stats)

pdf("MappingQC.pdf", onefile=TRUE)

ggplot(mapping.stats, aes(Numberofinputreads, Uniquelymappedreadsnumber)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of uniquely mapping reads")

ggplot(mapping.stats, aes(Numberofinputreads, Numberofsplices.Total)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of splices identified")

mapping.stats$Uniquelymappedreads. <- as.numeric(as.character(gsub("%", "", 
                                            mapping.stats$Uniquelymappedreads.)))

ggplot(mapping.stats, aes(Numberofinputreads, Uniquelymappedreads.)) +
  geom_point() +
  theme_bw()

ggplot(mapping.stats, aes(Numberofinputreads, Mismatchrateperbase..)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Mismatch rate per base")

ggplot(mapping.stats, aes(Numberofinputreads, Numberofreadsmappedtomultipleloci)) +
  geom_point() +
  theme_bw() +
  xlab(label="Number of reads") +
  ylab(label="Number of multi-mapping reads")

dev.off()
