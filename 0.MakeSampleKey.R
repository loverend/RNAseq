#!/usr/bin/env Rscript

# Make a sample-fastq key for use downstream

sample.key <- read.delim("/gpfs2/well/immune-rep/shared/PUBLISHED_DATASETS/BULK_RNA/EBV_LCL", stringsAsFactors=FALSE)
# This is provided by Core with each lane of data received and looks like this:
# FileName	Sample Name
# WTCHG_384296_201191	S97
# WTCHG_384296_202179	S328
# WTCHG_384296_203167	S81
# Needs updating with new sample IDs

mapping.key <- matrix(nrow=length(unique(sample.key[, 2])), ncol=2)
mapping.key[, 1] <- unique(sample.key[, 4])
for(i in 1:nrow(mapping.key)){
  files <- c(sample.key[, 1][sample.key[, 4] == mapping.key[i, 1]])
  mapping.key[i, 2] <- paste(files, collapse=",")
}

mapping.key <- data.frame(mapping.key)

write.table(mapping.key, "sampleinfo/mapping.info.txt", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=FALSE)

sample.key$Run <- c(rep("A", 144), rep("B", 144), rep("C", 144), rep("D", 144))
sample.key$Run <- as.factor(sample.key$Run)

# Check coverage
library(ggplot2)
ggplot(sample.key, aes(Sample, Reads, group=Run, fill=Run)) +
  geom_col() +
  theme_bw() +
  geom_hline(yintercept = 50e6, lty=2)

ggplot(sample.key, aes(Sample, Reads, group=Run, colour=Run)) +
  geom_line() +
  theme_bw()
