#!/usr/bin/env Rscript

# Make a list of the full paths to the final feature counts for all samples
basedir <- getwd()
myfiles <- list.files(paste0(basedir), recursive=TRUE, full.name=TRUE)
myfiles <- grep("counts.txt", myfiles, value=TRUE)
myfiles <- grep("summary", myfiles, value=TRUE, invert=TRUE)
sample.names <- basename(myfiles)
sample.names <- gsub(".counts.txt", "", sample.names)

# create a list to store single sample tables
DT <- list()

# read each file as array element of DT, keep gene id and count, and rename
# with sample name only
for (i in 1:length(myfiles) ) {
  tryCatch({
    DT[[myfiles[i]]] <- read.table(myfiles[i], header = T, stringsAsFactors = FALSE)
    colnames(DT[[myfiles[i]]]) <- c("ID", sample.names[i])
  }, error=function(e) print(myfiles[i]))
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# inspect
head(data)

# now add each additional table with the ID column as key
if(length(myfiles)>1){
	for (i in 2:length(myfiles)) {
	  y <- DT[[myfiles[i]]]
	  z <- merge(data, y, by = c("ID"))
	  data <- z
	}
}

# ID column becomes rownames
rownames(data) <- data$ID
genes <- rownames(data)
data <- data[, -1]

if(length(myfiles)==1){
	data <- data.frame(data)
	rownames(data) <- genes
	colnames(data) <- sample.names
}

# write data to file
write.table(data, "Full_count_data.txt", sep="\t", quote=FALSE)

##done