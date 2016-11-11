library("pasilla")
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- "~/Downloads/pasilla/inst/extdata/pasilla_sample_annotation.csv"
countData <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
colData <- read.csv(pasAnno, row.names=1)
colData <- colData[,c("condition","type")]

rownames(colData) <- sub("fb","",rownames(colData))
all(rownames(colData) %in% colnames(countData))

countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
