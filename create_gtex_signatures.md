---
title: Deriving tissue signatures from the GTEx dataset. 
output:
  md_document:
    variant: markdown_phpextra
    preserve_yaml: TRUE
---

The *GTEx-project*[^1] is a data resource providing high quality gene
expression data from multiple reference human tissues. In this document,
we identify differentially expressed genes for each tissue using
[*voom*](http://bioconductor.org/packages/release/bioc/html/limma.html)
and
[*limma*](http://bioconductor.org/packages/release/bioc/html/limma.html)
to build tissue-specific signatures.

These signatures can be used with
[*BioQC*](https://accio.github.io/BioQC/) to detect tissue heterogeneity
in gene expression data.

### Load the data {#load-the-data}

We read the raw counts from the `.gct` file and preprocessed tissue
annotations from `GTEx-UDISDataSetID5681-sampleAnnotation.txt`.

~~~~ r
# loading big files takes some time, not necessary if we have the cache
if(!file.exists("data/limma_model.RData")) {
  gct = read_gct_matrix("data/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct")
  annotation = read_tsv("data/gtex/GTEx-UDISDataSetID5681-sampleAnnotation.txt")
  annotation.tissue = read_tsv("data/gtex/v6_UDIS.txt")
  ids = colnames(gct)
  annotation.matched = matchColumn(ids, annotation, "Experiment name")
  annotation.tissue.matched = matchColumn(ids, annotation.tissue, "SAMPID")
  
  colData = cbind(annotation.matched, annotation.tissue.matched)
  rownames(colData) = colData[["Experiment name"]]
  
  #filter rows with NA
  # RIN = RNA quality index
  # SMTS = tissue
  validInds = which(!apply(is.na(colData[,c("RIN", "Gender", "SMTS")]), 1, sum) >= 1)
  gct = gct[,validInds]
  colData = colData[validInds,]
}
~~~~

### Perform differential analysis {#perform-differential-analysis}

We first filter the GTEx expression data to only include genes for which
at least 10 samples have a count-per-million (CPM) greater than 1.

We then perform a differntial expression analysis for tissue on the
filtered dataset with `Gender` and `RIN` (RNA quality index) as
covariates.

We call a differetially expressed if it meets both of the following
criteria:

-   fold-change of &gt;100 (log2: 6.64)
-   p-value &lt; 0.01 (Bonferroni corrected)

~~~~ r
# read file from cache if available. 
if(file.exists("data/limma_model.RData")) {
  load("data/limma_model.RData")
} else {
  # RIN = RNA quality index
  # SMTS = tissue
  dgList = DGEList(counts=gct, genes=rownames(gct))
  designMat = model.matrix(~colData$SMTS + colData$RIN + colData$Gender)
  
  # filter for cpm > 1 for at least 10 samples per gene
  CPM = cpm(dgList)
  cpmCnt = apply(CPM > 1, 1, sum)
  cpmInds = cpmCnt >= 10
  dgList.f = dgList[cpmInds,]
  
  dgList.f = calcNormFactors(dgList.f)
  
  v = voom(dgList.f, designMat, plot=FALSE)
  
  fit = lmFit(v, designMat)
  save(fit, file="limma_model.RData")
}

# apply cutoffs. 
fit.e = treat(fit, lfc=log2(100))
rslt = decideTests(fit.e, p.value=(0.01 / nrow(fit)))
~~~~

### Annotate Gene Symbols {#annotate-gene-symbols}

As the *BioQC* signatures use gene symbols, we need to convert the
ENSEMBL gene ids to HGNC gene symbols. We achieve this using the
[*biomaRt*](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
Bioconductor package.

~~~~ r
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")  
ensemble2hgnc = data.table(getBM(mart=ensembl, attributes=(c("ensembl_gene_id", "hgnc_symbol"))))
ensemble2hgnc = ensemble2hgnc[hgnc_symbol!='',]

# remove the .digit from the gene name. 
salitizeEnid = function(x) { 
  return(str_split(x, "\\.")[[1]][1])
}

# append sanitized EnsemblID to data table 
rslt.tab = data.table(rslt, keep.rownames=TRUE)
names = make.names(colnames(rslt.tab))
setnames(rslt.tab, names)
rn.2 = lapply(rslt.tab[,rn], salitizeEnid)
rslt.tab = rslt.tab[,rn2:=rn.2]

# append hgnc to data table 
rslt.matched = data.table(matchColumn(ensemble2hgnc$ensembl_gene_id, rslt.tab, 'rn2'))
rslt.matched = rslt.matched[,hgnc:=ensemble2hgnc[,hgnc_symbol]]
~~~~

### Create signatures {#create-signatures}

Finally, we create lists of signatures for each tissue by collating the
list of gene symbols that are marked as overexpressed by *limma*.

~~~~ r
signatures = list()
sigCols =  names[grepl("^colData", names)]
sigCols = sigCols[1:(length(sigCols)-2)] # remove RIN and Gender
for(tCol in sigCols) {
  signatures = append(signatures, list(rslt.matched[get(tCol)==1,hgnc]))
}

names(signatures) = str_match(sigCols, "colData.SMTS(.*)")[,2]
~~~~

For comparison, we create Gini-index-based signatures on the same
dataset. The gini-index for each gene has already been calculated in
`data/sigmat.all.gini`

~~~~ r
gini = data.table(read_tsv("data/gtex/sigmat.all.gini"))
~~~~

    ## Parsed with column specification:
    ## cols(
    ##   NAME = col_character(),
    ##   CATEGORY = col_character(),
    ##   VALUE = col_double(),
    ##   RANKING = col_integer(),
    ##   GINI_IDX = col_double()
    ## )

~~~~ r
gini = gini[RANKING <= 3,]
gini = gini[GINI_IDX >= .7,]
gini = gini[VALUE > 1,] # cpm filter
gini = gini[,CATEGORY:=str_match(CATEGORY, "(.*)\\.gct")[,2]]
signatures.gini = list()
sigCols.gini = unique(gini[,CATEGORY])
for(tCol in sigCols.gini) {
  signatures.gini = append(signatures.gini, list(gini[CATEGORY==tCol,NAME]))
}
names(signatures.gini) = sigCols.gini
~~~~

### Export to gmt {#export-to-gmt}

R Session Info {#r-session-info}
--------------

~~~~ r
sessionInfo()
~~~~

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] methods   stats     graphics  grDevices utils     datasets  base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.13        stringr_1.0.0     data.table_1.9.6  edgeR_3.14.0     
    ## [5] limma_3.28.14     biomaRt_2.28.0    ribiosIO_1.0-40   readr_1.0.0      
    ## [9] ribiosUtils_1.1-1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.5          AnnotationDbi_1.34.4 magrittr_1.5        
    ##  [4] IRanges_2.6.1        BiocGenerics_0.18.0  tools_3.3.1         
    ##  [7] parallel_3.3.1       Biobase_2.32.0       DBI_0.4-1           
    ## [10] htmltools_0.3.5      yaml_2.1.13          assertthat_0.1      
    ## [13] digest_0.6.9         tibble_1.1           formatR_1.4         
    ## [16] S4Vectors_0.10.2     bitops_1.0-6         RCurl_1.95-4.8      
    ## [19] evaluate_0.9         RSQLite_1.0.0        rmarkdown_1.0       
    ## [22] stringi_1.1.1        BiocInstaller_1.22.3 stats4_3.3.1        
    ## [25] XML_3.98-1.4         chron_2.3-47

[^1]: Lonsdale, J., Thomas, J., Salvatore, M., Phillips, R., Lo, E.,
    Shad, S., … Moore, H. F. (2013). The Genotype-Tissue Expression
    (GTEx) project. Nature Genetics, 45(6), 580–585.
    <http://doi.org/10.1038/ng.2653>
