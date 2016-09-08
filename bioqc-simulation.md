# BioQC-benchmark: Testing Efficiency, Sensitivity and Specificity of BioQC on simulated and real-world data
Supplementary Information for "Detect issue heterogenity in gene expression data with [*BioQC*](https://github.com/Accio/BioQC)" ([Jitao David Zhang](mailto:jitao_david.zhang@roche.com), Klas Hatje, Clemens Broger, Martin Ebeling and [Laura Badi](laura.badi@roche.com))




In this vignette, we perform simulations with both model-generated and real-world data using *BioQC*. We show that *BioQC* is a fast and sensitive method to detect tissue heterogeneity from high-throughput gene expression data. The source code to produce this document can be found in the github repository [BioQC-example](https://github.com/Accio/BioQC-example).

*BioQC* is a R/Bioconductor package to detect tissue heterogeneity from high-throughput gene expression profiling data. It implements an      efficient Wilcoxon-Mann-Whitney test, and offers tissue-specific gene signatures that are ready to use 'out of the box'.


Experiment setup
----------------
In this document, we perform three simulation studies with *BioQC*:

* **Time benchmark** tests the time-efficiency of the Wilcoxon test implemented in *BioQC*, compared with the native implementation in *R;
* **Sensitivity benchmark** tests the sensitivity and specificity of *BioQC* detecting tissue heterogeneity using model-generated simulated data;
* **Mixing benchmark** tests the sensitivity and specificity of *BioQC* using simulated contamination with real-world data.

All source code that is needed to reproduce the results can be found in the `.Rmd` file generating this document. 


```r
library(testthat)
library(BioQC)
library(hgu133plus2.db) ## to simulate an microarray expression dataset
library(lattice)
library(latticeExtra)
library(GEOquery)
library(xtable)
library(gplots)
library(rbenchmark)

pdf.options(family="ArialMT", useDingbats=FALSE)

set.seed(1887)

## list human genes
humanGenes <- unique(na.omit(unlist(as.list(hgu133plus2SYMBOL))))

## read tissue-specific gene signatures
gmtFile <- system.file("extdata/exp.tissuemark.affy.roche.symbols.gmt",
                       package="BioQC")
gmt <- readGmt(gmtFile)
```


Time benchmark
--------------
In the first experiment, we setup expression matrices of 20155 human protein-coding genes of 1, 5, 10, 50, or 100 samples. Genes  are $i.i.d$ distributed following $\mathcal{N}(0,1)$. The Wilcoxon-Mann-Whitney test implemented in *BioQC* and the native *R* implementation are applied to the matrices respectively.




The numeric results of both implementations, `bioqcNumRes` (from *BioQC*) and `rNumRes` (from *R*), are equivalent, as shown by the next     command.


```r
expect_equal(bioqcNumRes, rNumRes)
```

The *BioQC* implementation is more than 500 times much faster: while it takes about one second for BioQC to calculate enrichment scores of all 155 signatures in 100 samples, the native R implementation takes about 20 minutes: 



<div class="figure" style="text-align: center">
<img src="bioqc-simulation_files/figure-html/time_benchmark_vis-1.svg" alt="Time benchmark results of BioQC and R implementation of Wilcoxon-Mann-Whitney test. Left panel: elapsed time in seconds (logarithmic Y-axis). Right panel: ratio of elapsed time by two implementations. All results achieved by a single thread on in a RedHat Linux server." style="display:block; margin: auto" />
<p class="caption">Time benchmark results of BioQC and R implementation of Wilcoxon-Mann-Whitney test. Left panel: elapsed time in seconds (logarithmic Y-axis). Right panel: ratio of elapsed time by two implementations. All results achieved by a single thread on in a RedHat Linux server.</p>
</div>

The main reason underlying the low performance of R implementation is that the `wilcox.test` function sorts two numeric vectors that are to be compared. When the function is repeatedly applied to gene expression data, it performs many expensive sorting operations which are futile, because the sorting of genes outside of the gene set (*background genes*) does not change between samples. *BioQC* sorts the background genes     only once for each gene set, independent of how many samples are tested.

In addition, *BioQC* implements an approximate Wilcoxon test instead of the exact version, because the difference between the two is          negligible for high-throughput gene expression data. Last but not least, *BioQC* implements its core algorithm in C so as to maximize the     time- and memory-efficiency.

Putting these tweaks together, *BioQC* achieves identical results as the native implementation with two order of magnitude less time. This    renders *BioQC*~an highly efficient tool for quality control of large-scale high-throughput gene expression data.



Sensitivity Benchmark
---------------------
We next asked the question how sensitive is *BioQC* to expression changes of tissue signature genes. Similar to the previous simulation,      while keeping all other genes $i.i.d$ normally distributed following $\mathcal{N}(0,1)$, we dedicatedly increase the expression of genes in one randomly selected tissue signature (ovary, with 43 genes) by different amplitudes: these genes' expression levels are randomly drawn from different normal distributions with varying expectation and constant variance between $\mathcal{N}(0,1)$ and $\mathcal{N}(3,1)$. To test the robustness of the algorithm, 10 samples are generated for each mean expression difference value.

<div class="figure" style="text-align: center">
<img src="bioqc-simulation_files/figure-html/sensitivity_benchmark_fig-1.svg" alt="Sensitivity benchmark. Expression levels of genes in the ovary signature are dedicately sampled randomly from normal distributions with different mean values. Left panel: enrichment scores reported by *BioQC* for the ovary signature, plotted against the differences in mean expression values; Right panel: rank of ovary enrichment scores in all `r length(gmt)` signatures plotted against the difference in mean expression values." style="display:block; margin: auto" />
<p class="caption">Sensitivity benchmark. Expression levels of genes in the ovary signature are dedicately sampled randomly from normal distributions with different mean values. Left panel: enrichment scores reported by *BioQC* for the ovary signature, plotted against the differences in mean expression values; Right panel: rank of ovary enrichment scores in all `r length(gmt)` signatures plotted against the difference in mean expression values.</p>
</div>


The above figure visualizes the distribution of enrichment scores and their ranks dependent on the mean expression difference between ovary signature genes and background genes. As soon as the expression of signature genes increases by a very moderate ampltiude $1\sigma$, *BioQC* will identify the gene set as the highest-ranking signature. A even stronger difference in expression will lead to higher enrichment scores but no change in the rank.

The results suggest that *BioQC* is sensitive even to moderate changes in the average expression of a gene set.




Mixing Benchmark
----------------
The sensitivity benchmark above suffers from the limitation that the distributions of gene expression are not physiological. To overcome this, we      designed and performed a benchmark by *in silico* mixing expression profiles with weighted linear combination, thereby mimicking tissue contamination.

Given the expression profile of a sample of tissue A ($\mathbf{Y}_A$), and that of a sample of tissue B ($\mathbf{Y}_B$), the weighted linear mixing produces a new profile $\mathbf{Y}=\omega\mathbf{Y_A}+(1-\omega)\mathbf{Y_B}$, where $\omega \in [0,1]$. In essence the profiles of two tissue types are linearly mixed in different proportions, which simulates varying severities of contaminations. We asked whether BioQC could detect such mixings, and if so, how sensitive is the method.





### Dataset selection and quality control

In order to avoid over-fitting of signatures derived from human expression data, we decided to use a normal tissue expression dataset from a non-      human mammal species, because it has been shown that tissue-preferential expression patterns tend to be conserved between mammal species. We           identified a dataset of *Canis lupus familiaris* (dog), which is publicly available in Gene Expression Omnibus ([GDS4164](http://www.ncbi.nlm.nih.  gov/sites/GDSbrowser?acc=GDS4164)).

In this study, the authors examined 39 samples from 10 pathologically normal tissues (liver, kidney, heart, lung, brain, lymph node, spleen, jejunum,  pancreas, and skeletal muscle) of four dogs (with one pancreas sample missing). We downloaded the data, and performed minimal pre-processing: for multiple probesets that map to same genes, we kept the one with the highest average expression level and removed the rest. The resulting dataset contained expression of 16797 genes. BioQC was applied to the dataset to test whether there are major contamination issues. The results, including tissues reported by the authors, and the BioQC tissue signatures with the highest and second-highest scores, are reported in the following table: 


Table: Quality control of the mixing benchmark input data with *BioQC*. Four columns (f.l.t.r.): sample index; tissue         reported by the authors; the tissue signature with the highest enrichment score reported by *BioQC*; the tissue signature with the second-   highest enrichment score.

            Label            BioQC.best                BioQC.second        
----------  ---------------  ------------------------  --------------------
GSM502573   Brain            Spinal_cord               Nodose_nucleus      
GSM502574   Brain            Brain_Cortex_prefrontal   Brain_Amygdala      
GSM502575   Brain            Brain_Cortex_prefrontal   Brain_Amygdala      
GSM502576   Brain            Brain_Cortex_prefrontal   Brain_Amygdala      
GSM502577   Heart            Muscle_cardiac            Muscle_skeletal     
GSM502578   Heart            Muscle_cardiac            Muscle_skeletal     
GSM502579   Heart            Muscle_cardiac            Muscle_skeletal     
GSM502580   Heart            Muscle_cardiac            Muscle_skeletal     
GSM502581   Jejunum          Intestine_small           Intestine_Colon     
GSM502582   Jejunum          Intestine_small           Intestine_Colon     
GSM502583   Jejunum          Intestine_small           Intestine_Colon     
GSM502584   Jejunum          Intestine_small           Intestine_Colon     
GSM502585   Kidney           Kidney                    Kidney_Renal_Cortex 
GSM502586   Kidney           Kidney                    Kidney_Renal_Cortex 
GSM502587   Kidney           Kidney                    Kidney_Renal_Cortex 
GSM502588   Kidney           Kidney                    Kidney_Renal_Cortex 
GSM502589   Liver            Liver                     Liver               
GSM502590   Liver            Liver                     Liver               
GSM502591   Liver            Liver                     Liver               
GSM502592   Liver            Liver                     Liver               
GSM502593   Lung             Lung                      Monocytes           
GSM502594   Lung             Monocytes                 Lung                
GSM502595   Lung             Lung                      Monocytes           
GSM502596   Lung             Monocytes                 Lung                
GSM502597   LymphNode        Lymphocyte_B_FOLL         Lymphocytes_T_H     
GSM502598   LymphNode        Lymphocyte_B_FOLL         Lymphocytes_T_H     
GSM502599   LymphNode        Lymphocyte_B_FOLL         Lymphocytes_T_H     
GSM502600   LymphNode        Lymphocyte_B_FOLL         Lymphocytes_T_H     
GSM502601   Pancreas         Pancreas_Islets           Pancreas            
GSM502602   Pancreas         Pancreas_Islets           Pancreas            
GSM502603   Pancreas         Pancreas_Islets           Pancreas            
GSM502604   SkeletalMuscle   Muscle_skeletal           Muscle_cardiac      
GSM502605   SkeletalMuscle   Muscle_skeletal           Muscle_cardiac      
GSM502606   SkeletalMuscle   Muscle_skeletal           Muscle_cardiac      
GSM502607   SkeletalMuscle   Muscle_skeletal           Muscle_cardiac      
GSM502608   Spleen           Monocytes                 Lymphocyte_B_FOLL   
GSM502609   Spleen           Monocytes                 Lymphocyte_B_FOLL   
GSM502610   Spleen           Monocytes                 Erythroid_cells     
GSM502611   Spleen           Monocytes                 Myeloblast          

By comparing the tissue labels provided by the authors and the predictions of *BioQC*, we notice that in most cases the two match well        (despite of ontological differences). In three cases (sample ID GSM502573, GSM502594, and GSM502596) though there seem to be intriguing differences, which might be explained by different sampling procedures or immune cell infiltration. We will however in this vignette not further explore them. These three samples are removed from the simulation procedures.



R Session Info
----------------

```r
sessionInfo()
```

```
## R version 3.1.3 (2015-03-09)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux Server release 6.3 (Santiago)
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
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] rbenchmark_1.0.0     gplots_3.0.1         xtable_1.8-2        
##  [4] GEOquery_2.32.0      latticeExtra_0.6-28  RColorBrewer_1.1-2  
##  [7] lattice_0.20-33      hgu133plus2.db_3.0.0 org.Hs.eg.db_3.0.0  
## [10] RSQLite_1.0.0        DBI_0.4-1            AnnotationDbi_1.28.2
## [13] GenomeInfoDb_1.2.5   IRanges_2.0.1        S4Vectors_0.4.0     
## [16] BioQC_1.02.1         Biobase_2.26.0       BiocGenerics_0.12.1 
## [19] Rcpp_0.12.0          testthat_1.0.2       knitr_1.13          
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6       caTools_1.17.1     crayon_1.3.1      
##  [4] digest_0.6.9       evaluate_0.9       formatR_1.4       
##  [7] gdata_2.17.0       grid_3.1.3         gtools_3.5.0      
## [10] highr_0.6          htmltools_0.3.5    KernSmooth_2.23-15
## [13] magrittr_1.5       memoise_1.0.0      R6_2.1.2          
## [16] RCurl_1.95-4.8     rmarkdown_1.0      stringi_1.0-1     
## [19] stringr_1.0.0      tools_3.1.3        XML_3.98-1.3      
## [22] yaml_2.1.13
```
