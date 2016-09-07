# BioQC-kidney: The kidney expression example
Supplementary Information for "Detect issue heterogenity in gene expression data with [*BioQC*](https://github.com/Accio/BioQC)" ([Jitao David Zhang](mailto:jitao_david.zhang@roche.com), Klas Hatje, Clemens Broger, Martin Ebeling and [Laura Badi](laura.badi@roche.com))




In this vignette, we perform simulations with both model-generated and real-world data using *BioQC*. We show that *BioQC* is a fast and sensitive method to detect tissue heterogeneity from high-throughput gene expression data. The source code to produce this document can be found in the github repository[BioQC-example](https://github.com/Accio/BioQC-example).

*BioQC* is a R/Bioconductor package to detect tissue heterogeneity from high-throughput gene expression profiling data. It implements an      efficient Wilcoxon-Mann-Whitney test, and offers tissue-specific gene signatures that are ready to use 'out of the box'.


Experiment setup
----------------
In this document, we perform three simulation studies with *BioQC*:

* **Time benchmark** tests the time-efficiency of the Wilcoxon test implemented in *BioQC*, compared with the native implementation in *R;
* **Sensitivity benchmark** tests the sensitivity and specificity of *BioQC* detecting tissue heterogeneity using model-generated simulated data;
* **Mixing benchmark** tests the sensitivity and specificity of *BioQC* using simulated contamination with real-world data.

