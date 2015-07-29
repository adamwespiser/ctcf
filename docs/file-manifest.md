### File Manifest
I am the sole contributer of analysis and code in this project

`ctcf/analysis` contains the R code to for the entire pipeline.    
`enhancer_pred/data`  contains the data files for exons, alternative splicing measures, and CTCF binding metrics.    
`ctcf/lots` has plots for explority analysis between CTCF binding and measures of exon splicing    

### /ctcf/analysis
`main.R` contains information about where resources are located, a set of helper functions to facilitate input/output, source statements for R program files, and the R packages needed to run this analysis.
`exprPeak.R` is responsible for loading downloaded and cleaned data, applying statistical tests, then plotting the results. 
`/downloadData/getENCODE_ChIP_Seq.R` handles all data preprocessing. First, a set of experiments is determined from ENCODE metadata. For ChIP-seq data, these experiments are consolidated into a download script, which is transfered to our High Performance Computing Cluster(HPC), then executed. Next, the ChIP-seq analysis script is created, transfered, then queued into the job scheduler on the HPC. The results are transfered back to my local server. The exon splicing statistics are directly downloaded from an authors site, and transformed into usable form. Given the exons locations, their downstream regions can be intersected with CTCF binding sites using bedtools. Again, a meta-programming approach is taken, wherein a script is created as an R string, written to file, then executed. Finally, the downstream CTCF peak values are combined with the exon splicing measures for all cell types, giving the final processed data.     


### /ctcf/plots
`CoSI-FoldChange` directory contains plots about the relationship between CoSI, a measure of how complete a primary transcript is spliced, and CTCF binding. 
`PSI-FoldChange` directory contains plots about the relationship between PSI, a measure of exon retention in mature transcript, and CTCF binding
