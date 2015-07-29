##CTCF and RNA splicing

### Background:
DNA -> RNA -> protein
DNA(gene) -> RNA primary transcript -> RNA mature transcript -> protein.
Splicing is the process of converting a primary to mature RNA transcript.
Each gene consists of two categories of sequences, Exons, which are kept in the mature transcript, and Introns which are disregarded. 


### Objective
To determine if the presence of CTCF protein near a gene changes the RNA structure transcribed from the gene. 

###  Approach:
CTCF binding, and RNA structure changes can be quantified from ChIP-Seq and RNA-seq data. For humans, we will look at the correlation between these two measures over all the cell types where both experiments are done. Since there are thousands of genes, the p-value cutoff for statistical significance must be adjusted using FDR analysis. The list of significant results represent genes where binding of CTCF is associated with changes in RNA structure.   

###  Data
The two primary data sources are RNA-seq and CTCF ChIP-seq. To quantify changes in RNA structure of a gene, there are two measurements collected for each exon of a gene: PSI,which measures the percent of exon inclusion in mature transcripts, and CoSI, which is a measure of how efficiently splicing occurs. To measure protein binding, near the exon, the strength of downstream CTCF peaks in measured from ChIP-seq data. In total, there are 15 cell types that contain both data sources.  

###  Analysis Pipeline
There are two main parts to this analysis: Downloading and processing the data, then running statistical tests on the refined data. To download and process the data, lists of ChIP-seq and RNA-seq experiments must be compared to find cell types with both experiments. For the selected cell types, PSI and CoSI values are found online for each exon. The downstream region for each exon is calculated. Chip_Seq data is downloaded to the HPC, processed into peaks with MACS2, then intersected with regions downstream regions of the exons using bed tools. Finally, each exon has an array of CoSI or PSI measurements over celltypes, and peak Fold change over celltypes. At this point, collecting the data is done, we have all the exons with all their measurements over the cell types. 
		Next, we will look at the correlation between PSI and CTCF binding peak strength.  Of 165,352 total exons identified by PSI, only 7,251 have a correlation value, and 194 have a correlation p-value below the FDR correction. ".  So what's going on here? NA, or missing values, are dominating the PSI measurements, so very few exons have the required 3 points to compute correlation. However, 194 exons pass the FDR threshold. This is great...untill we look a little closer and notice that most of the points that pass the FDR cutoff have all PSI values of 1, with one value less than one, and CTCF peak fold expression 0, which one value greater than 1. This suggests the results are an artifact of the distribution. To test this, we shuffled downstream exon regions and re-ran the analysis: This time there were 112 with signification p-value. To reduce this nasty tendency, I treated peak fold change values of 0 as NA. The FDR analysis contained only one point passing a FDR cutoff of 0.25, SEMA-3B. Thus, when numerical properties of the distribution are considered, there apears to be no significant correlation between exon inclusion (PSI) and CTCF peakk binding.  


###  Conclusions
CTCF and exon measurements don't seem to show correlation after artifacts are removed. This is not enough to rule out this effect from happening, as many cell yps were not tested, and deeper sequence coverage may be required to eulicate the trend. 
[FDR plot for PSI vs. CTCF peak, all artifacts removed](https://www.dropbox.com/s/slco6xm6tkkbije/PSI-noPeakNa-pValVsCutoff-line.png)

###  Lessons Learned
- When statistical measurements are not available of a majority of points and center around a single value for those that do, artificial correlation may be observed. To determine if your observation is determined by the numerical distribution, not some underlying effect, perform a shuffle test where the distributions for the variables will remain the same, but lose meaningful interpretation.  
- While still in the planning phase, knowledge of a measures highly skewed distribution or overwhelming amount of missing data can lead suitability tests and research into algorithms designed to handle non-normal data sources with missing values.  
- FDR corrections are simple to make, and greatly improve the reliability of your findings when conducting many experimental tests. 

