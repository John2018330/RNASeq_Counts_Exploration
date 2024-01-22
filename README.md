# RNASeq Counts Exploration- a RShiny applicaiton
BF591 Final Project

Created by Johnathan Zhang

### Project Overview
This RShiny application is my final project for the class "BF591 R for Biological Sciences". The purpose of the app allows the user to explore a differential expression dataset comparing post-mortem Huntington’s Disease prefrontal cortex tissue and neurologically healthy tissue. By feeding in different data files from the `/data/` folder into the tabs of the application, one can explore, visualize, and analyze the data to gain insights into the differential expression data.

The input data files include:
- `GSE64810_series_matrix.txt` - a text file containing metadata of the samples, publication and author notes, sequencing platform information, etc...

- `GSE64810_raw_counts.tsv` - a tsv of raw counts generated from RNASeq; samples as columns and genes as rows

- `GSE64810_norm_counts_TPM.tsv` - a tsv of the raw counts normalized by transcripts per million (TPM)

- `GSE64810_norm_counts_FPKM.tsv` - a tsv of the raw counts normalized by fragments per kilobase of transcript per million reads (FPKM)

- `GSE64810_DESeq2.txt` - a text file containing the results of running DESeq2 on normalized counts

- `h.all.v2023.2.Hs.symbols.gmt` - a file containing the hallmark gene sets to use in gene set enrichment analysis

### Useful Links
The data used for this project can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810) and the original paper is [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4670106/), titled "mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals".

Analyses used include [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html).


## Tab Descriptions
### Tab 1: Sample Metadata
After uploading the series matrix metadata file, the app selects for columns only relevant to the differential expression dataset. The app outputs three tabs that allows one to view and explore sample information.

#### Outputs
1) Summary information about each of the selected columns, such as mean & standard deviation for numeric columns, and factors representing experimental groups. Note that many of the numeric columns (onset age, duration, etc...) are specific to samples from Huntington's Disease patients.
2) A sortable & searchable table consisting of the selected columns.
3) Density plots of the selected numeric column. For variables that are shared by both diagnoses, density curves are constructed independently and then overlapped.

### Tab 2: Counts Exploration
> This tab requires the uploaded series matrix file in Tab 1 and one of the normalized counts matrix files to be uploaded. 

The purpose of the counts exploration tab is to see the structure and variance of counts across samples and genes, and to provide two metrics for filtering through the data. These metrics are (1) percentile variance, i.e. how variable the gene counts are across all samples, and (2) a non-zero count, i.e. how many samples must have a count greater than 0 for a gene. The user may use the sliders to change the threshold for these metrics and assess the effects of the filtering in the output tabs. 

#### Outputs
1) A simple table showing the numbers of genes and samples filtered out based on the selected threshold values in the input tab
2) Two scatter plots showing the distribution of the ranked median count for each gene versus variance across samples on a logarithmic scale (left) and number of samples with zero counts for the gene (right). The median counts of all genes are *ranked* (opposed to plotting their actual values) to generate a more even distribution across the x-axis. The scatter plots are also reactively color-coded based on the input filtering thresholds.
3) A heatmap of the genes and samples that pass the input filtering thresholds with clustering for samples. I am unclear of how to interpret this heatmap.
4) A PCA scatterplot of samples where the user may choose which principal components to plot against each other and see their percent variance explained. Note that **ALL** genes were used in the PCA, not the filtered ones, as was specified in the instructions of the project. This leads to no clear grouping between Huntington's disease and neurological normal diagnoses being seen. I believe using the filtered genes and samples may yield a more interesting scatterplot and will add this feature in the future.

### Tab 3: Differential Expression
Differential gene expression analysis was performed on the normalized counts using the bioconductor package DESeq2, and the results are included in the `data` folder. In this tab, upload the DESeq2 results file to view a sortable & searchable table of the results and a customizable volcano plot. By default, *log 2 fold change* is plotted against the *adjusted p-values*, and the threshold for coloring genes by p-adjusted values is 1e-25. These values and colors are changeable and the plot is updated once the plot button is clicked.


### Tab 4: Gene Set Enrichment Analysis
> This tab requires the DESeq results from the previous tab and a `.gmt` genesets file

Using the 'fgsea' bioconductor package, we can perform gene set enrichment analysis on the differential expression dataset to identify groups of genes that have been positively or negatively enriched in Huntington's disease diagnoses. The hallmark genesets are provided in a `.gmt` file in the `data` folder, though any collection of genesets should be able to be used. Using the provided sliders in the output tabs, the app will filter out or color points based on the adjusted p-values threshold set. 

#### Output
1) A horizontal barplot of genesets that pass the adjusted p-value filtering and their NES value
2) A sortable and searchable table of the fgsea results. The table can also be downloaded as a csv.
3) A scatterplot of geneset NES values against -log10 adjusted p-values. Adjust the slider to change point colors. 

