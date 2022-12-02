# Final Project Outline

## Title
Differential Gene Expression in TCGA within stage 1-4 Ovarian Serous Cystadenocarcinoma comparing radiation treatment vs. non-radiation treatment using DeSEQ2

## Author
Shanna Badhesha

## Overview of the Project
I will identify differentially expressed genes between Ovarian Serous Cystadenocarcinoma radiation treatment using carboplatin vs. non-radiation treatment without carboplatin. This analysis will utilize the package DeSEQ2 and the http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html vignette. For this analysis, I will be using the TCGA-OV cohort. I have identified a total of 40 STAR-counts files with 20 radiation treatment and 20 non-radiation treatment files. 

## Data
I will use the data from the following website: https://portal.gdc.cancer.gov/repository. There are 585 STAR-count files with 450 radiation treatment and 135 non-radiation treatment files in total. I will be using 40 STAR-counts files with 20 samples per group. The specific files are available [here](https://portal.gdc.cancer.gov/repository?facetTab=files&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22ovary%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D).

## Milestone 1
I have downloaded the data and matched the file names to the corresponding TCGA IDs. I have isolated the necessary unstranded column from each file and merged them together for the analysis. I have created another spreadsheet identifying the TCGA IDs for each of the two treatment groups.

## Milestone 2
I will do an initial run of the entire vignette listed above for analysis. This will include having data loaded in (completed during milestone 1) and run through the entire vignette. I will be seeking feedback during this step. I will begin the process of uploading documentation and defining my results.

## Deliverable
Below is the completed vignette with my 40 samples. 

## Organzing Sample Files

The STAR-Count files provide a lot of information that I will not be needing for the analysis. For the purpose of this vignette, we will extract the unstranded gene count column from each file and merge them together.

Create a ```gene.id``` file that includes all the gene id's from one of the sample files: 
	```awk '{print $1}’ TCGA-04-1331.tsv > TCGA-04-1331.txt```

Extract the unstranded column (column 4) from each tsv file:
	```awk '{print $4}’ TCGA-04-1331.tsv > TCGA-04-1331.txt```

Remove the first 6 lines from each tsv file since they will not be needed for our analysis:
	```sed -i .bak '1,6d' 'TCGA-04-1332.txt'```

Add a header to each of the file using ```vi```. The header should be the sample ID. 

Merge the files together using the ```paste``` command: 
	```paste gene.id *.txt > merged_files.txt```

### Sample Information Table 

Create a .txt file using Microsoft Excel with the Sample IDs and their condition. The Sample IDs should be in the same order as the merged_file.txt otherwise DESeq2 will not accept the file. 

## DESeq2 Vignette

### Count matrix input
Load in DESeq2. The function DESeqDataSetFromMatrix can be used if you already have a matrix of read counts prepared from another source. Another method for quickly producing count matrices from alignment files is the featureCounts function (Liao, Smyth, and Shi 2013) in the Rsubread package. To use DESeqDataSetFromMatrix, the user should provide the counts matrix, the information about the samples (the columns of the count matrix) as a DataFrame or data.frame, and the design formula.

```r
library("DESeq2")


file <- read.delim(file= "Desktop/Final_project/merged_files.txt", header = TRUE, sep = "\t", row.names = "X")
cts <- as.matrix(file)                 
coldata <- read.delim(file= "Desktop/Final_project/TCGA_groups.txt", header = TRUE, sep = "\t", row.names=1) 
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[,rownames(coldata)]
all(rownames(coldata) == colnames(cts))
```

We examine the count matrix and column data to see if they are consistent in terms of sample order.

```r
head(cts,2)
```
                   TCGA.04.1331 TCGA.04.1332 TCGA.04.1338 TCGA.04.1341 TCGA.04.1343 TCGA.04.1347 TCGA.04.1350
ENSG00000000003.15         6932         2495         1262         3957         3726         5462         4687
ENSG00000000005.6            10           12            1            8          185           13            0
                   TCGA.04.1356 TCGA.04.1357 TCGA.04.1361 TCGA.04.1362 TCGA.04.1364 TCGA.04.1365 TCGA.04.1514
ENSG00000000003.15         5423         2130         7111        12071         5130         5787        18147
ENSG00000000005.6             2            3           42            6            5            6          158
                   TCGA.04.1519 TCGA.04.1530 TCGA.04.1536 TCGA.04.1542 TCGA.04.1648 TCGA.04.1651 TCGA.04.1655
ENSG00000000003.15         7414         1519         1401         3181         9861         8024         6103
ENSG00000000005.6             1            9            7            9            5           49           55
                   TCGA.09.0364 TCGA.09.0366 TCGA.09.1659 TCGA.09.2048 TCGA.09.2056 TCGA.10.0928 TCGA.13.0730
ENSG00000000003.15         6595        13998         5339        13953         3810         4327         3218
ENSG00000000005.6            15           16            6           14            0           14         1244
                   TCGA.13.0762 TCGA.13.0765 TCGA.13.0766 TCGA.13.0800 TCGA.13.0916 TCGA.13.0920 TCGA.13.0923
ENSG00000000003.15         5838         4273         2143         5905        15855         3566         2008
ENSG00000000005.6             1            1            5           68          823            7            2
                   TCGA.13.1403 TCGA.13.1505 TCGA.13.1507 TCGA.13.2060 TCGA.23.1029
ENSG00000000003.15         2875         5474         2561         6388         9557
ENSG00000000005.6             1           49            9           37          153

```coldata```

             condition
TCGA.04.1331       rad
TCGA.04.1332       rad
TCGA.04.1338       rad
TCGA.04.1341   non-rad
TCGA.04.1343       rad
TCGA.04.1347       rad
TCGA.04.1350       rad
TCGA.04.1356       rad
TCGA.04.1357   non-rad
TCGA.04.1361       rad
TCGA.04.1362       rad
TCGA.04.1364       rad
TCGA.04.1365       rad
TCGA.04.1514       rad
TCGA.04.1519   non-rad
TCGA.04.1530       rad
TCGA.04.1536       rad
TCGA.04.1542       rad
TCGA.04.1648       rad
TCGA.04.1651       rad
TCGA.04.1655       rad
TCGA.09.0364       rad
TCGA.09.0366       rad
TCGA.09.1659   non-rad
TCGA.09.2048   non-rad
TCGA.09.2056   non-rad
TCGA.10.0928   non-rad
TCGA.13.0730   non-rad
TCGA.13.0762   non-rad
TCGA.13.0765   non-rad
TCGA.13.0766   non-rad
TCGA.13.0800   non-rad
TCGA.13.0916   non-rad
TCGA.13.0920   non-rad
TCGA.13.0923   non-rad
TCGA.13.1403   non-rad
TCGA.13.1505   non-rad
TCGA.13.1507   non-rad
TCGA.13.2060   non-rad
TCGA.23.1029   non-rad

### Building the DESeqDataSet
With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:

```r
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

`dds
```

class: DESeqDataSet 
dim: 60660 40 
metadata(1): version
assays(1): counts
rownames(60660): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288674.1 ENSG00000288675.1
rowData names(0):
colnames(40): TCGA.04.1331 TCGA.04.1332 ... TCGA.13.2060 TCGA.23.1029
colData names(1): condition

### Pre-Filtering
While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted.

Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total. Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function.

```r
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

### Factor
By default, R will choose a reference level for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the levels. There are two solutions: you can either explicitly tell results which comparison to make using the contrast argument (this will be shown later), or you can explicitly set the factors levels. In order to see the change of reference levels reflected in the results names, you need to either run DESeq or nbinomWaldTest/nbinomLRT after the re-leveling operation. 

```r
dds$condition <- factor(dds$condition, levels = c("rad","non-rad"))
```

### Differential expression analysis
The standard differential expression analysis steps are wrapped into a single function, DESeq. The estimation steps performed by this function are described below, in the manual page for ```?DESeq``` and in the Methods section of the DESeq2 publication (Love, Huber, and Anders 2014).

Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values. With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level (see previous note on factor levels). However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.

Details about the comparison are printed to the console, directly above the results table. The text, ```condition non-rad vs. rad```, tells you that the estimates are of the logarithmic fold change log2(rad/non-rad).

```r
dds <- DESeq(dds)
res <- results(dds)
res
```

log2 fold change (MLE): condition non.rad vs rad 
Wald test p-value: condition non.rad vs rad 
DataFrame with 45641 rows and 6 columns
                     baseMean log2FoldChange     lfcSE       stat    pvalue      padj
                    <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
ENSG00000000003.15  5253.9122     -0.0505555  0.196991 -0.2566384 0.7974579  0.945787
ENSG00000000005.6     25.2649      0.9278340  0.608227  1.5254727 0.1271412  0.536655
ENSG00000000419.13  3237.4165      0.0142636  0.192966  0.0739176 0.9410760  0.984920
ENSG00000000457.14   676.9305      0.0751432  0.164663  0.4563460 0.6481412  0.891204
ENSG00000000460.17   457.1319      0.4220209  0.217847  1.9372332 0.0527168  0.393108
...                       ...            ...       ...        ...       ...       ...
ENSG00000288663.1   54.255850       0.122303  0.243648   0.501965 0.6156920  0.876257
ENSG00000288667.1    0.805876      -1.523881  0.743250  -2.050295 0.0403357        NA
ENSG00000288670.1  272.987292      -0.167227  0.240384  -0.695663 0.4866399  0.812998
ENSG00000288674.1   10.466399       0.499068  0.342389   1.457603 0.1449500  0.557620
ENSG00000288675.1   37.388481       0.446755  0.245943   1.816497 0.0692941  0.436438

### Log fold change shrinkage for visualization and ranking
Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the ```dds``` object to the function lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.

We provide the ```dds``` object and the name or number of the coefficient we want to shrink, where the number refers to the order of the coefficient as it appears in ```resultsNames(dds)```.

```r
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_non.rad_vs_rad", type="apeglm")
resLFC
```
log2 fold change (MAP): condition non.rad vs rad 
Wald test p-value: condition non.rad vs rad 
DataFrame with 45641 rows and 5 columns
                     baseMean log2FoldChange      lfcSE    pvalue      padj
                    <numeric>      <numeric>  <numeric> <numeric> <numeric>
ENSG00000000003.15  5253.9122   -1.34247e-06 0.00144266 0.7974579  0.945787
ENSG00000000005.6     25.2649    2.15299e+00 0.66318587 0.1271412  0.536655
ENSG00000000419.13  3237.4165    4.11597e-07 0.00144265 0.9410760  0.984920
ENSG00000000457.14   676.9305    3.87099e-06 0.00144264 0.6481412  0.891204
ENSG00000000460.17   457.1319    1.68436e-05 0.00144271 0.0527168  0.393108
...                       ...            ...        ...       ...       ...
ENSG00000288663.1   54.255850    4.92564e-07 0.00144267 0.6156920  0.876257
ENSG00000288667.1    0.805876   -2.98944e-06 0.00144269 0.0403357        NA
ENSG00000288670.1  272.987292   -8.32522e-06 0.00144268 0.4866399  0.812998
ENSG00000288674.1   10.466399    4.63382e-06 0.00144269 0.1449500  0.557620
ENSG00000288675.1   37.388481    3.58464e-06 0.00144267 0.0692941  0.436438

### Speed-up and parallelization thoughts
```r
library("BiocParallel")
register(MulticoreParam(4))
```

### p-values and adjusted p-values
```r
resOrdered <- res[order(res$pvalue),]
summary(res)
```

out of 45628 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 222, 0.49%
LFC < 0 (down)     : 179, 0.39%
outliers [1]       : 0, 0%
low counts [2]     : 15052, 33%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

```r
sum(res$padj < 0.1, na.rm=TRUE)
```
[1] 401

```r
res05 <- results(dds, alpha=0.05)
summary(res05)
```
out of 45628 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 80, 0.18%
LFC < 0 (down)     : 52, 0.11%
outliers [1]       : 0, 0%
low counts [2]     : 9744, 21%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

```r
sum(res05$padj < 0.05, na.rm=TRUE)
```
[1] 132

### Exploring and exporting results

### MA-plot
In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

```r
plotMA(res, ylim=c(-2,2))
```
GRAPH_1
```r
plotMA(resLFC, ylim=c(-2,2))
```
GRAPH_2

After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:
```r
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
```

### Alternative shrinkage estimators

The moderated log fold changes proposed by Love, Huber, and Anders (2014) use a normal prior distribution, centered on zero and with a scale that is fit to the data. The shrunken log fold changes are useful for ranking and visualization, without the need for arbitrary filters on low count genes. The normal prior can sometimes produce too strong of shrinkage for certain datasets. In DESeq2 version 1.18, we include two additional adaptive shrinkage estimators, available via the type argument of lfcShrink. For more details, see ?lfcShrink

The options for type are:

apeglm is the adaptive t prior shrinkage estimator from the apeglm package (Zhu, Ibrahim, and Love 2018). As of version 1.28.0, it is the default estimator.
ashr is the adaptive shrinkage estimator from the ashr package (Stephens 2016). Here DESeq2 uses the ashr option to fit a mixture of Normal distributions to form the prior, with method="shrinkage".
normal is the the original DESeq2 shrinkage estimator, an adaptive Normal distribution as prior.
If the shrinkage estimator apeglm is used in published research, please cite:

Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. 10.1093/bioinformatics/bty895

If the shrinkage estimator ashr is used in published research, please cite:

Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. 10.1093/biostatistics/kxw041

In the LFC shrinkage code above, we specified coef="condition_non.rad_vs_rad". We can also just specify the coefficient by the order that it appears in resultsNames(dds), in this case coef=2. For more details explaining how the shrinkage estimators differ, and what kinds of designs, contrasts and output is provided by each, see the extended section on shrinkage estimators.

```r
resultsNames(dds)
```

  ## [1] "Intercept"                      "condition_non.rad_vs_untreated"

```r
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
```

```r
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```

### Plot Counts
It can also be useful to examine the counts of reads for a single gene across the groups. A simple function for making this plot is plotCounts, which normalizes counts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables in intgroup, where more than one variable can be specified. Here we specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by rowname or by numeric index.

```r
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```

For customized plotting, an argument returnData specifies that the function should only return a data.frame for plotting with ggplot.

```r
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```

### Exporting results to CSV files 
A plain-text file of the results can be exported using the base R functions write.csv or write.delim. We suggest using a descriptive file name indicating the variable and levels which were tested.

```r
write.csv(as.data.frame(resOrdered), 
          file="condition_rad_results.csv")
```


Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, followed by the write.csv function.

```r
resSig <- subset(resOrdered, padj < 0.1)
resSig
```

### Data transformations and visualization
### Extracting transformed values
These transformation functions return an object of class DESeqTransform which is a subclass of RangedSummarizedExperiment. For ~20 samples, running on a newly created DESeqDataSet, rlog may take 30 seconds, while vst takes less than 1 second. The running times are shorter when using blind=FALSE and if the function DESeq has already been run, because then it is not necessary to re-estimate the dispersion values. The assay function is used to extract the matrix of normalized values.

```r
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
```

### Effects of transformations on the variance
The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation, the regularized log transformation and the variance stabilizing transformation. The shifted logarithm has elevated standard deviation in the lower count range, and the regularized log to a lesser extent, while for the variance stabilized data the standard deviation is roughly constant along the whole dynamic range.

Note that the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions. While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions.

```r
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
```

```r
meanSdPlot(assay(vsd))
```

```r
meanSdPlot(assay(rld))
```

### Data quality assessment by sample clustering and visualization
Data quality assessment and quality control (i.e. the removal of insufficiently good data) are essential steps of any data analysis. These steps should typically be performed very early in the analysis of a new data set, preceding or in parallel to the differential expression testing.

### Heatmap of the count matrix
To explore a count matrix, it is often instructive to look at it as a heatmap. Below we show how to produce such a heatmap for various transformations of the data.

```r
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

```r
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

```r
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

### Heatmap of the sample-to-sample distances
Another use of the transformed data is sample clustering. Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.

```r
sampleDists <- dist(t(assay(vsd)))
```

A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. We have to provide a hierarchical clustering hc to the heatmap function based on the sample distances, or else the heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix.

```r
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

### Principal component plot of the samples
Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

```r
plotPCA(vsd, intgroup=c("condition", "type"))
```

It is also possible to customize the PCA plot using the ggplot function.

```r
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```

### Likelihood ratio test
DESeq2 offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model.

The likelihood ratio test can be performed by specifying test="LRT" when using the DESeq function, and providing a reduced design formula, e.g. one in which a number of terms from design(dds) are removed. The degrees of freedom for the test is obtained from the difference between the number of parameters in the two models. A simple likelihood ratio test, if the full design was ~condition would look like:

```r
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)
res
```

### Extended section on shrinkage estimators
```r
resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
```

### Dispersion plot and fitting alternatives
Plotting the dispersion estimates is a useful diagnostic. The dispersion plot below is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the fitted value, (this outlier detection is described in the manual page for estimateDispersionsMAP). The amount of shrinkage can be more or less than seen here, depending on the sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.

```r
plotDispEsts(dds)
```

### Independent filtering of results
The results function of the DESeq2 package performs independent filtering by default using the mean of normalized counts as a filter statistic. A threshold on the filter statistic is found which optimizes the number of adjusted p values lower than a significance level alpha (we use the standard variable name for significance level, though it is unrelated to the dispersion parameter α). The theory behind independent filtering is discussed in greater detail below. The adjusted p values for the genes which do not pass the filter threshold are set to NA.

The default independent filtering is performed using the filtered_p function of the genefilter package, and all of the arguments of filtered_p can be passed to the results function. The filter threshold value and the number of rejections at each quantile of the filter statistic are available as metadata of the object returned by results.

```r
metadata(res)$alpha
```

```r
metadata(res)$filterThreshold
```

```r
plot(metadata(res)$filterNumRej, 
     type="b",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
```

```r
resNoFilt <- results(dds, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))
```

### Access to all calculated values
All row-wise calculated values (intermediate dispersion calculations, coefficients, standard errors, etc.) are stored in the DESeqDataSet object, e.g. dds in this vignette. These values are accessible by calling mcols on dds. Descriptions of the columns are accessible by two calls to mcols. Note that the call to substr below is only for display purposes.

```r
mcols(dds,use.names=TRUE)[1:4,1:4]
```

```r
substr(names(mcols(dds)),1,10) 
```
```r
mcols(mcols(dds), use.names=TRUE)[1:4,]
```
```r
head(assays(dds)[["mu"]])
```
```r
head(assays(dds)[["cooks"]])
```
```r
head(dispersions(dds))
```
```r
head(mcols(dds)$dispersion)
```
```r
sizeFactors(dds)
```
```r
head(coef(dds))
```
```r
attr(dds, "betaPriorVar")
```
```r
priorInfo(resLFC)
```
```r
priorInfo(resNorm)
```
```r
priorInfo(resAsh)
```
```r
dispersionFunction(dds)
```
```r
attr(dispersionFunction(dds), "dispPriorVar")
```
```r
metadata(dds)[["version"]]
```

## Known Issues
The original dataset included over 500 samples. I narrowed it down to 300 and then further narrowed it down to 40 samples for the analysis. Since there were no filters avaliable, I chose 20 samples for each group at random to speed up the process of the analysis. 

## Conclusions 




