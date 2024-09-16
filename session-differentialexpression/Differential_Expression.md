Differential Expression
================

Created by: Ahmed Mahfouz Edited by: Mohammed Charrout, Lieke Michielsen

# Overview

In this tutorial we will explore different methods to perform
differential expression analysis on scRNA-seq data. The exercises are
based on Seurat’s differential expression testing
[vignette](https://satijalab.org/seurat/articles/de_vignette.html).

Load required packages:

``` r
require(Seurat)
require(scran)
require(scater)
require(pheatmap)
```

We will continue with our PBMC dataset.

``` r
pbmc <- readRDS(file = "../session-clustering/pbmc3k.rds")
```

## Differential expression testing in Seurat

In Seurat, differential expression analysis can be performed using the
`FindMarkers` function. As a default, Seurat performs differential
expression based on the non-parameteric Wilcoxon rank sum test.
Differential expression is performed between groups of cells. To test
for differential expression between two specific groups of cells,
specify the `ident.1` and `ident.2` parameters. The function will
automatically retrieve the cluster identities from the Seurat object
using the `Idents()` function.

Before applying the function, we first have to change the identities to
the original `celltype` column, as we have changed them in the
clustering lab.

``` r
levels(pbmc)
```

    ## [1] "NK and T cells" "Monocytes"      "B cells"

``` r
Idents(pbmc) <- 'celltype'
levels(pbmc)
```

    ## [1] "Monocyte"       "B cell"         "CD8 T cell"     "CD4 T cell"    
    ## [5] "NK cell"        "Dendritic cell"

``` r
# Find differentially expressed features between CD8 and CD4 T-cells
tcell.de.markers <- FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell")
# View results
head(tcell.de.markers)
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 7.901233e-61   7.965983 0.741 0.027 2.649915e-56
    ## CD8A 2.501225e-55   5.527868 0.695 0.020 8.388608e-51
    ## CTSW 1.230991e-54   3.030683 0.925 0.314 4.128497e-50
    ## CCL5 2.573606e-50   4.168283 0.764 0.109 8.631360e-46
    ## GZMA 2.840748e-45   4.292352 0.707 0.109 9.527300e-41
    ## CST7 4.187390e-41   5.035465 0.609 0.058 1.404367e-36

The results data frame has the following columns:

- `p_val`: Unadjusted p-value  
- `avg_log2FC`: log fold-change of the average expression between the
  two groups. Positive values indicate that the feature is more highly
  expressed in the first group.  
- `pct.1`: The percentage of cells where the feature is detected in the
  first group.  
- `pct.2`: The percentage of cells where the feature is detected in the
  second group.  
- `p_val_adj`: Adjusted p-value, based on Bonferroni correction using
  all features in the dataset.

If the ident.2 parameter is omitted or set to NULL, `FindMarkers` will
test for differentially expressed features between the group specified
by ident.1 and all other cells.

``` r
# Find differentially expressed features between CD8 T-cells and all other cells, only
# search for positive markers
tcell.de.markers <- FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = NULL, only.pos = TRUE)
# view results
head(tcell.de.markers)
```

    ##               p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## CD8A  1.068397e-126   5.226774 0.695 0.022 3.583191e-122
    ## CD8B   4.247617e-98   6.303952 0.523 0.011  1.424566e-93
    ## CCL5   6.657331e-89   3.177861 0.764 0.101  2.232736e-84
    ## GZMK   9.705468e-88   4.543261 0.615 0.050  3.255020e-83
    ## CTSW   6.102467e-87   2.020273 0.925 0.184  2.046645e-82
    ## TRGC2  3.555947e-72   4.712553 0.546 0.053  1.192593e-67

To increase the speed of marker discovery, particularly for large
datasets, Seurat allows for pre-filtering of features or cells. For
example, features that are very infrequently detected in either group of
cells, or features that are expressed at similar average levels, are
unlikely to be differentially expressed. Example use cases of the
`min.pct`, `logfc.threshold`, `min.diff.pct`, and `max.cells.per.ident`
parameters are demonstrated below.

``` r
# Pre-filter features that are detected at <50% frequency in either CD8 T-cells or CD4 T-cells. 
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", min.pct = 0.5))
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 7.901233e-61   7.965983 0.741 0.027 2.649915e-56
    ## CD8A 2.501225e-55   5.527868 0.695 0.020 8.388608e-51
    ## CTSW 1.230991e-54   3.030683 0.925 0.314 4.128497e-50
    ## CCL5 2.573606e-50   4.168283 0.764 0.109 8.631360e-46
    ## GZMA 2.840748e-45   4.292352 0.707 0.109 9.527300e-41
    ## CST7 4.187390e-41   5.035465 0.609 0.058 1.404367e-36

``` r
# Pre-filter features that have less than a two-fold change between the average expression of
# CD8 T-cells vs CD4 T-cells.
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", logfc.threshold = log(2)))
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 7.901233e-61   7.965983 0.741 0.027 2.649915e-56
    ## CD8A 2.501225e-55   5.527868 0.695 0.020 8.388608e-51
    ## CTSW 1.230991e-54   3.030683 0.925 0.314 4.128497e-50
    ## CCL5 2.573606e-50   4.168283 0.764 0.109 8.631360e-46
    ## GZMA 2.840748e-45   4.292352 0.707 0.109 9.527300e-41
    ## CST7 4.187390e-41   5.035465 0.609 0.058 1.404367e-36

``` r
# Pre-filter features whose detection percentages across the two groups are similar (within 0.25)
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", min.diff.pct = 0.25))
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 7.901233e-61   7.965983 0.741 0.027 2.649915e-56
    ## CD8A 2.501225e-55   5.527868 0.695 0.020 8.388608e-51
    ## CTSW 1.230991e-54   3.030683 0.925 0.314 4.128497e-50
    ## CCL5 2.573606e-50   4.168283 0.764 0.109 8.631360e-46
    ## GZMA 2.840748e-45   4.292352 0.707 0.109 9.527300e-41
    ## CST7 4.187390e-41   5.035465 0.609 0.058 1.404367e-36

Finally, you can also identify all cluster markers in one go using
`FindAllMarkers`.

``` r
head(FindAllMarkers(pbmc, logfc.threshold = log(2), min.pct = 0.5, min.diff.pct = 0.25))
```

    ## Calculating cluster Monocyte

    ## Calculating cluster B cell

    ## Calculating cluster CD8 T cell

    ## Calculating cluster CD4 T cell

    ## Calculating cluster NK cell

    ## Calculating cluster Dendritic cell

    ##                  p_val avg_log2FC pct.1 pct.2     p_val_adj  cluster     gene
    ## FCN1     2.665083e-214   8.526285 0.994 0.012 8.938154e-210 Monocyte     FCN1
    ## SERPINA1 3.087393e-212   8.888234 0.982 0.008 1.035450e-207 Monocyte SERPINA1
    ## FGL2     3.210232e-203   6.326080 0.991 0.042 1.076648e-198 Monocyte     FGL2
    ## MNDA     9.642581e-203   6.877179 0.979 0.032 3.233929e-198 Monocyte     MNDA
    ## CSTA     1.058587e-200   6.461783 0.955 0.012 3.550290e-196 Monocyte     CSTA
    ## VCAN     2.667989e-199   9.188847 0.943 0.012 8.947902e-195 Monocyte     VCAN

### Alternative DE tests in Seurat

The following differential expression tests are currently supported by
Seurat:

- `wilcox`: Wilcoxon rank sum test (default)  
- `bimod`: Likelihood-ratio test for single cell feature expression,
  (McDavid et al., Bioinformatics, 2013)  
- `roc`: Standard AUC classifier  
- `t`: Student’s t-test  
- `poisson`: Likelihood ratio test assuming an underlying negative
  binomial distribution.  
- `negbinom`: Likelihood ratio test assuming an underlying negative
  binomial distribution.  
- `LR`: Uses a logistic regression framework to determine differentially
  expressed genes. Constructs a logistic regression model predicting
  group membership based on each feature individually and compares this
  to a null model with a likelihood ratio test.  
- `MAST`: GLM-framework that treates cellular detection rate as a
  covariate (Finak et al, Genome Biology, 2015)  
- `DESeq2`: DE based on a model using the negative binomial distribution
  (Love et al, Genome Biology, 2014)

For MAST and DESeq2 please ensure that these packages are installed
separately in order to use them as part of Seurat. Once installed, the
`test.use` parameter can be used to specify which DE test to use.

``` r
# Test for DE features using the MAST package
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", test.use = "wilcox"))
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 7.901233e-61   7.965983 0.741 0.027 2.649915e-56
    ## CD8A 2.501225e-55   5.527868 0.695 0.020 8.388608e-51
    ## CTSW 1.230991e-54   3.030683 0.925 0.314 4.128497e-50
    ## CCL5 2.573606e-50   4.168283 0.764 0.109 8.631360e-46
    ## GZMA 2.840748e-45   4.292352 0.707 0.109 9.527300e-41
    ## CST7 4.187390e-41   5.035465 0.609 0.058 1.404367e-36

``` r
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", test.use = "bimod"))
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 1.000273e-70   7.965983 0.741 0.027 3.354715e-66
    ## CTSW 5.034059e-67   3.030683 0.925 0.314 1.688323e-62
    ## CD8A 1.780653e-59   5.527868 0.695 0.020 5.971953e-55
    ## CCL5 1.357176e-53   4.168283 0.764 0.109 4.551698e-49
    ## PRF1 5.107729e-51   4.395039 0.609 0.085 1.713030e-46
    ## GZMA 1.399478e-50   4.292352 0.707 0.109 4.693569e-46

``` r
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", test.use = "t"))
```

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## CTSW 2.750562e-60   3.030683 0.925 0.314 9.224834e-56
    ## CCL5 8.151017e-44   4.168283 0.764 0.109 2.733688e-39
    ## NKG7 1.645896e-43   7.965983 0.741 0.027 5.520007e-39
    ## CD8A 1.194938e-38   5.527868 0.695 0.020 4.007582e-34
    ## GZMA 2.669118e-37   4.292352 0.707 0.109 8.951686e-33
    ## HCST 1.592907e-31   1.287539 0.971 0.816 5.342291e-27

``` r
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", test.use = "poisson"))
```

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 13. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 25. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 60. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 69. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 72. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 184. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 333. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 353. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 378. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 395. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 411. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 494. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 504. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 510. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 542. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 543. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 598. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 615. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 632. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 649. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 656. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 658. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 696. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 720. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 750. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 767. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 775. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 788. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 791. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 817. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 819. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 821. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 840. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 845. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 871. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 943. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 947. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 964. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 984. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 992. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1002. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1041. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1073. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1116. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1127. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1129. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1174. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1178. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1201. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1255. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1258. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1325. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1326. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1356. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1376. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1395. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1408. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1444. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1514. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1520. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1563. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1569. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1588. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1590. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1607. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1633. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1670. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1688. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1729. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1736. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1740. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1848. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1903. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1921. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1928. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1936. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1938. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1943. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 1960. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2007. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2012. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2071. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2090. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2146. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2177. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2197. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2222. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2275. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2283. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2284. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2288. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2352. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2356. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2359. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2390. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2405. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2419. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2499. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2509. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2521. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2541. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2557. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2602. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2617. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2645. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2649. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2653. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2703. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2716. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2735. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2754. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2757. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2769. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2862. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2928. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2931. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 2991. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3021. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3032. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3060. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3062. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3074. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3081. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3086. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3160. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3191. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3198. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3258. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3282. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3290. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3292. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3297. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3359. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3368. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3415. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3479. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3504. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3524. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3534. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3553. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3587. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3589. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3651. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3652. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3653. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3669. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3670. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3680. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3707. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3733. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3742. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3809. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3811. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3817. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3846. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3860. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3868. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3875. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3877. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3909. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3911. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 3974. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4050. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4060. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4072. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4085. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4094. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4109. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4220. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4230. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4232. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4237. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4243. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4301. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4348. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4360. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4385. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4401. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4468. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4485. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4594. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4661. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4721. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4738. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4742. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4755. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4802. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4835. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4865. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4879. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4886. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4889. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4892. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4893. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4898. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4938. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 4953. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5042. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5043. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5044. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5070. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5086. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5129. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5135. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5192. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5216. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5241. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5242. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5264. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5279. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5288. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5365. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5427. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5547. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5548. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5553. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5561. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5577. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5585. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5594. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5636. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5652. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5658. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5667. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5669. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5699. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5742. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5823. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5826. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5855. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5857. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5885. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5902. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5936. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5944. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5967. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 5983. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6021. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6051. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6151. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6173. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6174. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6205. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6218. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6265. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6280. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6314. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6340. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6358. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6387. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6392. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6409. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6415. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6418. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6424. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6430. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6520. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6531. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6542. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6581. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6587. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6599. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6635. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6640. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6647. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6677. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6721. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6743. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6812. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6832. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6838. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6903. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 6910. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7002. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7077. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7094. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7142. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7180. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7194. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7266. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7284. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7345. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7370. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7395. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7398. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7401. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7403. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7425. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7473. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7474. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7484. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7488. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7510. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7517. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7550. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7553. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7556. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7566. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7645. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7691. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7710. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7735. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7760. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7766. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7793. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7794. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7806. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7829. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7865. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7868. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 7903. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8012. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8027. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8032. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8036. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8047. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8051. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8060. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8095. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8097. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8116. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8135. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8139. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8140. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8144. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8194. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8199. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8209. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8211. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8220. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8263. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8274. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8296. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8310. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8350. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8360. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8371. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8376. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8425. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8439. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8459. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8501. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8512. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8514. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8589. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8602. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8613. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8657. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8685. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8698. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8712. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8726. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8753. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8778. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8788. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8805. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8828. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8848. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8860. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8897. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8901. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8962. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8976. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 8983. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9002. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9025. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9058. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9068. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9078. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9102. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9126. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9133. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9143. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9156. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9159. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9161. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9174. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9176. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9182. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9190. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9204. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9215. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9281. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9291. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9293. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9347. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9353. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9366. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9384. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9405. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9435. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9456. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9489. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9491. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9503. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9509. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9524. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9558. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9605. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9606. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9616. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9628. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9631. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9659. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9703. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9705. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9707. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9712. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9714. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9727. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9759. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9760. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9837. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9860. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9861. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9873. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9927. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 9938. Fewer than 3 cells in both
    ## clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10005. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10018. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10032. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10100. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10105. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10106. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10134. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10205. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10242. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10267. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10269. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10297. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10298. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10329. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10335. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10338. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10352. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10371. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10388. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10456. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10462. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10476. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10478. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10506. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10512. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10514. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10539. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10541. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10571. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10575. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10587. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10596. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10616. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10618. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10639. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10665. Fewer than 3 cells in
    ## both clusters.

    ## Warning in FUN(X[[i]], ...): Skipping gene --- 10692. Fewer than 3 cells in
    ## both clusters.

    ##                p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## CCL5   1.705622e-222  4.1682833 0.764 0.109 5.720313e-218
    ## TPT1   1.283944e-208 -0.4763124 1.000 1.000 4.306092e-204
    ## RPS27  5.031472e-155 -0.3240157 1.000 1.000 1.687455e-150
    ## KLRB1  6.476062e-142  2.1906689 0.511 0.270 2.171942e-137
    ## EEF1A1 3.732156e-129 -0.2574449 1.000 1.000 1.251690e-124
    ## MT-CO1 5.296595e-125  0.4367286 1.000 1.000 1.776372e-120

``` r
head(FindMarkers(pbmc, ident.1 = "CD8 T cell", ident.2 = "CD4 T cell", test.use = "MAST"))
```

    ## 
    ## Done!

    ## Combining coefficients and standard errors

    ## Calculating log-fold changes

    ## Calculating likelihood ratio tests

    ## Refitting on reduced model...

    ## 
    ## Done!

    ##             p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## NKG7 1.506208e-70   7.965983 0.741 0.027 5.051522e-66
    ## CTSW 7.684912e-63   3.030683 0.925 0.314 2.577366e-58
    ## CD8A 1.525308e-60   5.527868 0.695 0.020 5.115579e-56
    ## CCL5 1.135653e-54   4.168283 0.764 0.109 3.808754e-50
    ## GZMA 1.462596e-50   4.292352 0.707 0.109 4.905256e-46
    ## CST7 8.473749e-44   5.035465 0.609 0.058 2.841926e-39

## Differential expression analysis using scran

The `findMarkers()` function in scran uses a different approach to
identify marker genes compared to Seurat. While in Seurat the default is
to perform one vs all comparisons, `findMarkers()` in scran performs
pairwise comparisons between clusters for each gene. The default test in
`findMarkers()` is the Welch t-test.

Scran intentionally uses pairwise comparisons between clusters rather
than comparing each cluster to the average of all other cells. The
latter approach is sensitive to the population composition, potentially
resulting in substantially different sets of markers when cell type
abundances change in different contexts. In the worst case, the presence
of a single dominant subpopulation will drive the selection of top
markers for every other cluster, pushing out useful genes that can
resolve the various minor subpopulations.

First, let’s convert our Seurat object to a SingleCellExperiment object.

``` r
pbmc.sce <- as.SingleCellExperiment(pbmc)
```

`findMarkers()` returns a list of data frames containing ranked
candidate markers for each cluster.

``` r
markers.pbmc <- findMarkers(pbmc.sce, groups=pbmc.sce$ident)
```

You can then choose one data frame (in this example, corresponding to
CD8 T-cells). This data frame contains log2-fold changes of expression
in the chosen cluster over each other cluster as well as several
statistics obtained by combining p-values across the pairwise
comparisons involving the cluster of interest.

``` r
chosen <- "CD8 T cell"
interesting <- markers.pbmc[[chosen]]
interesting[1:10,1:4]
```

    ## DataFrame with 10 rows and 4 columns
    ##                Top      p.value          FDR summary.logFC
    ##          <integer>    <numeric>    <numeric>     <numeric>
    ## HLA-DRB1         1 6.03248e-176 1.12398e-172      -2.95127
    ## CTSW             1  2.43221e-72  4.77027e-70       1.65306
    ## CD3D             1  5.83999e-96  1.99859e-93       1.83222
    ## IL32             1  1.13157e-97  4.12505e-95       2.47244
    ## CST3             1  0.00000e+00  0.00000e+00      -2.86031
    ## LCK              2  2.74067e-75  5.81750e-73       1.44093
    ## FCER1G           2 7.94116e-310 1.33165e-305      -2.58399
    ## CD3G             2  4.02478e-61  5.79326e-59       1.29761
    ## CCL5             2  1.35929e-47  1.36490e-45       2.03777
    ## CD79A            2 6.39837e-171 1.12941e-167      -2.84229

The `summary.logFC` field provides a summary of the direction and effect
size for each gene. `logFC` is defined here as the log-fold change from
the comparison with the lowest p-value. The `p.value` field contains the
combined p-value that is obtained by applying Simes’ method to the
pairwise p-values for each gene. Of particular interest is the `Top`
field. The set of genes with `Top` $≤ X$ is the union of the top $X$
genes (ranked by p-value) from each pairwise comparison involving the
cluster of interest.

Let’s plot a heatmap of the top 5 genes for CD8 T-cells.

``` r
best.set <- interesting[interesting$Top <= 5,]
logFCs <- getMarkerEffects(best.set)
pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))
```

![](Differential_Expression_files/figure-gfm/scran_e-1.png)<!-- -->

### Wilcoxon vs t-test

Also in scran, you can use different DE tests. Beside the default Welch
t-test, you can also use a Wilcoxon rank-sum test or a binomial test.

``` r
markers.pbmc.wrs <- findMarkers(pbmc.sce, groups=pbmc.sce$ident, test="wilcox")
interesting.wrs <- markers.pbmc.wrs[[chosen]]
interesting.wrs[1:10,1:4]
```

    ## DataFrame with 10 rows and 4 columns
    ##              Top      p.value         FDR summary.AUC
    ##        <integer>    <numeric>   <numeric>   <numeric>
    ## FCER1G         1  1.16092e-78 3.89351e-75 1.72586e-05
    ## MARCH1         1  2.87753e-64 9.94913e-62 5.70571e-02
    ## CD3E           1 3.22427e-103 1.08136e-98 9.97101e-01
    ## CD79A          1  8.03182e-66 3.28502e-63 0.00000e+00
    ## NKG7           1  3.95062e-60 1.01920e-57 8.66100e-01
    ## CD8A           2  1.13117e-65 4.57073e-63 8.47701e-01
    ## KLRF1          2  2.29872e-32 1.59947e-30 4.91071e-02
    ## RNASE6         2  9.78777e-40 9.71190e-38 2.94118e-02
    ## IL32           2  1.48340e-99 1.72695e-95 9.88092e-01
    ## TGFBI          3  1.60505e-62 5.03085e-60 6.32270e-02

One advantage of the Wilcoxon rank-sum test over the Welch t-test is
that it is symmetric with respect to differences in the size of the
groups being compared. In other words, it is less affected by the number
of cells in each group. On the other hand, the t-test will favor genes
where the larger group has the higher relative variance as this
increases the estimated degrees of freedom and decreases the resulting
p-value.

To illustrate this we will use an example from [“Orchestrating
Single-Cell Analysis with
Bioconductor”](https://osca.bioconductor.org/marker-detection.html). In
this example, we will compare alpha and gamma cells in the human
pancreas data set from Lawlor et al. (2017)

``` r
sce.lawlor <- readRDS(file = "sce_lawlor.rds")
marker.lawlor.t <- findMarkers(sce.lawlor, groups=sce.lawlor$`cell type`, 
                               direction="up", restrict=c("Alpha", "Gamma/PP"))
marker.lawlor.w <- findMarkers(sce.lawlor, groups=sce.lawlor$`cell type`, 
                               direction="up", restrict=c("Alpha", "Gamma/PP"), test.type="wilcox")
# Upregulated in alpha:
marker.alpha.t <- marker.lawlor.t$Alpha
marker.alpha.w <- marker.lawlor.w$Alpha
chosen.alpha.t <- rownames(marker.alpha.t)[1:5]
chosen.alpha.w <- rownames(marker.alpha.w)[1:5]
u.alpha.t <- setdiff(chosen.alpha.t, chosen.alpha.w)
u.alpha.w <- setdiff(chosen.alpha.w, chosen.alpha.t)
# Upregulated in gamma:
marker.gamma.t <- marker.lawlor.t$`Gamma/PP`
marker.gamma.w <- marker.lawlor.w$`Gamma/PP`
chosen.gamma.t <- rownames(marker.gamma.t)[1:5]
chosen.gamma.w <- rownames(marker.gamma.w)[1:5]
u.gamma.t <- setdiff(chosen.gamma.t, chosen.gamma.w)
u.gamma.w <- setdiff(chosen.gamma.w, chosen.gamma.t)
# Examining all uniquely detected markers in each direction.
subset <- sce.lawlor[,sce.lawlor$`cell type` %in% c("Alpha", "Gamma/PP")]
gridExtra::grid.arrange(
  plotExpression(subset, x="cell type", features=u.alpha.t, ncol=2) +
    ggtitle("Upregulated in alpha, t-test-only"),
  plotExpression(subset, x="cell type", features=u.alpha.w, ncol=2) +
    ggtitle("Upregulated in alpha, WMW-test-only"),
  plotExpression(subset, x="cell type", features=u.gamma.t, ncol=2) +
    ggtitle("Upregulated in gamma, t-test-only"),
  plotExpression(subset, x="cell type", features=u.gamma.w, ncol=2) +
    ggtitle("Upregulated in gamma, WMW-test-only"),
  ncol=2
)
```

![](Differential_Expression_files/figure-gfm/scran_g-1.png)<!-- -->

Can you observe the effects of the tests in the resulting genes?

## DE testing for integrated data

Nowadays, it is common to work with multiple scRNA-seq datasets. As we
have seen in the integration practical, several strategies exist to
integrate multiple datasets and perform cell-based analysis
(e.g. clustering or trajectory inference) using the integrated data
(i.e. data corrected for batch effects).

But what about gene-based analysis? Can we use the integrated
(i.e. corrected) data for differential expression analysis? In general,
this is not recommended. The reason is that arbitrary correction
algorithms do not preserve the magnitude or the direction of differences
in per-gene expression when attempting to align multiple batches.
Further, the correction can introduce artificial agreement across
batches. For a good discussion of these implications, check [chapter
13.7](https://osca.bioconductor.org/integrating-datasets.html#using-corrected-values)
in [“Orchestrating Single-Cell Analysis with
Bioconductor”](https://osca.bioconductor.org/marker-detection.html).

For these reasons, it is preferred to perform DE testing using the
uncorrected expression data. There are two strategies to handle the
batch effects in this case.

### Identifying conserved markers (Seurat)

To identify canonical cell type marker genes that are conserved across
batches, Seurat provides the `FindConservedMarkers` function. This
function performs differential gene expression testing for each
dataset/group/batch separately and combines the p-values using
meta-analysis methods from the `MetaDE` R package.

To illustrate this, we will again load in the other PBMC datasets and
merge them with our PBMC dataset. The cell type labels from the
integration practical will be used again to identify the cell
populations.

``` r
# Load and process v2.1k dataset
v2.1k <- Read10X_h5("../session-qc-normalization/pbmc_1k_v2_filtered_feature_bc_matrix.h5")
meta.data <- read.table('../session-integration/celltypes_1k_v2.tsv', sep = '\t', header = TRUE, row.names = 1)
v2.1k <- CreateSeuratObject(v2.1k, meta.data = meta.data, project = 'v2.1k')
v2.1k <- NormalizeData(v2.1k)
```

    ## Normalizing layer: counts

``` r
# Load and process p3.1k dataset
p3.1k <- Read10X_h5("../session-qc-normalization/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5")
```

    ## Genome matrix has multiple modalities, returning a list of matrices for this genome

``` r
meta.data <- read.table('../session-integration/celltypes_1k_protein.tsv', sep = '\t', header = TRUE, row.names = 1)
p3.1k <- CreateSeuratObject(p3.1k$`Gene Expression`, meta.data = meta.data, project = 'p3.1k')
p3.1k <- NormalizeData(p3.1k)
```

    ## Normalizing layer: counts

``` r
# Merge
pbmc.all <- merge(x = pbmc, y = c(v2.1k, p3.1k), add.cell.ids = c("v3.1k", "v2.1k", "p3.1k"))
pbmc.all <- AddMetaData(pbmc.all, sub("_.*", "", colnames(pbmc.all)), "batch")
Idents(pbmc.all) <- 'celltype'
table(pbmc.all[[c("batch", "celltype")]])
```

    ##        celltype
    ## batch   B cell CD4 T cell CD8 T cell Dendritic cell Monocyte NK cell Other
    ##   p3.1k     80        145         97             11      245      31   104
    ##   v2.1k    142        300        137             16      322      65    14
    ##   v3.1k    181        293        174             17      333      56     0

Now find the conserved markers for CD8 T-cells by batch.

``` r
pbmc.all <- JoinLayers(pbmc.all)
cd8.markers <- FindConservedMarkers(pbmc.all, ident.1 = "CD8 T cell", grouping.var = "batch", 
                                    logfc.threshold = 0.1, min.pct = 0.25, min.diff.pct = 0.25)
```

    ## Testing group p3.1k: (CD8 T cell) vs (B cell, CD4 T cell, Monocyte, Dendritic cell, Other, NK cell)

    ## Testing group v2.1k: (CD8 T cell) vs (Monocyte, CD4 T cell, NK cell, Other, B cell, Dendritic cell)

    ## Testing group v3.1k: (CD8 T cell) vs (Monocyte, B cell, CD4 T cell, NK cell, Dendritic cell)

``` r
head(cd8.markers, n=10)
```

    ##        p3.1k_p_val p3.1k_avg_log2FC p3.1k_pct.1 p3.1k_pct.2 p3.1k_p_val_adj
    ## CD8A  5.606455e-46         3.780649       0.474       0.029    1.880293e-41
    ## CD8B  9.689188e-33         4.355596       0.320       0.016    3.249560e-28
    ## CCL5  2.644436e-44         3.469557       0.670       0.112    8.868911e-40
    ## GZMK  1.630749e-43         3.630495       0.485       0.037    5.469207e-39
    ## CTSW  7.522387e-28         1.855259       0.619       0.138    2.522858e-23
    ## TRGC2 1.445072e-32         3.968144       0.412       0.042    4.846484e-28
    ## KLRG1 7.329561e-27         2.558725       0.526       0.104    2.458188e-22
    ## GZMA  5.359275e-31         2.077133       0.515       0.073    1.797394e-26
    ## NKG7  6.793357e-26         1.547488       0.557       0.119    2.278356e-21
    ## IL32  8.959317e-27         1.906753       0.773       0.255    3.004776e-22
    ##        v2.1k_p_val v2.1k_avg_log2FC v2.1k_pct.1 v2.1k_pct.2 v2.1k_p_val_adj
    ## CD8A  5.924812e-84         5.226761       0.533       0.021    1.987063e-79
    ## CD8B  5.977072e-90         6.127611       0.526       0.014    2.004590e-85
    ## CCL5  1.482409e-67         2.623989       0.774       0.135    4.971704e-63
    ## GZMK  7.670258e-71         4.625634       0.496       0.027    2.572451e-66
    ## CTSW  3.568178e-54         1.855343       0.788       0.169    1.196695e-49
    ## TRGC2 7.867071e-55         5.059770       0.401       0.024    2.638458e-50
    ## KLRG1 1.306959e-39         3.023087       0.445       0.063    4.383279e-35
    ## GZMA  1.139184e-39         1.805020       0.620       0.123    3.820595e-35
    ## NKG7  3.273629e-52         1.284913       0.745       0.144    1.097910e-47
    ## IL32  8.653273e-56         2.041696       0.978       0.336    2.902135e-51
    ##         v3.1k_p_val v3.1k_avg_log2FC v3.1k_pct.1 v3.1k_pct.2 v3.1k_p_val_adj
    ## CD8A  1.068397e-126         5.226774       0.695       0.022   3.583191e-122
    ## CD8B   4.247617e-98         6.303952       0.523       0.011    1.424566e-93
    ## CCL5   6.657331e-89         3.177861       0.764       0.101    2.232736e-84
    ## GZMK   9.705468e-88         4.543261       0.615       0.050    3.255020e-83
    ## CTSW   6.102467e-87         2.020273       0.925       0.184    2.046645e-82
    ## TRGC2  3.555947e-72         4.712553       0.546       0.053    1.192593e-67
    ## KLRG1  5.363311e-72         3.292885       0.690       0.116    1.798747e-67
    ## GZMA   1.635706e-69         2.297356       0.707       0.105    5.485831e-65
    ## NKG7   2.496780e-69         1.851690       0.741       0.127    8.373701e-65
    ## IL32   9.307544e-66         2.093265       0.977       0.334    3.121564e-61
    ##           max_pval minimump_p_val
    ## CD8A  5.606455e-46  3.205192e-126
    ## CD8B  9.689188e-33   1.274285e-97
    ## CCL5  2.644436e-44   1.997199e-88
    ## GZMK  1.630749e-43   2.911640e-87
    ## CTSW  7.522387e-28   1.830740e-86
    ## TRGC2 1.445072e-32   1.066784e-71
    ## KLRG1 7.329561e-27   1.608993e-71
    ## GZMA  5.359275e-31   4.907119e-69
    ## NKG7  6.793357e-26   7.490340e-69
    ## IL32  8.959317e-27   2.792263e-65

### Identifying consistent markers (scran)

Alternatively, we can perform DE tests on the uncorrected data with
blocking on the batch. The rational here is that we expect true DE
between clusters to be present in a within-batch comparison where batch
effects are absent. This strategy penalizes genes that exhibit
inconsistent DE across batches.

``` r
pbmc.all.sce <- as.SingleCellExperiment(pbmc.all)
m.out <- findMarkers(pbmc.all.sce, group = pbmc.all.sce$celltype, block = pbmc.all.sce$batch,
    direction="up", lfc=1)
m.out$`CD8 T cell`[1:20, c("Top", "p.value", "FDR")]
```

    ## DataFrame with 20 rows and 3 columns
    ##             Top      p.value          FDR
    ##       <integer>    <numeric>    <numeric>
    ## RPS12         1  1.40573e-11  2.77325e-08
    ## CD3D          1  1.37286e-72  1.53477e-68
    ## IL32          1 5.88491e-125 1.97368e-120
    ## CCL5          1  1.28965e-44  7.20871e-41
    ## CD3E          2  5.53481e-79  9.28131e-75
    ## ...         ...          ...          ...
    ## RPS8          5  3.74089e-09  5.97438e-06
    ## GZMK          5  7.82838e-04  3.05289e-01
    ## TRBC2         5  5.66708e-26  1.68322e-22
    ## RPS27         6  3.98428e-09  6.07386e-06
    ## CD3G          6  6.02648e-05  2.80717e-02

``` r
plotExpression(pbmc.all.sce, x=I(factor(pbmc.all.sce$celltype)), 
    features="CD3D", colour_by="batch") + facet_wrap(~colour_by) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](Differential_Expression_files/figure-gfm/integration_c-1.png)<!-- -->

### Session info

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/Madrid
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] pheatmap_1.0.12             scater_1.32.1              
    ##  [3] ggplot2_3.5.1               scran_1.32.0               
    ##  [5] scuttle_1.14.0              SingleCellExperiment_1.26.0
    ##  [7] SummarizedExperiment_1.34.0 Biobase_2.64.0             
    ##  [9] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
    ## [11] IRanges_2.38.1              S4Vectors_0.42.1           
    ## [13] BiocGenerics_0.50.0         MatrixGenerics_1.16.0      
    ## [15] matrixStats_1.4.1           Seurat_5.1.0               
    ## [17] SeuratObject_5.0.2          sp_2.1-4                   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.22          splines_4.4.1            
    ##   [3] later_1.3.2               tibble_3.2.1             
    ##   [5] polyclip_1.10-7           fastDummies_1.7.4        
    ##   [7] lifecycle_1.0.4           Rdpack_2.6.1             
    ##   [9] edgeR_4.2.1               hdf5r_1.3.11             
    ##  [11] globals_0.16.3            lattice_0.22-6           
    ##  [13] MASS_7.3-61               MAST_1.30.0              
    ##  [15] magrittr_2.0.3            limma_3.60.4             
    ##  [17] plotly_4.10.4             rmarkdown_2.28           
    ##  [19] plotrix_3.8-4             yaml_2.3.10              
    ##  [21] qqconf_1.3.2              metapod_1.12.0           
    ##  [23] httpuv_1.6.15             sn_2.1.1                 
    ##  [25] sctransform_0.4.1         spam_2.10-0              
    ##  [27] spatstat.sparse_3.1-0     reticulate_1.39.0        
    ##  [29] cowplot_1.1.3             pbapply_1.7-2            
    ##  [31] RColorBrewer_1.1-3        multcomp_1.4-26          
    ##  [33] abind_1.4-8               zlibbioc_1.50.0          
    ##  [35] Rtsne_0.17                purrr_1.0.2              
    ##  [37] presto_1.0.0              TH.data_1.1-2            
    ##  [39] sandwich_3.1-1            GenomeInfoDbData_1.2.12  
    ##  [41] ggrepel_0.9.6             irlba_2.3.5.1            
    ##  [43] listenv_0.9.1             spatstat.utils_3.1-0     
    ##  [45] TFisher_0.2.0             goftest_1.2-3            
    ##  [47] RSpectra_0.16-2           spatstat.random_3.3-1    
    ##  [49] dqrng_0.4.1               fitdistrplus_1.2-1       
    ##  [51] parallelly_1.38.0         DelayedMatrixStats_1.26.0
    ##  [53] leiden_0.4.3.1            codetools_0.2-20         
    ##  [55] DelayedArray_0.30.1       tidyselect_1.2.1         
    ##  [57] farver_2.1.2              UCSC.utils_1.0.0         
    ##  [59] ScaledMatrix_1.12.0       viridis_0.6.5            
    ##  [61] spatstat.explore_3.3-2    mathjaxr_1.6-0           
    ##  [63] jsonlite_1.8.8            multtest_2.60.0          
    ##  [65] BiocNeighbors_1.22.0      progressr_0.14.0         
    ##  [67] ggridges_0.5.6            survival_3.7-0           
    ##  [69] progress_1.2.3            tools_4.4.1              
    ##  [71] ica_1.0-3                 Rcpp_1.0.13              
    ##  [73] glue_1.7.0                mnormt_2.1.1             
    ##  [75] gridExtra_2.3             SparseArray_1.4.8        
    ##  [77] metap_1.11                xfun_0.47                
    ##  [79] dplyr_1.1.4               numDeriv_2016.8-1.1      
    ##  [81] withr_3.0.1               fastmap_1.2.0            
    ##  [83] bluster_1.14.0            fansi_1.0.6              
    ##  [85] digest_0.6.37             rsvd_1.0.5               
    ##  [87] R6_2.5.1                  mime_0.12                
    ##  [89] colorspace_2.1-1          scattermore_1.2          
    ##  [91] tensor_1.5                spatstat.data_3.1-2      
    ##  [93] utf8_1.2.4                tidyr_1.3.1              
    ##  [95] generics_0.1.3            data.table_1.16.0        
    ##  [97] prettyunits_1.2.0         httr_1.4.7               
    ##  [99] htmlwidgets_1.6.4         S4Arrays_1.4.1           
    ## [101] uwot_0.2.2                pkgconfig_2.0.3          
    ## [103] gtable_0.3.5              lmtest_0.9-40            
    ## [105] XVector_0.44.0            htmltools_0.5.8.1        
    ## [107] dotCall64_1.1-1           scales_1.3.0             
    ## [109] png_0.1-8                 spatstat.univar_3.0-1    
    ## [111] knitr_1.48                reshape2_1.4.4           
    ## [113] nlme_3.1-166              zoo_1.8-12               
    ## [115] stringr_1.5.1             KernSmooth_2.23-24       
    ## [117] parallel_4.4.1            miniUI_0.1.1.1           
    ## [119] vipor_0.4.7               pillar_1.9.0             
    ## [121] grid_4.4.1                vctrs_0.6.5              
    ## [123] RANN_2.6.2                promises_1.3.0           
    ## [125] BiocSingular_1.20.0       beachmat_2.20.0          
    ## [127] xtable_1.8-4              cluster_2.1.6            
    ## [129] beeswarm_0.4.0            evaluate_0.24.0          
    ## [131] mvtnorm_1.3-1             cli_3.6.3                
    ## [133] locfit_1.5-9.10           compiler_4.4.1           
    ## [135] rlang_1.1.4               crayon_1.5.3             
    ## [137] mutoss_0.1-13             future.apply_1.11.2      
    ## [139] labeling_0.4.3            plyr_1.8.9               
    ## [141] ggbeeswarm_0.7.2          stringi_1.8.4            
    ## [143] viridisLite_0.4.2         deldir_2.0-4             
    ## [145] BiocParallel_1.38.0       munsell_0.5.1            
    ## [147] lazyeval_0.2.2            spatstat.geom_3.3-2      
    ## [149] Matrix_1.7-0              RcppHNSW_0.6.0           
    ## [151] hms_1.1.3                 patchwork_1.2.0          
    ## [153] bit64_4.0.5               sparseMatrixStats_1.16.0 
    ## [155] future_1.34.0             statmod_1.5.0            
    ## [157] shiny_1.9.1               highr_0.11               
    ## [159] rbibutils_2.2.16          ROCR_1.0-11              
    ## [161] igraph_2.0.3              bit_4.0.5
