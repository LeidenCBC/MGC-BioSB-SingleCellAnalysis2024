Trajectory inference
================
Qirong Mao

Acknowledgment: This tutorial of Dimensionality Reduction is based on
the Slingshot tutorial on Bioconductor
([Link](https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html))
and the single cell RNA-seq workshop from Broad Institute
([Link](https://broadinstitute.github.io/2020_scWorkshop/trajectory-analysis.html))

``` r
suppressPackageStartupMessages({
library(SingleCellExperiment)
library(Seurat)
library(ggbeeswarm)
library(ggthemes)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scater)
library(scran)
library(tradeSeq)
library(slingshot)
})
```

## Data loading

We will use a nice SMART-Seq2 single cell RNA-seq data from Single-Cell
RNA-Seq Reveals Dynamic, Random Monoallelic Gene Expression in Mammalian
Cells. Here is one relevant detail from their paper: “To investigate
allele-specific gene expression at single-cell resolution, we isolated
269 individual cells dissociated from in vivo F1 embryos (CAST/EiJ ×
C57BL/6J, hereafter abbreviated as CAST and C57, respectively) from
oocyte to blastocyst stages of mouse preimplantation development (PD)”

``` r
mydir <- "~/MGC-BioSB-SingleCellAnalysis2023-main/session-trajectories/data/"
path.deng <- paste0(mydir, "deng-reads.rds")
deng_SCE <- readRDS(path.deng)
```

``` r
table(deng_SCE$cell_type2)
```

    ## 
    ##     16cell      4cell      8cell early2cell earlyblast  late2cell  lateblast 
    ##         50         14         37          8         43         10         30 
    ##   mid2cell   midblast         zy 
    ##         12         60          4

## Principle Components Analysis

Let us take a first look at the Deng data. One simple approach to
ordering cells in pseudotime is to use PCA. By carrying out PCA and
labeling the cells by the stage at which they were collected, we can see
how well the principal components separate cells along a differentiation
trajectory.

``` r
##
deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy", "early2cell", "mid2cell", "late2cell", 
                                         "4cell", "8cell", "16cell", "earlyblast",
                                         "midblast", "lateblast"))
```

``` r
# Run PCA on Deng data. Use the runPCA function from the SingleCellExperiment package.
deng_SCE <- runPCA(deng_SCE, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results. 
pca <- reducedDim(deng_SCE, "PCA")

# Describe how the PCA is stored in a matrix. Why does it have this structure?
head(pca)
```

    ##                PC1       PC2        PC3         PC4        PC5         PC6
    ## 16cell   -24.79868 -62.20826  -8.035201 -2.07391816 -2.1297390 -14.0930954
    ## 16cell.1 -28.77121 -50.35974 -13.607012  0.08664449 -0.9454185  -3.5987880
    ## 16cell.2 -26.67901 -61.03875  -7.239352 -6.60967794  1.0775002 -11.8876579
    ## 16cell.3 -29.01151 -56.03620  -6.433057  2.85332708  4.2885083   0.1488504
    ## 16cell.4 -26.38026 -58.09265  -4.671850  7.99519397 -9.8077416  -2.0570042
    ## 16cell.5 -24.90566 -60.77897  -5.632497 -3.80156587 -9.8835527 -11.9028394
    ##                 PC7        PC8        PC9      PC10      PC11        PC12
    ## 16cell   -2.4645020  1.6350660  -7.202260  9.862212 10.660702  -0.6401721
    ## 16cell.1 -2.1726663 -3.3481641  -8.967394  6.664942 14.493227 -11.7471565
    ## 16cell.2  7.9007309  0.3368756  -6.032645  5.295515 15.384993  -4.2930696
    ## 16cell.3  4.3727592 -1.1582470  -1.520145 -8.789699 19.386866   0.4999047
    ## 16cell.4  0.6031572 -3.6743278  -5.793753 10.823787  7.613724  -4.7288640
    ## 16cell.5  4.3269009  3.8968881 -11.805221  9.798854 11.016137 -19.1535086
    ##               PC13      PC14      PC15       PC16       PC17       PC18
    ## 16cell   -5.716841  6.544614  6.652210  -3.458346  -4.499013 -11.360753
    ## 16cell.1 13.284708 -4.206404  8.721043  -7.926277  -0.703508  -5.418131
    ## 16cell.2  9.633173  1.672498  9.609001  -9.302794 -10.219743  -5.763834
    ## 16cell.3 14.177687 -8.509097  6.978210  10.771078  -6.188808   6.504081
    ## 16cell.4  3.106382 -4.078414 10.739979 -12.032452  -6.239499   2.331292
    ## 16cell.5  9.544362 -2.255400  8.614958  -2.832196  -1.798584   2.321082
    ##                 PC19      PC20       PC21       PC22       PC23      PC24
    ## 16cell     2.2617345 -2.456274 -11.227414 -1.7122827  -8.418641  4.254968
    ## 16cell.1 -11.8613891  4.069530  -9.320831 -0.5802347 -11.878096 -6.412425
    ## 16cell.2  -3.3460356  4.165813  -2.031473  2.1106373  -1.762218 -1.135134
    ## 16cell.3  -0.6042649  6.008176  -9.982856 -9.4888653   2.822138 12.871921
    ## 16cell.4   3.9402029 -0.298227 -10.773722  0.6374236   4.730329  4.670391
    ## 16cell.5  -2.0280791  5.050525   3.252243  7.1527175  -9.923140 -1.791511
    ##               PC25      PC26       PC27      PC28        PC29       PC30
    ## 16cell   -4.049629  4.133374 -0.6235391 -3.381254 13.94917609  -8.217824
    ## 16cell.1 -8.052083  8.334263 -0.5815629  4.592214  1.32417854   5.266909
    ## 16cell.2 -2.326133  3.775858 -2.3388745 -6.947394  0.08121559  -2.942813
    ## 16cell.3 -5.860750  1.869659  7.0402429  5.092207 -2.53575943 -18.529304
    ## 16cell.4 -4.291113 13.005331  3.2802102  4.606226 -3.52531994  -3.599833
    ## 16cell.5  4.708265  5.717693  1.1023767  9.761377 -4.57312078 -12.138646
    ##                PC31       PC32        PC33      PC34       PC35       PC36
    ## 16cell    -6.897320  -5.675943   8.6076039 -3.713348 -0.9099737  4.7467546
    ## 16cell.1  -4.538307   9.166969  -9.4525575 -8.848231 -2.0782319  7.4318993
    ## 16cell.2   3.082470  -2.207176   0.5365986 -3.895378  7.4493361  0.7465149
    ## 16cell.3   1.680117  -3.839556 -13.3156066 -6.257479 -4.1112596  0.2780589
    ## 16cell.4 -13.314741  -1.453554   0.1334034  2.941487 -0.8162660 -2.9940693
    ## 16cell.5  -4.608498 -12.180530   5.8667454  6.645273  1.0224859  0.8960299
    ##               PC37       PC38       PC39      PC40      PC41       PC42
    ## 16cell   -9.063470 -5.2765051  1.1758453  9.474215  3.559391  4.7781174
    ## 16cell.1 -6.217009  1.0216459  0.5798035 21.705585 -3.570104 -2.3279922
    ## 16cell.2 -6.227582  3.0863112  8.6153521 -1.401230  2.266017 -0.8150665
    ## 16cell.3 -8.411600  3.7169411 -0.7050601  2.959623 -3.123082 -1.0916370
    ## 16cell.4  2.871774 -4.2664023 -7.4894594 -8.207422  4.223035  1.4763578
    ## 16cell.5 10.169730  0.3923632 -9.3346900  8.114487 11.186021  4.5635674
    ##                PC43      PC44       PC45      PC46        PC47      PC48
    ## 16cell   -7.9228097  8.558202  7.0589601  3.058213  -0.5723815  4.675859
    ## 16cell.1  5.6006755 -8.717056 -6.4809594  8.554812 -13.1868751  3.397683
    ## 16cell.2  5.2532880  5.803788  2.7268218  1.241770   7.4824426 -4.088268
    ## 16cell.3 -0.0513553  2.181424  2.4047798  8.691231   8.9700019 -3.713547
    ## 16cell.4  1.5501973 -4.946841  0.5207532 -3.068228  10.7801143  5.167618
    ## 16cell.5 -9.9821175  8.759947 -3.7277578 -9.064882  -1.7524451 -3.306568
    ##               PC49       PC50
    ## 16cell   -2.937621 -3.5882779
    ## 16cell.1 -3.419868  3.4878920
    ## 16cell.2  4.445723  0.2352798
    ## 16cell.3 -5.179552 -9.7180911
    ## 16cell.4 -1.077760 -3.0507728
    ## 16cell.5  5.018972 -0.9282575

``` r
dim(pca)
```

    ## [1] 268  50

``` r
plotReducedDim(deng_SCE, dimred = "PCA",colour_by="cell_type2",ncomponents = (1:2))
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

You could find PC1 and PC2 captured variances based on different
development stages, let’s take a further look:

``` r
# designed to capture differentiation processes. As a simple measure of pseudotime 
# we can use the coordinates of PC1.
# Plot PC1 vs cell_type2. 
deng_SCE$PC1 <- pca[,1]
deng_SCE$PC2 <- pca[,2]
deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  # rank cells by their PC1 score
ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2, 
                                             colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("PC1") + ylab("Timepoint") +
    ggtitle("Cells ordered by first principal component")
```

    ## Orientation inferred to be along y-axis; override with
    ## `position_quasirandom(orientation = 'x')`

![](session-trajectories_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# Try separating the cell types using other PCs. How does the separation look?
plotReducedDim(deng_SCE, dimred = "PCA",colour_by="cell_type2",ncomponents = (2:3))
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Slingshot map pseudotime

Let us see how another advance trajectory inference method, Slingshot,
performs at placing cells along the expected differentiation trajectory.

``` r
sce <- slingshot(deng_SCE, reducedDim = 'PCA')  # no clusters provided
```

    ## No cluster labels provided. Continuing with one cluster.

``` r
# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
Let’s order cells based on the pseudotime from Slingshot:

``` r
slingshot_df <- data.frame(colData(sce)[, names(colData(sce)) != 'slingshot', drop=FALSE])

# Plot Slingshot pseudotime vs cell stage. 
ggplot(slingshot_df, aes(x = sce$slingPseudotime_1, y = cell_type2, 
                              colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime")
```

    ## Orientation inferred to be along y-axis; override with
    ## `position_quasirandom(orientation = 'x')`

![](session-trajectories_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Next, we are gonna re-clustered the data and run Slingshot again to see
if they still align with the known cluster labels in the data.

``` r
# Cluster cells using the Seurat workflow below.
gcdata <- CreateSeuratObject(counts = counts(deng_SCE), project = "slingshot")
```

    ## Warning: Feature names cannot have pipe characters ('|'), replacing with dashes
    ## ('-')

    ## Warning: Data is of class matrix. Coercing to dgCMatrix.

``` r
gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", scale.factor = 10000)
```

    ## Normalizing layer: counts

``` r
gcdata <- FindVariableFeatures(gcdata, selection.method = "vst", nfeatures = 2000)
```

    ## Finding variable features for layer counts

``` r
gcdata <- ScaleData(object = gcdata, do.center = T, do.scale = F)
```

    ## Centering data matrix

``` r
gcdata <- RunPCA(gcdata, features = VariableFeatures(gcdata), npcs = 40, ndims.print = 1:5, nfeatures.print = 5)
```

    ## PC_ 1 
    ## Positive:  Actb, Fabp3, Psap, Akr1b8, Krt18 
    ## Negative:  Zbed3, C86187, Klf17, Btg4, Ccdc6 
    ## PC_ 2 
    ## Positive:  Krt18, Id2, Akr1b8, BC053393, Fabp3 
    ## Negative:  Gm11517, Alppl2, Obox6, Pdxk, Trim43b 
    ## PC_ 3 
    ## Positive:  Id2, Krt18, Tspan8, BC053393, Krt8 
    ## Negative:  Gm11517, Alppl2, Ypel5, Pdxk, Fam46c 
    ## PC_ 4 
    ## Positive:  Alppl2, Dab2, Gm11517, Krt18, Tspan8 
    ## Negative:  Upp1, Tdgf1, Spp1, Zfp57, Tat 
    ## PC_ 5 
    ## Positive:  Klf17, Ddx24, Bod1l, Tor1b, Gm1995 
    ## Negative:  Alppl2, Gm4340, Gm11756, Gm8300, Gm5039

``` r
# Cluster the cells using the first twenty principal components.
gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:20, k.param = 20)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
gcdata <- FindClusters(gcdata, resolution = 0.8, algorithm = 1, random.seed = 100)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 268
    ## Number of edges: 6814
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7647
    ## Number of communities: 6
    ## Elapsed time: 0 seconds

``` r
# Add clustering information from Seurat to the deng_SCE object
colData(deng_SCE)$Seurat_clusters <- as.character(Idents(gcdata))
```

``` r
plotReducedDim(deng_SCE, dimred = "PCA",color_by="cell_type2",shape_by="Seurat_clusters",ncomponents = (1:2))
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# Then run Slingshot using these cluster assignments.
sce2 <- slingshot(deng_SCE, clusterLabels = 'Seurat_clusters', reducedDim = 'PCA')

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce2)$PCA, col = colors[cut(sce2$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce2), lwd=2, type = 'lineages', col = 'black')
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
## Checking the linkage information
SlingshotDataSet(sce2)
```

    ## class: SlingshotDataSet 
    ## 
    ##  Samples Dimensions
    ##      268         50
    ## 
    ## lineages: 1 
    ## Lineage1: 3  2  4  0  1  5  
    ## 
    ## curves: 1 
    ## Curve1: Length: 392.79   Samples: 268

``` r
slingshot_df2 <- data.frame(colData(sce2)[, names(colData(sce2)) != 'slingshot', drop=FALSE])

slingshot_df2$Seurat_clusters <- factor(slingshot_df2$Seurat_clusters,
                              levels = c("3", "2", "4", "0", 
                                         "1", "5"))

# Plot Slingshot pseudotime vs cell stage. 
ggplot(slingshot_df2, aes(x = slingPseudotime_1, y = Seurat_clusters, 
                              colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime")
```

    ## Orientation inferred to be along y-axis; override with
    ## `position_quasirandom(orientation = 'x')`

![](session-trajectories_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy", "early2cell", "mid2cell", "late2cell", 
                                         "4cell", "8cell", "16cell", "earlyblast",
                                         "midblast", "lateblast"))

slingshot_df2 <- data.frame(colData(sce2)[, names(colData(sce2)) != 'slingshot', drop=FALSE])

# Plot Slingshot pseudotime vs cell stage. 
ggplot(slingshot_df2, aes(x = slingPseudotime_1, y = cell_type2, 
                              colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime")
```

    ## Orientation inferred to be along y-axis; override with
    ## `position_quasirandom(orientation = 'x')`

![](session-trajectories_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

## Identifying temporally dynamic genes

After running slingshot, we are often interested in finding genes that
change their expression over the course of development. We will
demonstrate this type of analysis using the tradeSeq package ([Van den
Berge et al. 2020](https://www.nature.com/articles/s41467-020-14766-3)).

For each gene, we will fit a general additive model (GAM) using a
negative binomial noise distribution to model the (potentially
nonlinear) relationships between gene expression and pseudotime. We will
then test for significant associations between expression and pseudotime
using the associationTest.

Here, we are only chosing ~160 highly variable genes to run the GAM due
to the computational time. You should use all the genes during your
analysis.

``` r
## Selecting highly variable genes
hvg <- modelGeneVar(sce)
chosen <- getTopHVGs(hvg, prop=0.02)
sce <- sce[chosen,]

# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)
```

Next, we visualize these genes in a heatmap with the order of the
pseudotime:

``` r
topgenes <- rownames(ATres[order(ATres$pvalue), ])
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$cell_type2[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Based on the heatmap, you could find temporally dynamic genes across
diffrent development stage, let’s visualize some genes as example:

``` r
p1=plotExpression(deng_SCE, features = "Ccne1", x = "cell_type2",colour_by="cell_type2",exprs_values = "logcounts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2=plotReducedDim(deng_SCE, dimred = "PCA",color_by="Ccne1",ncomponents = (1:2))
p1+p2
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
p1=plotExpression(deng_SCE, features = "Stac2", x = "cell_type2",colour_by="cell_type2",exprs_values = "logcounts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2=plotReducedDim(deng_SCE, dimred = "PCA",color_by="Stac2",ncomponents = (1:2))
p1+p2
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
p1=plotExpression(deng_SCE, features = "Tmem180", x = "cell_type2",colour_by="cell_type2",exprs_values = "logcounts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2=plotReducedDim(deng_SCE, dimred = "PCA",color_by="Tmem180",ncomponents = (1:2))
p1+p2
```

![](session-trajectories_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->
Here, we are using a simple dataset with only 268 cells, which is why we
only got one linkage in the dataset, if you want to explore larger
dataset with multiple linkages with Slingshot, please check
[here](https://nbisweden.github.io/workshop-scRNAseq/labs/trajectory/slingshot.html)

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
    ##  [1] slingshot_2.12.0            TrajectoryUtils_1.12.0     
    ##  [3] princurve_2.1.6             tradeSeq_1.18.0            
    ##  [5] scran_1.32.0                scater_1.32.1              
    ##  [7] scuttle_1.14.0              RColorBrewer_1.1-3         
    ##  [9] cowplot_1.1.3               dplyr_1.1.4                
    ## [11] ggthemes_5.1.0              ggbeeswarm_0.7.2           
    ## [13] ggplot2_3.5.1               Seurat_5.1.0               
    ## [15] SeuratObject_5.0.2          sp_2.1-4                   
    ## [17] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
    ## [19] Biobase_2.64.0              GenomicRanges_1.56.1       
    ## [21] GenomeInfoDb_1.40.1         IRanges_2.38.1             
    ## [23] S4Vectors_0.42.1            BiocGenerics_0.50.0        
    ## [25] MatrixGenerics_1.16.0       matrixStats_1.4.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.22          splines_4.4.1            
    ##   [3] later_1.3.2               tibble_3.2.1             
    ##   [5] polyclip_1.10-7           fastDummies_1.7.4        
    ##   [7] lifecycle_1.0.4           edgeR_4.2.1              
    ##   [9] globals_0.16.3            lattice_0.22-6           
    ##  [11] MASS_7.3-61               magrittr_2.0.3           
    ##  [13] limma_3.60.4              plotly_4.10.4            
    ##  [15] rmarkdown_2.28            yaml_2.3.10              
    ##  [17] metapod_1.12.0            httpuv_1.6.15            
    ##  [19] sctransform_0.4.1         spam_2.10-0              
    ##  [21] spatstat.sparse_3.1-0     reticulate_1.39.0        
    ##  [23] pbapply_1.7-2             abind_1.4-8              
    ##  [25] zlibbioc_1.50.0           Rtsne_0.17               
    ##  [27] purrr_1.0.2               GenomeInfoDbData_1.2.12  
    ##  [29] ggrepel_0.9.6             irlba_2.3.5.1            
    ##  [31] listenv_0.9.1             spatstat.utils_3.1-0     
    ##  [33] goftest_1.2-3             RSpectra_0.16-2          
    ##  [35] spatstat.random_3.3-1     dqrng_0.4.1              
    ##  [37] fitdistrplus_1.2-1        parallelly_1.38.0        
    ##  [39] DelayedMatrixStats_1.26.0 leiden_0.4.3.1           
    ##  [41] codetools_0.2-20          DelayedArray_0.30.1      
    ##  [43] tidyselect_1.2.1          farver_2.1.2             
    ##  [45] UCSC.utils_1.0.0          ScaledMatrix_1.12.0      
    ##  [47] viridis_0.6.5             spatstat.explore_3.3-2   
    ##  [49] jsonlite_1.8.8            BiocNeighbors_1.22.0     
    ##  [51] progressr_0.14.0          ggridges_0.5.6           
    ##  [53] survival_3.7-0            tools_4.4.1              
    ##  [55] ica_1.0-3                 Rcpp_1.0.13              
    ##  [57] glue_1.7.0                gridExtra_2.3            
    ##  [59] SparseArray_1.4.8         xfun_0.47                
    ##  [61] mgcv_1.9-1                withr_3.0.1              
    ##  [63] fastmap_1.2.0             bluster_1.14.0           
    ##  [65] fansi_1.0.6               digest_0.6.37            
    ##  [67] rsvd_1.0.5                R6_2.5.1                 
    ##  [69] mime_0.12                 colorspace_2.1-1         
    ##  [71] scattermore_1.2           tensor_1.5               
    ##  [73] spatstat.data_3.1-2       utf8_1.2.4               
    ##  [75] tidyr_1.3.1               generics_0.1.3           
    ##  [77] data.table_1.16.0         httr_1.4.7               
    ##  [79] htmlwidgets_1.6.4         S4Arrays_1.4.1           
    ##  [81] uwot_0.2.2                pkgconfig_2.0.3          
    ##  [83] gtable_0.3.5              lmtest_0.9-40            
    ##  [85] XVector_0.44.0            htmltools_0.5.8.1        
    ##  [87] dotCall64_1.1-1           scales_1.3.0             
    ##  [89] png_0.1-8                 spatstat.univar_3.0-1    
    ##  [91] knitr_1.48                reshape2_1.4.4           
    ##  [93] nlme_3.1-166              zoo_1.8-12               
    ##  [95] stringr_1.5.1             KernSmooth_2.23-24       
    ##  [97] parallel_4.4.1            miniUI_0.1.1.1           
    ##  [99] vipor_0.4.7               pillar_1.9.0             
    ## [101] grid_4.4.1                vctrs_0.6.5              
    ## [103] RANN_2.6.2                promises_1.3.0           
    ## [105] BiocSingular_1.20.0       beachmat_2.20.0          
    ## [107] xtable_1.8-4              cluster_2.1.6            
    ## [109] beeswarm_0.4.0            evaluate_0.24.0          
    ## [111] cli_3.6.3                 locfit_1.5-9.10          
    ## [113] compiler_4.4.1            rlang_1.1.4              
    ## [115] crayon_1.5.3              future.apply_1.11.2      
    ## [117] labeling_0.4.3            plyr_1.8.9               
    ## [119] stringi_1.8.4             viridisLite_0.4.2        
    ## [121] deldir_2.0-4              BiocParallel_1.38.0      
    ## [123] munsell_0.5.1             lazyeval_0.2.2           
    ## [125] spatstat.geom_3.3-2       Matrix_1.7-0             
    ## [127] RcppHNSW_0.6.0            patchwork_1.2.0          
    ## [129] sparseMatrixStats_1.16.0  future_1.34.0            
    ## [131] statmod_1.5.0             shiny_1.9.1              
    ## [133] highr_0.11                ROCR_1.0-11              
    ## [135] igraph_2.0.3
