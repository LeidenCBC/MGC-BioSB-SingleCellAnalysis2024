Spatial Transcriptomics
================

Created by: Qirong Mao

Edited by: Mohammed Charrout, Claudio Novella-Rausell

# Overview

In this tutorial we will focus on how to analysis the spatial
transcriptomics data. In this exercise we will cover the following
steps:

- Quality control
- Normalization
- Dimentional reduction and clustering
- Identifying spatially variable genes
- Integretion with single-cell RNA-seq data
- Working with multiple slices

This exercise is based on the [Seurat
Vignette](https://satijalab.org/seurat/articles/spatial_vignette.html)
of analysing spatial transcriptomics.

## Datasets

For this tutorial we will use the dataset from the sagittal mouse brain
slices genarated by Visium, including the anterior section and the
matched posterior section. First we will focus of matched anterior and
posterior slides from the sagittal mouse brain. generated Visium
platform. First we will focus on the anterior section.

Load required packages:

``` r
suppressMessages(require(Seurat))
suppressMessages(require(SeuratData))
suppressMessages(require(ggplot2))
suppressMessages(require(patchwork))
suppressMessages(require(dplyr))
```

Load the dataset:

``` r
brain <- LoadData("stxBrain", type = "anterior1")
```

    ## Validating object structure

    ## Updating object slots

    ## Ensuring keys are in the proper structure
    ## Ensuring keys are in the proper structure

    ## Ensuring feature names don't have underscores or pipes

    ## Updating slots in Spatial

    ## Updating slots in anterior1

    ## Validating object structure for Assay5 'Spatial'

    ## Validating object structure for VisiumV2 'anterior1'

    ## Object representation is consistent with the most current Seurat version

``` r
brain
```

    ## An object of class Seurat 
    ## 31053 features across 2696 samples within 1 assay 
    ## Active assay: Spatial (31053 features, 0 variable features)
    ##  1 layer present: counts
    ##  1 spatial field of view present: anterior1

Now we have the assay called “Spatial” instead of “RNA”.

## Quality control

Similar to scRNA-seq analysis pipeline, we can use statistics like
numbers of count, numbers of feature, the percentage of mitochondria and
hemoglobin gene for QC. Please notify that: in the single-cell data, the
unit of these statistics is per single cell. In the Visium data, the
unit is per spot. Each spot contains UMI of 5-20 cells instead of
single-cells.

``` r
# Calculate the percentage of mitochondrial and hemoglobin genes within each spot
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "percent_mito")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")
# Plot the QC-features as violin plots
VlnPlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito",
    "percent_hb"), pt.size = 0.1, ncol = 2) + NoLegend()
```

    ## Warning: Default search for "data" layer in "Spatial" assay yielded no results;
    ## utilizing "counts" layer instead.

![](Spatialtranscriptomics_files/figure-gfm/violin%20plot-1.png)<!-- -->

We can also visualize these statistics on the tissue slide.

``` r
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito","percent_hb"))
```

![](Spatialtranscriptomics_files/figure-gfm/visualize%20tissue%20slide-1.png)<!-- -->

Then we filter out spots with low UMI count/gene number and high
proportion of mitochondrial/hemoglobin genes. You can also choose
different filtering criteria by your own based on your judgement.

``` r
brain = brain[, brain$nCount_Spatial > 2500 & brain$nFeature_Spatial > 500 & 
brain$percent_mito < 25 & brain$percent_hb < 10]
```

Now let us check the data again after filtering:

``` r
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito","percent_hb"))
```

![](Spatialtranscriptomics_files/figure-gfm/visualize%20tissue%20slide%20after%20QC-1.png)<!-- -->

## Normalization: SCTransform

In this tutorial, we use SCTransform method for data normalization.
SCTransform builds regularized negative binomial models of gene
expression in order to account for technical effects while preserving
biological variance. For more details, please check the
[paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1).

``` r
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
```

We can try plot gene expression onto the tissue section

``` r
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
```

![](Spatialtranscriptomics_files/figure-gfm/plotting%20genes-1.png)<!-- -->

## Dimentionality reduction and clustering

After data normalization, we can run dimensionality reduction and
clustering the expression data using the same workflow for analyzing
scRNA-seq:

``` r
#Perform hierarchical clustering using complete linkage & euclidean distance

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
```

    ## 14:47:37 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:47:37 Read 2560 rows and found 30 numeric columns

    ## 14:47:37 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:47:37 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:47:37 Writing NN index file to temp file /var/folders/kq/_ntfn5zx1b199xp4vyd82t840000gn/T//RtmpqlyX5i/file82f6760099c
    ## 14:47:37 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:47:37 Annoy recall = 100%
    ## 14:47:37 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:47:38 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 14:47:38 Commencing optimization for 500 epochs, with 99564 positive edges
    ## 14:47:40 Optimization finished

We can visualize the clustering results on UMAP or onto the tissue
section:

![](Spatialtranscriptomics_files/figure-gfm/hierarchical_eucledian_complete_pcaplot-1.png)<!-- -->

We can also plot each cluster separately for better discrimination

``` r
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 4, 5, 1,
    0, 10)), facet.highlight = TRUE, ncol = 3)
```

![](Spatialtranscriptomics_files/figure-gfm/kmeans-1.png)<!-- -->

### Discussion

Here, we focus on clustering depends on the expression data only. But
for spatial transcriptomics, we can also make use of extra spatial
information and/or tissue slice H&E image. We listed other clustering
methods for spatial transcriptomics:

- [BayesSpace](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763026/)
- [SpaGCN](https://www.nature.com/articles/s41592-021-01255-8)

## Spatially variable genes

Seurat offers two workflows to identify molecular features that
correlate with spatial location within a tissue, so called spatially
variable genes. In this tutorial, we perform differential expression
based on pre-annotated anatomical regions within the tissue, which may
be determined either from unsupervised clustering or prior knowledge.
This strategy works will in this case, as the clusters above exhibit
clear spatial restriction.

``` r
de_markers <- FindMarkers(brain, ident.1 = 1, ident.2 = 0)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[0:6], alpha = c(0.1, 1), ncol = 3)
```

![](Spatialtranscriptomics_files/figure-gfm/detecting%20SVG-1.png)<!-- -->

If you interested more in detecting spatially variable genes, you can
alsoc check other methods like
[SpatialDE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6350895/) and
[SPARK](https://www.nature.com/articles/s41592-019-0701-7) .

## Deconvolution

Spots from the visium assay will encompass the expression profiles of
multiple cells. For the growing list of systems where scRNA-seq data is
available, users may be interested to ‘deconvolute’ each of the spatial
voxels to predict the underlying composition of cell types. In this
tutorial,we use a reference scRNA-seq dataset of ~14,000 adult mouse
cortical cell taxonomy from the Allen Institute, generated with the
SMART-Seq2 protocol.

Since the single cell reference data for deconvolution is from the
cortex only, we need to subset the cortex region from the spatial data
first and then renormalized the subset data:

``` r
# Subset the cortical cells based on cluster labels
cortex <- subset(brain, idents = c(1, 2, 4, 6, 5))
```

``` r
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
p1 + p2
```

![](Spatialtranscriptomics_files/figure-gfm/plot%20the%20deconvolution%20results-1.png)<!-- -->

After processing of spatial data, then we load the pre-processed
single-cell reference data of brain cortex:

``` r
## For downloading the reference dataset: https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
allen_reference <- readRDS("allen_cortex.rds")
```

``` r
allen_reference
```

    ## An object of class Seurat 
    ## 34617 features across 14249 samples within 1 assay 
    ## Active assay: RNA (34617 features, 0 variable features)
    ##  2 layers present: counts, data

``` r
allen_reference@meta.data['subclass']
```

    ##                       subclass
    ## F1S4_160108_001_A01        Vip
    ## F1S4_160108_001_B01      Lamp5
    ## F1S4_160108_001_C01      Lamp5
    ## F1S4_160108_001_D01        Vip
    ## F1S4_160108_001_E01      Lamp5
    ## F1S4_160108_001_F01        Sst
    ## F1S4_160108_001_G01       Sncg
    ## F1S4_160108_001_H01      Lamp5
    ## F1S4_160108_002_A01        Vip
    ## F1S4_160108_002_B01      Lamp5
    ## F1S4_160108_002_C01        Vip
    ## F1S4_160108_002_D01        Vip
    ## F1S4_160108_002_E01        Vip
    ## F1S4_160108_002_F01        Vip
    ## F1S4_160108_002_G01        Vip
    ## F1S4_160108_002_H01        Vip
    ## F1S4_160108_003_A01      Lamp5
    ## F1S4_160108_003_B01      Lamp5
    ## F1S4_160108_003_C01      Lamp5
    ## F1S4_160108_003_D01      Lamp5
    ## F1S4_160108_003_E01       Sncg
    ## F1S4_160108_003_F01        Vip
    ## F1S4_160108_003_G01   Serpinf1
    ## F1S4_160108_003_H01        Vip
    ## F1S4_160108_004_A01      Lamp5
    ## F1S4_160108_004_B01      Lamp5
    ## F1S4_160108_004_C01        Vip
    ## F1S4_160108_004_D01        Vip
    ## F1S4_160108_004_E01        Vip
    ## F1S4_160108_004_F01        Vip
    ## F1S4_160108_004_G01        Vip
    ## F1S4_160108_004_H01      Pvalb
    ## F1S4_160108_005_A01      Lamp5
    ## F1S4_160108_005_B01      Lamp5
    ## F1S4_160108_005_C01      Lamp5
    ## F1S4_160108_005_D01      Lamp5
    ## F1S4_160108_005_E01      Lamp5
    ## F1S4_160108_005_F01      Lamp5
    ## F1S4_160108_005_G01      Lamp5
    ## F1S4_160108_005_H01      Lamp5
    ## F1S4_160108_006_A01      Lamp5
    ## F1S4_160108_006_B01      Lamp5
    ## F1S4_160108_006_C01      Lamp5
    ## F1S4_160108_006_D01        Vip
    ## F1S4_160108_006_E01        Vip
    ## F1S4_160108_006_F01      Lamp5
    ## F1S4_160108_006_G01      Lamp5
    ## F1S4_160108_006_H01        Vip
    ## F1S4_160108_007_A01      Lamp5
    ## F1S4_160108_007_B01      Lamp5
    ## F1S4_160108_007_C01        Vip
    ## F1S4_160108_007_D01        Vip
    ## F1S4_160108_007_E01        Sst
    ## F1S4_160108_007_F01      Pvalb
    ## F1S4_160108_007_G01        Sst
    ## F1S4_160108_007_H01        Vip
    ## F1S4_160108_008_A01        Sst
    ## F1S4_160108_008_B01        Sst
    ## F1S4_160108_008_C01      Lamp5
    ## F1S4_160108_008_D01        Sst
    ## F1S4_160108_008_E01      Pvalb
    ## F1S4_160108_008_F01        Sst
    ## F1S4_160108_008_G01      Pvalb
    ## F1S4_160108_008_H01        Sst
    ## F1S4_160108_009_A01      Pvalb
    ## F1S4_160108_009_B01      Pvalb
    ## F1S4_160108_009_C01      Pvalb
    ## F1S4_160108_009_D01       Endo
    ## F1S4_160108_009_E01        Vip
    ## F1S4_160108_009_F01        Sst
    ## F1S4_160108_009_G01        Sst
    ## F1S4_160108_009_H01        Sst
    ## F1S4_160108_010_A01      Lamp5
    ## F1S4_160108_010_B01      Pvalb
    ## F1S4_160108_010_C01      Lamp5
    ## F1S4_160108_010_D01      Pvalb
    ## F1S4_160108_010_E01        Vip
    ## F1S4_160108_010_F01        Sst
    ## F1S4_160108_010_G01        Vip
    ## F1S4_160108_010_H01        Vip
    ## F1S4_160108_011_A01        Sst
    ## F1S4_160108_011_B01        Sst
    ## F1S4_160108_011_C01        Sst
    ## F1S4_160108_011_D01        Vip
    ## F1S4_160108_011_E01        Sst
    ## F1S4_160108_011_F01      Lamp5
    ## F1S4_160108_011_G01      Pvalb
    ## F1S4_160108_011_H01        Sst
    ## F1S4_160108_012_A01        Vip
    ## F1S4_160108_012_B01        Sst
    ## F1S4_160108_012_C01      Pvalb
    ## F1S4_160108_012_E01        Sst
    ## F1S4_160108_012_F01        Vip
    ## F1S4_160108_012_G01      Lamp5
    ## F1S4_160108_012_H01        Sst
    ## F1S4_160108_013_A01        Vip
    ## F1S4_160108_013_B01        Sst
    ## F1S4_160108_013_C01        Sst
    ## F1S4_160108_013_D01      Pvalb
    ## F1S4_160108_013_E01      Lamp5
    ## F1S4_160108_013_F01        Sst
    ## F1S4_160108_013_G01      Lamp5
    ## F1S4_160108_013_H01       Sncg
    ## F1S4_160108_014_A01      Pvalb
    ## F1S4_160108_014_B01       Peri
    ## F1S4_160108_014_C01        Sst
    ## F1S4_160108_014_D01      Pvalb
    ## F1S4_160108_014_E01        Sst
    ## F1S4_160108_014_F01        Vip
    ## F1S4_160108_014_G01        Sst
    ## F1S4_160108_014_H01      Pvalb
    ## F1S4_161209_001_A01      L6 CT
    ## F1S4_161209_001_B01      L6 CT
    ## F1S4_161209_001_C01      L6 CT
    ## F1S4_161209_001_D01      L6 CT
    ## F1S4_161209_001_E01      L6 CT
    ## F1S4_161209_001_F01      L6 CT
    ## F1S4_161209_001_G01      L6 CT
    ## F1S4_161209_001_H01      L6 CT
    ## F1S4_161209_002_A01      L6 CT
    ## F1S4_161209_002_B01      L6 CT
    ## F1S4_161209_002_C01      L6 CT
    ## F1S4_161209_002_D01      L6 CT
    ## F1S4_161209_002_E01      L6 CT
    ## F1S4_161209_002_F01      L6 CT
    ## F1S4_161209_002_G01      L6 CT
    ## F1S4_161209_002_H01      L6 CT
    ## F1S4_161209_003_A01      L6 CT
    ## F1S4_161209_003_B01      L6 CT
    ## F1S4_161209_003_C01      L6 CT
    ## F1S4_161209_003_D01      L6 CT
    ## F1S4_161209_003_E01      L6 CT
    ## F1S4_161209_003_F01      L6 CT
    ## F1S4_161209_003_G01      L6 CT
    ## F1S4_161209_003_H01      L6 CT
    ## F1S4_161209_004_A01      L6 CT
    ## F1S4_161209_004_B01      L6 CT
    ## F1S4_161209_004_C01      L6 CT
    ## F1S4_161209_004_D01      L6 CT
    ## F1S4_161209_004_E01      L6 CT
    ## F1S4_161209_004_F01      L6 CT
    ## F1S4_161209_004_G01      L6 CT
    ## F1S4_161209_004_H01        L6b
    ## F1S4_161209_005_A01      L6 CT
    ## F1S4_161209_005_B01      L6 CT
    ## F1S4_161209_005_C01      L6 CT
    ## F1S4_161209_005_D01      L6 CT
    ## F1S4_161209_005_E01      L6 CT
    ## F1S4_161209_005_F01      L6 CT
    ## F1S4_161209_005_H01      L6 CT
    ## F1S4_161209_006_A01      L6 CT
    ## F1S4_161209_006_B01      L6 CT
    ## F1S4_161209_006_C01      L6 CT
    ## F1S4_161209_006_D01      L6 CT
    ## F1S4_161209_006_E01      L6 CT
    ## F1S4_161209_006_F01      L6 CT
    ## F1S4_161209_006_G01      L6 CT
    ## F1S4_161209_006_H01      L6 CT
    ## F1S4_161209_007_A01      L6 CT
    ## F1S4_161209_007_B01      L6 CT
    ## F1S4_161209_007_C01      L6 CT
    ## F1S4_161209_007_D01      L6 CT
    ## F1S4_161209_007_E01      L6 CT
    ## F1S4_161209_007_F01      L6 CT
    ## F1S4_161209_007_G01      L6 CT
    ## F1S4_161209_007_H01      L6 CT
    ## F1S4_161209_008_A01      L6 CT
    ## F1S4_161209_008_B01      L6 CT
    ## F1S4_161209_008_C01      L6 CT
    ## F1S4_161209_008_D01      L6 CT
    ## F1S4_161209_008_E01      L6 CT
    ## F1S4_161209_008_F01      L6 CT
    ## F1S4_161209_008_G01      L6 CT
    ## F1S4_161209_008_H01      L6 CT
    ## F1S4_161209_009_A01      L6 CT
    ## F1S4_161209_009_B01      L6 CT
    ## F1S4_161209_009_C01      L6 CT
    ## F1S4_161209_009_D01      L6 CT
    ## F1S4_161209_009_E01      L6 CT
    ## F1S4_161209_009_F01      L6 CT
    ## F1S4_161209_009_G01      L6 CT
    ## F1S4_161209_009_H01      L6 CT
    ## F1S4_161209_010_A01      L6 CT
    ## F1S4_161209_010_B01      L6 CT
    ## F1S4_161209_010_C01        L6b
    ## F1S4_161209_010_D01      L6 CT
    ## F1S4_161209_010_E01      L6 CT
    ## F1S4_161209_010_G01      L6 CT
    ## F1S4_161209_010_H01      L6 CT
    ## F1S4_161209_011_A01      L6 CT
    ## F1S4_161209_011_B01      L6 CT
    ## F1S4_161209_011_C01      L6 CT
    ## F1S4_161209_011_D01      L6 CT
    ## F1S4_161209_011_E01      L6 CT
    ## F1S4_161209_011_F01      L6 CT
    ## F1S4_161209_011_G01      L6 CT
    ## F1S4_161209_011_H01      L6 CT
    ## F1S4_161209_012_A01      L6 CT
    ## F1S4_161209_012_B01      L6 CT
    ## F1S4_161209_012_C01      L6 CT
    ## F1S4_161209_012_D01      L6 CT
    ## F1S4_161209_012_E01      L6 CT
    ## F1S4_161209_012_F01      L6 CT
    ## F1S4_161209_012_H01      L6 CT
    ## F1S4_161209_013_E01      L6 CT
    ## F1S4_161209_013_F01      L6 CT
    ## F1S4_161209_013_G01      L6 CT
    ## F1S4_161209_013_H01      L6 CT
    ## F1S4_170426_005_A01      L6 IT
    ## F1S4_170426_005_B01      L6 IT
    ## F1S4_170426_005_C01      L6 IT
    ## F1S4_170426_005_D01      L6 IT
    ## F1S4_170426_005_E01      L6 IT
    ## F1S4_170426_005_F01      L6 IT
    ## F1S4_170426_005_G01      L6 IT
    ## F1S4_170426_005_H01      L6 IT
    ## F1S4_170426_006_A01      L6 IT
    ## F1S4_170426_006_B01      L6 IT
    ## F1S4_170426_006_C01      L6 IT
    ## F1S4_170426_006_D01      L6 IT
    ## F1S4_170426_006_E01        Sst
    ## F1S4_170426_006_F01      L6 IT
    ## F1S4_170426_006_G01      L6 IT
    ## F1S4_170426_006_H01        Sst
    ## F1S4_170426_009_A01        Sst
    ## F1S4_170426_009_B01        Vip
    ## F1S4_170426_009_C01      Lamp5
    ## F1S4_170426_009_D01        Vip
    ## F1S4_170426_009_E01        Sst
    ## F1S4_170426_009_F01        Vip
    ## F1S4_170426_009_H01    L2/3 IT
    ## F1S4_170426_010_A01      Lamp5
    ## F1S4_170426_010_B01        Vip
    ## F1S4_170426_010_C01    L2/3 IT
    ## F1S4_170426_010_D01        Sst
    ## F1S4_170426_010_E01        Vip
    ## F1S4_170426_010_F01      Lamp5
    ## F1S4_170426_010_G01    L2/3 IT
    ## F1S4_170426_010_H01      Lamp5
    ## F1S4_170426_011_A01        Vip
    ## F1S4_170426_011_B01    L2/3 IT
    ## F1S4_170426_011_C01      Lamp5
    ## F1S4_170426_011_D01        Vip
    ## F1S4_170426_011_E01        Vip
    ## F1S4_170426_011_F01      Lamp5
    ## F1S4_170426_011_G01        Vip
    ## F1S4_170426_011_H01      Lamp5
    ## F1S4_170426_012_A01        Vip
    ## F1S4_170426_012_B01      Pvalb
    ## F1S4_170426_012_C01        Sst
    ## F1S4_170426_012_D01        Vip
    ## F1S4_170426_012_E01         CR
    ## F1S4_170426_012_F01      Lamp5
    ## F1S4_170426_012_G01    L2/3 IT
    ## F1S4_170426_012_H01        Vip
    ## F1S4_170426_013_A01        Vip
    ## F1S4_170426_013_C01    L2/3 IT
    ## F1S4_170426_013_D01      Pvalb
    ## F1S4_170426_013_E01        Vip
    ## F1S4_170426_013_F01      Pvalb
    ## F1S4_170426_013_G01        Vip
    ## F1S4_170426_013_H01        Vip
    ## F1S4_170426_014_A01    L2/3 IT
    ## F1S4_170426_014_B01        Vip
    ## F1S4_170426_014_C01        Vip
    ## F1S4_170426_014_D01    L2/3 IT
    ## F1S4_170426_014_E01        Vip
    ## F1S4_170426_014_G01      Lamp5
    ## F1S4_170426_014_H01        Sst
    ## F1S4_170426_015_A01    L2/3 IT
    ## F1S4_170426_015_B01        Vip
    ## F1S4_170426_015_C01      Pvalb
    ## F1S4_170426_015_D01        Vip
    ## F1S4_170426_015_E01        Vip
    ## F1S4_170426_015_F01    L2/3 IT
    ## F1S4_170426_015_G01      Lamp5
    ## F1S4_170426_015_H01        Vip
    ## F1S4_170426_016_A01      Lamp5
    ## F1S4_170426_016_B01      Pvalb
    ## F1S4_170426_016_C01      Lamp5
    ## F1S4_170426_016_D01      Pvalb
    ## F1S4_170426_016_E01      Pvalb
    ## F1S4_170426_016_F01    L2/3 IT
    ## F1S4_170426_016_G01        Vip
    ## F1S4_170426_016_H01    L2/3 IT
    ## F1S4_170426_017_A01      Pvalb
    ## F1S4_170426_017_B01        Vip
    ## F1S4_170426_017_C01      Lamp5
    ## F1S4_170426_017_D01    L2/3 IT
    ## F1S4_170426_017_E01    L2/3 IT
    ## F1S4_170426_017_F01        Vip
    ## F1S4_170426_017_G01      Pvalb
    ## F1S4_170426_017_H01        Vip
    ## F1S4_170426_022_A01        Sst
    ## F1S4_170426_022_B01      L6 IT
    ## F1S4_170426_022_C01      L6 IT
    ## F1S4_170426_022_E01      L6 IT
    ## F1S4_170426_022_F01      L6 IT
    ## F1S4_170426_022_G01      L6 IT
    ## F1S4_170426_022_H01      L6 IT
    ## F1S4_170426_023_A01      L6 IT
    ## F1S4_170426_023_B01      L6 IT
    ## F1S4_170426_023_C01      L6 IT
    ## F1S4_170426_023_D01      L6 IT
    ## F1S4_170426_023_E01      L6 IT
    ## F1S4_170426_023_F01      L6 IT
    ## F1S4_170426_023_G01      L6 IT
    ## F1S4_170426_024_A01      L6 IT
    ## F1S4_170426_024_C01      L6 IT
    ## F1S4_170426_024_D01      L6 IT
    ## F1S4_170426_024_E01      L6 IT
    ## F1S4_170426_024_F01      L6 IT
    ## F1S4_170426_024_G01      L6 IT
    ## F1S4_170426_024_H01      L6 IT
    ## F1S4_170428_001_A01      L6 IT
    ## F1S4_180108_301_A01      Lamp5
    ## F1S4_180108_301_B01      Lamp5
    ## F1S4_180108_301_C01      Lamp5
    ## F1S4_180108_301_D01        Sst
    ## F1S4_180108_301_E01      Lamp5
    ## F1S4_180108_301_F01        Sst
    ## F1S4_180108_301_G01        Sst
    ## F1S4_180108_301_H01      Lamp5
    ## F1S4_180108_302_A01        Sst
    ## F1S4_180108_302_B01      Lamp5
    ## F1S4_180108_302_C01        Sst
    ## F1S4_180108_302_D01      Lamp5
    ## F1S4_180108_302_E01      Lamp5
    ## F1S4_180108_302_F01        Vip
    ## F1S4_180108_302_G01        Vip
    ## F1S4_180108_302_H01        Sst
    ## F1S4_180108_303_A01        Vip
    ## F1S4_180108_303_B01       Endo
    ## F1S4_180108_303_D01      Lamp5
    ## F1S4_180108_303_E01        Vip
    ## F1S4_180108_303_F01        Sst
    ## F1S4_180108_303_G01      Pvalb
    ## F1S4_180108_303_H01        Sst
    ## F1S4_180108_304_A01        Vip
    ## F1S4_180108_304_B01        Vip
    ## F1S4_180108_304_C01      Lamp5
    ## F1S4_180108_304_D01        Sst
    ## F1S4_180108_304_E01      Lamp5
    ## F1S4_180108_304_F01        Sst
    ## F1S4_180108_304_G01        Vip
    ## F1S4_180108_304_H01      Lamp5
    ## F1S4_180124_301_B01      Lamp5
    ## F1S4_180124_301_C01    L2/3 IT
    ## F1S4_180124_301_D01    L2/3 IT
    ## F1S4_180124_301_E01      Lamp5
    ## F1S4_180124_301_F01      Lamp5
    ## F1S4_180124_301_G01    L2/3 IT
    ## F1S4_180124_302_A01    L2/3 IT
    ## F1S4_180124_302_B01      Lamp5
    ## F1S4_180124_302_C01        Sst
    ## F1S4_180124_302_D01    L2/3 IT
    ## F1S4_180124_302_E01    L2/3 IT
    ## F1S4_180124_302_F01      Lamp5
    ## F1S4_180124_302_G01    L2/3 IT
    ## F1S4_180124_302_H01    L2/3 IT
    ## F1S4_180124_303_A01    L2/3 IT
    ## F1S4_180124_303_B01      Lamp5
    ## F1S4_180124_303_D01      Lamp5
    ## F1S4_180124_303_E01    L2/3 IT
    ## F1S4_180124_303_F01    L2/3 IT
    ## F1S4_180124_303_G01      Lamp5
    ## F1S4_180124_303_H01      Lamp5
    ## F1S4_180124_307_A01      Lamp5
    ## F1S4_180124_307_B01      Lamp5
    ## F1S4_180124_307_C01    L2/3 IT
    ## F1S4_180124_307_D01      Lamp5
    ## F1S4_180124_307_E01    L2/3 IT
    ## F1S4_180124_307_F01    L2/3 IT
    ## F1S4_180124_307_G01       Sncg
    ## F1S4_180124_307_H01       Sncg
    ## F1S4_180124_308_A01      Lamp5
    ## F1S4_180124_308_B01        Vip
    ## F1S4_180124_308_C01       Sncg
    ## F1S4_180124_308_D01        Sst
    ## F1S4_180124_308_E01       Sncg
    ## F1S4_180124_308_F01       Sncg
    ## F1S4_180124_308_G01        Sst
    ## F1S4_180124_308_H01      L6 IT
    ## F1S4_180124_309_A01      Lamp5
    ## F1S4_180124_309_B01       Sncg
    ## F1S4_180124_309_C01        Sst
    ## F1S4_180124_309_D01       Sncg
    ## F1S4_180124_309_F01      Lamp5
    ## F1S4_180124_309_G01      Lamp5
    ## F1S4_180124_309_H01       Sncg
    ## F1S4_180124_310_A01      Lamp5
    ## F1S4_180124_310_B01      Pvalb
    ## F1S4_180124_310_C01      Lamp5
    ## F1S4_180124_310_D01      Lamp5
    ## F1S4_180124_310_E01      Lamp5
    ## F1S4_180124_310_F01      Lamp5
    ## F1S4_180124_310_G01      Lamp5
    ## F1S4_180124_310_H01    L2/3 IT
    ## F1S4_180124_311_A01    L2/3 IT
    ## F1S4_180124_311_B01    L2/3 IT
    ## F1S4_180124_311_C01    L2/3 IT
    ## F1S4_180124_311_D01    L2/3 IT
    ## F1S4_180124_311_E01    L2/3 IT
    ## F1S4_180124_311_F01    L2/3 IT
    ## F1S4_180124_311_G01    L2/3 IT
    ## F1S4_180124_311_H01      Lamp5
    ## F1S4_180124_312_A01    L2/3 IT
    ## F1S4_180124_312_B01    L2/3 IT
    ## F1S4_180124_312_C01    L2/3 IT
    ## F1S4_180124_312_D01      Lamp5
    ## F1S4_180124_312_E01    L2/3 IT
    ## F1S4_180124_312_F01      Lamp5
    ## F1S4_180124_312_H01    L2/3 IT
    ## F1S4_180124_313_A01    L2/3 IT
    ## F1S4_180124_313_B01    L2/3 IT
    ## F1S4_180124_313_C01      Lamp5
    ## F1S4_180124_313_D01      Lamp5
    ## F1S4_180124_313_E01    L2/3 IT
    ## F1S4_180124_313_F01       Sncg
    ## F1S4_180124_313_G01      Lamp5
    ## F1S4_180124_313_H01    L2/3 IT
    ## F1S4_180124_314_B01      Lamp5
    ## F1S4_180124_314_C01    L2/3 IT
    ## F1S4_180124_314_D01    L2/3 IT
    ## F1S4_180124_314_E01    L2/3 IT
    ## F1S4_180124_314_F01    L2/3 IT
    ## F1S4_180124_314_G01    L2/3 IT
    ## F1S4_180124_314_H01    L2/3 IT
    ## F1S4_180129_301_A01      Pvalb
    ## F1S4_180129_301_B01      Pvalb
    ## F1S4_180129_301_C01      Pvalb
    ## F1S4_180129_301_D01      Pvalb
    ## F1S4_180129_301_E01      Pvalb
    ## F1S4_180129_301_F01      Pvalb
    ## F1S4_180129_301_G01      Pvalb
    ## F1S4_180129_301_H01        Sst
    ## F1S4_180129_302_A01      Pvalb
    ## F1S4_180129_302_B01      L5 PT
    ## F1S4_180129_302_C01      Pvalb
    ## F1S4_180129_302_D01      Pvalb
    ## F1S4_180129_302_E01      Pvalb
    ## F1S4_180129_302_F01      Pvalb
    ## F1S4_180129_302_G01        Sst
    ## F1S4_180129_302_H01      Pvalb
    ## F1S4_180129_303_A01      Pvalb
    ## F1S4_180129_303_B01      Pvalb
    ## F1S4_180129_303_C01      Pvalb
    ## F1S4_180129_303_E01      Pvalb
    ## F1S4_180129_303_F01      Pvalb
    ## F1S4_180129_303_G01      Pvalb
    ## F1S4_180129_303_H01        Sst
    ## F1S4_180129_304_A01        Sst
    ## F1S4_180129_304_B01      Pvalb
    ## F1S4_180129_304_C01      Pvalb
    ## F1S4_180129_304_D01        Sst
    ## F1S4_180129_304_E01        Sst
    ## F1S4_180129_304_F01      Pvalb
    ## F1S4_180129_304_G01      Pvalb
    ## F1S4_180129_304_H01        Sst
    ## F1S4_180129_305_A01      Pvalb
    ## F1S4_180129_305_B01        Sst
    ## F1S4_180129_305_C01      Pvalb
    ## F1S4_180129_305_D01      Pvalb
    ## F1S4_180129_305_E01      Pvalb
    ## F1S4_180129_305_F01        Sst
    ## F1S4_180129_305_G01      Pvalb
    ## F1S4_180129_305_H01      Pvalb
    ## F1S4_180129_306_A01      Pvalb
    ## F1S4_180129_306_B01      Pvalb
    ## F1S4_180129_306_C01        Sst
    ## F1S4_180129_306_E01        Sst
    ## F1S4_180129_306_F01      Pvalb
    ## F1S4_180129_306_G01        Sst
    ## F1S4_180129_306_H01      Pvalb
    ## F1S4_180129_307_A01        Sst
    ## F1S4_180129_307_B01      Pvalb
    ## F1S4_180129_307_C01        Sst
    ## F1S4_180129_307_D01        Sst
    ## F1S4_180129_307_E01        Sst
    ## F1S4_180129_307_F01        Sst
    ## F1S4_180129_307_G01        Sst
    ## F1S4_180129_307_H01      Pvalb
    ## F1S4_180129_308_A01      Pvalb
    ## F1S4_180129_308_B01        Sst
    ## F1S4_180129_308_C01        Sst
    ## F1S4_180129_308_D01      Pvalb
    ## F1S4_180129_308_E01      Pvalb
    ## F1S4_180129_308_F01        Sst
    ## F1S4_180129_308_G01        Sst
    ## F1S4_180129_308_H01        Sst
    ## F1S4_180129_309_A01      Pvalb
    ## F1S4_180129_309_B01      Pvalb
    ## F1S4_180129_309_C01      Pvalb
    ## F1S4_180129_309_D01      Pvalb
    ## F1S4_180129_309_E01        Sst
    ## F1S4_180129_309_F01      Pvalb
    ## F1S4_180129_309_G01      Pvalb
    ## F1S4_180129_309_H01        Sst
    ## F1S4_180129_310_A01      Pvalb
    ## F1S4_180129_310_B01      Pvalb
    ## F1S4_180129_310_C01        Sst
    ## F1S4_180129_310_D01      Pvalb
    ## F1S4_180129_310_E01      Pvalb
    ## F1S4_180129_310_F01        Sst
    ## F1S4_180129_310_G01        Sst
    ## F1S4_180129_310_H01        Sst
    ## F1S4_180129_311_A01      Pvalb
    ## F1S4_180129_311_B01      Pvalb
    ## F1S4_180129_311_C01        Sst
    ## F1S4_180129_311_D01        Sst
    ## F1S4_180129_311_E01      Pvalb
    ## F1S4_180129_311_F01        Sst
    ## F1S4_180129_311_G01      Pvalb
    ## F1S4_180129_311_H01      Pvalb
    ## F1S4_180129_312_A01        Sst
    ## F1S4_180129_312_B01        Sst
    ## F1S4_180129_312_C01      Pvalb
    ## F1S4_180129_312_D01        Sst
    ## F1S4_180129_312_E01        Sst
    ## F1S4_180129_312_F01      Pvalb
    ## F1S4_180129_312_G01        Sst
    ## F1S4_180129_312_H01      Pvalb
    ## F1S4_180129_313_A01      Pvalb
    ## F1S4_180129_313_B01      Pvalb
    ## F1S4_180129_313_C01      Pvalb
    ## F1S4_180129_313_D01      Pvalb
    ## F1S4_180129_313_E01      Pvalb
    ## F1S4_180129_313_F01        Sst
    ## F1S4_180129_313_G01      Pvalb
    ## F1S4_180129_313_H01        Sst
    ## F1S4_180129_314_A01      Pvalb
    ## F1S4_180129_314_B01      Pvalb
    ## F1S4_180129_314_C01      Pvalb
    ## F1S4_180129_314_D01      Pvalb
    ## F1S4_180129_314_E01        Sst
    ## F1S4_180129_314_F01      Pvalb
    ## F1S4_180129_314_G01      Pvalb
    ## F1S4_180129_314_H01      Pvalb
    ## F1S4_180129_315_A01        Sst
    ## F1S4_180129_315_B01      Pvalb
    ## F1S4_180129_315_C01        Sst
    ## F1S4_180129_315_D01      Pvalb
    ## F1S4_180129_315_E01      Pvalb
    ## F1S4_180129_315_F01        Sst
    ## F1S4_180129_315_G01        Sst
    ## F1S4_180129_315_H01        Sst
    ## F1S4_180129_316_A01        Sst
    ## F1S4_180129_316_B01        Sst
    ## F1S4_180129_316_C01      Pvalb
    ## F1S4_180129_316_D01      Pvalb
    ## F1S4_180129_316_E01      Pvalb
    ## F1S4_180129_316_F01      Pvalb
    ## F1S4_180129_316_G01      Pvalb
    ## F1S4_180129_316_H01      Pvalb
    ## F1S4_180130_309_A01      Pvalb
    ## F1S4_180130_309_B01        Vip
    ## F1S4_180130_309_C01        Vip
    ## F1S4_180130_309_D01        Vip
    ## F1S4_180130_309_E01        Sst
    ## F1S4_180130_309_F01        Vip
    ## F1S4_180130_309_G01      Lamp5
    ## F1S4_180130_309_H01        Vip
    ## F1S4_180130_310_A01        Vip
    ## F1S4_180130_310_B01      Lamp5
    ## F1S4_180130_310_C01      Lamp5
    ## F1S4_180130_310_D01      Lamp5
    ## F1S4_180130_310_E01      Lamp5
    ## F1S4_180130_310_F01      Pvalb
    ## F1S4_180130_310_G01      Pvalb
    ## F1S4_180130_310_H01      Pvalb
    ## F1S4_180130_311_A01      Pvalb
    ## F1S4_180130_311_B01        Sst
    ## F1S4_180130_311_C01        Sst
    ## F1S4_180130_311_D01        Vip
    ## F1S4_180130_311_E01      Lamp5
    ## F1S4_180130_311_F01      Lamp5
    ## F1S4_180130_311_G01      Lamp5
    ## F1S4_180130_311_H01      Lamp5
    ## F1S4_180130_330_A01      Lamp5
    ## F1S4_180130_330_B01      Lamp5
    ## F1S4_180130_330_C01      Lamp5
    ## F1S4_180130_330_D01      Lamp5
    ## F1S4_180130_330_E01      Lamp5
    ## F1S4_180130_330_F01        Vip
    ## F1S4_180130_330_G01        Vip
    ## F1S4_180130_330_H01      Lamp5
    ## F1S4_180130_331_A01      Lamp5
    ## F1S4_180130_331_B01        Vip
    ## F1S4_180130_331_C01        Sst
    ## F1S4_180130_331_D01      Pvalb
    ## F1S4_180130_331_E01        Sst
    ## F1S4_180130_331_F01        Vip
    ## F1S4_180130_331_G01      Lamp5
    ## F1S4_180130_331_H01      Lamp5
    ## F1S4_180130_332_A01      Lamp5
    ## F1S4_180130_332_B01        Sst
    ## F1S4_180130_332_C01      Lamp5
    ## F1S4_180130_332_D01      Lamp5
    ## F1S4_180130_332_E01        Sst
    ## F1S4_180130_332_F01        Vip
    ## F1S4_180130_332_G01      Lamp5
    ## F1S4_180130_332_H01        Vip
    ## F1S4_180130_333_A01      L6 IT
    ## F1S4_180130_333_B01      L6 IT
    ## F1S4_180130_333_C01      L6 IT
    ## F1S4_180130_333_D01      L6 CT
    ## F1S4_180130_333_E01         NP
    ## F1S4_180130_333_F01      L6 IT
    ## F1S4_180130_333_G01        Sst
    ## F1S4_180130_333_H01      L6 CT
    ## F1S4_180130_334_A01      L6 CT
    ## F1S4_180130_334_B01        Sst
    ## F1S4_180130_334_C01        Sst
    ## F1S4_180130_334_E01      L6 IT
    ## F1S4_180130_334_F01      L6 IT
    ## F1S4_180130_334_G01      L6 IT
    ## F1S4_180130_334_H01        Sst
    ## F1S4_180130_335_A01      L6 IT
    ## F1S4_180130_335_B01      L6 CT
    ## F1S4_180130_335_C01      L6 IT
    ## F1S4_180130_335_D01      Lamp5
    ## F1S4_180130_335_E01      L6 IT
    ## F1S4_180130_335_F01        Sst
    ## F1S4_180130_335_G01        Sst
    ## F1S4_180130_335_H01      L6 IT
    ## F1S4_180130_354_A01      Lamp5
    ## F1S4_180130_354_B01      L6 CT
    ## F1S4_180130_354_C01      L6 IT
    ## F1S4_180130_354_D01      L6 IT
    ## F1S4_180130_354_E01      Lamp5
    ## F1S4_180130_354_F01      L5 PT
    ## F1S4_180130_354_G01         NP
    ## F1S4_180130_354_H01      L6 CT
    ## F1S4_180130_355_A01        Sst
    ## F1S4_180130_355_B01        Sst
    ## F1S4_180130_355_C01      L6 IT
    ## F1S4_180130_355_D01      L6 IT
    ## F1S4_180130_355_E01        Sst
    ## F1S4_180130_355_F01      L6 IT
    ## F1S4_180130_355_G01         NP
    ## F1S4_180130_355_H01      L6 CT
    ## F1S4_180130_356_A01      L5 PT
    ## F1S4_180130_356_B01        Sst
    ## F1S4_180130_356_C01      L6 IT
    ## F1S4_180130_356_D01        Sst
    ## F1S4_180130_356_E01      Lamp5
    ## F1S4_180130_356_F01        Sst
    ## F1S4_180130_356_G01      L6 IT
    ## F1S4_180130_356_H01      Lamp5
    ## F2S4_151217_005_B01      Pvalb
    ## F2S4_151217_005_C01         L4
    ## F2S4_151217_005_E01         L4
    ## F2S4_151217_005_F01         L4
    ## F2S4_151217_005_G01         L4
    ## F2S4_151217_005_H01         L4
    ## F2S4_151217_006_A01         L4
    ## F2S4_151217_006_B01        Vip
    ## F2S4_151217_006_C01    L2/3 IT
    ## F2S4_151217_006_D01         L4
    ## F2S4_151217_006_E01         L4
    ## F2S4_151217_006_F01         L4
    ## F2S4_151217_006_G01         L4
    ## F2S4_151217_006_H01      Lamp5
    ## F2S4_151217_007_A01    L2/3 IT
    ## F2S4_151217_007_B01         L4
    ## F2S4_151217_007_C01         L4
    ## F2S4_151217_007_D01         L4
    ## F2S4_151217_007_E01         L4
    ## F2S4_151217_007_F01         L4
    ## F2S4_151217_007_G01         L4
    ## F2S4_151217_007_H01         L4
    ## F2S4_151217_008_A01         L4
    ## F2S4_151217_008_B01         L4
    ## F2S4_151217_008_C01         L4
    ## F2S4_151217_008_E01         L4
    ## F2S4_151217_008_F01         L4
    ## F2S4_151217_008_G01         L4
    ## F2S4_151217_008_H01         L4
    ## F2S4_151217_009_A01         NP
    ## F2S4_151217_009_B01        Sst
    ## F2S4_151217_009_C01        Vip
    ## F2S4_151217_009_D01      L5 IT
    ## F2S4_151217_009_E01        Vip
    ## F2S4_151217_009_F01        Vip
    ## F2S4_151217_009_G01        Vip
    ## F2S4_151217_009_H01        Vip
    ## F2S4_151217_010_A01        Vip
    ## F2S4_151217_010_B01         NP
    ## F2S4_151217_010_C01      Oligo
    ## F2S4_151217_010_D01         L4
    ## F2S4_151217_010_E01        Vip
    ## F2S4_151217_010_F01         NP
    ## F2S4_151217_010_G01      L5 IT
    ## F2S4_151217_010_H01        Sst
    ## F2S4_151217_011_A01      L6 IT
    ## F2S4_151217_011_B01    L2/3 IT
    ## F2S4_151217_011_C01      L6 IT
    ## F2S4_151217_011_D01      L6 CT
    ## F2S4_151217_011_E01      L6 IT
    ## F2S4_151217_011_G01      L6 IT
    ## F2S4_151217_011_H01      Pvalb
    ## F2S4_151217_012_A01      L6 IT
    ## F2S4_151217_012_B01      Pvalb
    ## F2S4_151217_012_C01      L6 IT
    ## F2S4_151217_012_D01      L6 IT
    ## F2S4_151217_012_E01      L6 IT
    ## F2S4_151217_012_F01      Oligo
    ## F2S4_151217_012_G01         NP
    ## F2S4_151217_012_H01         NP
    ## F2S4_151217_013_A01        Vip
    ## F2S4_151217_013_B01      L6 CT
    ## F2S4_151217_013_C01      L6 IT
    ## F2S4_151217_013_D01      L6 CT
    ## F2S4_151217_013_E01      Lamp5
    ## F2S4_151217_013_F01      L6 IT
    ## F2S4_151217_013_G01         NP
    ## F2S4_151217_013_H01      L6 IT
    ## F2S4_151217_014_A01      L6 IT
    ## F2S4_151217_014_B01      L6 IT
    ## F2S4_151217_014_C01      L6 IT
    ## F2S4_151217_014_D01      L6 IT
    ## F2S4_151217_014_E01      Oligo
    ## F2S4_151217_014_F01        Sst
    ## F2S4_151217_014_G01      L6 IT
    ## F2S4_151217_014_H01      L6 IT
    ## F2S4_151223_001_A01        Vip
    ## F2S4_151223_001_B01        Vip
    ## F2S4_151223_001_C01        Vip
    ## F2S4_151223_001_D01      Lamp5
    ## F2S4_151223_001_E01      Lamp5
    ## F2S4_151223_001_F01      Lamp5
    ## F2S4_151223_001_G01      Lamp5
    ## F2S4_151223_001_H01      Lamp5
    ## F2S4_151223_002_A01        Vip
    ## F2S4_151223_002_C01      Lamp5
    ## F2S4_151223_002_D01        Vip
    ## F2S4_151223_002_E01        Vip
    ## F2S4_151223_002_F01        Vip
    ## F2S4_151223_002_G01        Vip
    ## F2S4_151223_002_H01      Lamp5
    ## F2S4_151223_003_A01      Lamp5
    ## F2S4_151223_003_B01      Lamp5
    ## F2S4_151223_003_C01        Vip
    ## F2S4_151223_003_D01      Lamp5
    ## F2S4_151223_003_E01      Lamp5
    ## F2S4_151223_003_F01        Vip
    ## F2S4_151223_003_G01        Vip
    ## F2S4_151223_003_H01        Vip
    ## F2S4_151223_004_A01      Lamp5
    ## F2S4_151223_004_B01        Vip
    ## F2S4_151223_004_C01      Lamp5
    ## F2S4_151223_004_D01        Vip
    ## F2S4_151223_004_E01        Vip
    ## F2S4_151223_004_F01      Lamp5
    ## F2S4_151223_004_G01      Lamp5
    ## F2S4_151223_004_H01        Vip
    ## F2S4_151223_005_A01        Vip
    ## F2S4_151223_005_B01        Vip
    ## F2S4_151223_005_C01        Vip
    ## F2S4_151223_005_D01        Vip
    ## F2S4_151223_005_E01      Lamp5
    ## F2S4_151223_005_F01      Pvalb
    ## F2S4_151223_005_G01        Vip
    ## F2S4_151223_006_A01        Vip
    ## F2S4_151223_006_B01        Vip
    ## F2S4_151223_006_C01      Lamp5
    ## F2S4_151223_006_D01      Lamp5
    ## F2S4_151223_006_E01        Vip
    ## F2S4_151223_006_F01        Vip
    ## F2S4_151223_006_G01        Vip
    ## F2S4_151223_006_H01        Vip
    ## F2S4_151223_007_A01        Vip
    ## F2S4_151223_007_B01      Lamp5
    ## F2S4_151223_007_C01        Vip
    ## F2S4_151223_007_D01        Vip
    ## F2S4_151223_007_E01        Sst
    ## F2S4_151223_007_F01        Vip
    ## F2S4_151223_007_G01        Vip
    ## F2S4_151223_007_H01        Vip
    ## F2S4_151223_008_A01        Sst
    ## F2S4_151223_008_B01        Sst
    ## F2S4_151223_008_C01   Serpinf1
    ## F2S4_151223_008_D01      Pvalb
    ## F2S4_151223_008_E01        Sst
    ## F2S4_151223_008_F01        Sst
    ## F2S4_151223_008_G01        Vip
    ## F2S4_151223_008_H01        Vip
    ## F2S4_151223_009_A01        Sst
    ## F2S4_151223_009_C01      Pvalb
    ## F2S4_151223_009_D01      Pvalb
    ## F2S4_151223_009_E01      Pvalb
    ## F2S4_151223_009_F01      Pvalb
    ## F2S4_151223_009_G01        Sst
    ## F2S4_151223_009_H01        Vip
    ## F2S4_151223_010_A01      Pvalb
    ## F2S4_151223_010_B01      Pvalb
    ## F2S4_151223_010_C01        Vip
    ## F2S4_151223_010_D01        Sst
    ## F2S4_151223_010_E01      Pvalb
    ## F2S4_151223_010_F01        Vip
    ## F2S4_151223_010_G01      Lamp5
    ## F2S4_151223_011_A01      Pvalb
    ## F2S4_151223_011_B01        Vip
    ## F2S4_151223_011_C01      Pvalb
    ## F2S4_151223_011_D01      Pvalb
    ## F2S4_151223_011_E01      Lamp5
    ## F2S4_151223_011_F01      Lamp5
    ## F2S4_151223_011_G01        Sst
    ## F2S4_151223_011_H01        Sst
    ## F2S4_151223_012_A01        Vip
    ## F2S4_151223_012_B01        Vip
    ## F2S4_151223_012_C01        Sst
    ## F2S4_151223_012_D01        Sst
    ## F2S4_151223_012_E01       Sncg
    ## F2S4_151223_012_F01        Vip
    ## F2S4_151223_012_G01        Sst
    ## F2S4_151223_012_H01        Sst
    ## F2S4_151223_013_A01      Pvalb
    ## F2S4_151223_013_B01      Lamp5
    ## F2S4_151223_013_C01        Sst
    ## F2S4_151223_013_D01        Sst
    ## F2S4_151223_013_E01        Sst
    ## F2S4_151223_013_F01      Oligo
    ## F2S4_151223_013_G01      Pvalb
    ## F2S4_151223_013_H01       Sncg
    ## F2S4_151223_014_A01      Lamp5
    ## F2S4_151223_014_B01      Lamp5
    ## F2S4_151223_014_C01      Pvalb
    ## F2S4_151223_014_D01   Serpinf1
    ## F2S4_151223_014_E01      Lamp5
    ## F2S4_151223_014_F01      Lamp5
    ## F2S4_151223_014_G01        Vip
    ## F2S4_151223_014_H01      Lamp5
    ## F2S4_151223_015_A01        Vip
    ## F2S4_151223_015_B01      Lamp5
    ## F2S4_151223_015_C01      Lamp5
    ## F2S4_151223_015_D01        Vip
    ## F2S4_151223_015_E01      Lamp5
    ## F2S4_151223_015_F01      Lamp5
    ## F2S4_151223_015_G01        Vip
    ## F2S4_151223_015_H01      Lamp5
    ## F2S4_151223_016_A01        Vip
    ## F2S4_151223_016_B01      Lamp5
    ## F2S4_151223_016_C01      Lamp5
    ## F2S4_151223_016_D01        Vip
    ## F2S4_151223_016_E01        Vip
    ## F2S4_151223_016_F01      Lamp5
    ## F2S4_151223_016_G01      Lamp5
    ## F2S4_151223_016_H01      Lamp5
    ## F2S4_151223_017_A01      Lamp5
    ## F2S4_151223_017_B01       Sncg
    ## F2S4_151223_017_C01      Lamp5
    ## F2S4_151223_017_D01       Sncg
    ## F2S4_151223_017_E01      Pvalb
    ## F2S4_151223_017_F01      Lamp5
    ## F2S4_151223_017_G01      Lamp5
    ## F2S4_151223_017_H01      Lamp5
    ## F2S4_151223_018_A01      Lamp5
    ## F2S4_151223_018_B01        Vip
    ## F2S4_151223_018_C01      Lamp5
    ## F2S4_151223_018_D01      Lamp5
    ## F2S4_151223_018_E01      Lamp5
    ## F2S4_151223_018_F01      Lamp5
    ## F2S4_151223_018_G01      Lamp5
    ## F2S4_151223_018_H01        Vip
    ## F2S4_160113_001_A01    L2/3 IT
    ## F2S4_160113_001_B01    L2/3 IT
    ## F2S4_160113_001_C01    L2/3 IT
    ## F2S4_160113_001_D01    L2/3 IT
    ## F2S4_160113_001_E01    L2/3 IT
    ## F2S4_160113_001_F01        Vip
    ## F2S4_160113_001_G01    L2/3 IT
    ## F2S4_160113_001_H01        Vip
    ## F2S4_160113_002_A01    L2/3 IT
    ## F2S4_160113_002_B01    L2/3 IT
    ## F2S4_160113_002_C01    L2/3 IT
    ## F2S4_160113_002_D01    L2/3 IT
    ## F2S4_160113_002_E01    L2/3 IT
    ## F2S4_160113_002_F01    L2/3 IT
    ## F2S4_160113_002_G01    L2/3 IT
    ## F2S4_160113_002_H01        Vip
    ## F2S4_160113_003_A01    L2/3 IT
    ## F2S4_160113_003_B01    L2/3 IT
    ## F2S4_160113_003_C01    L2/3 IT
    ## F2S4_160113_003_D01    L2/3 IT
    ## F2S4_160113_003_E01    L2/3 IT
    ## F2S4_160113_003_F01    L2/3 IT
    ## F2S4_160113_003_G01      Lamp5
    ## F2S4_160113_003_H01    L2/3 IT
    ## F2S4_160113_004_A01    L2/3 IT
    ## F2S4_160113_004_B01    L2/3 IT
    ## F2S4_160113_004_C01    L2/3 IT
    ## F2S4_160113_004_D01    L2/3 IT
    ## F2S4_160113_004_E01    L2/3 IT
    ## F2S4_160113_004_F01    L2/3 IT
    ## F2S4_160113_004_G01    L2/3 IT
    ## F2S4_160113_004_H01      Lamp5
    ## F2S4_160113_005_A01    L2/3 IT
    ## F2S4_160113_005_B01    L2/3 IT
    ## F2S4_160113_005_C01    L2/3 IT
    ## F2S4_160113_005_D01    L2/3 IT
    ## F2S4_160113_005_E01    L2/3 IT
    ## F2S4_160113_005_F01    L2/3 IT
    ## F2S4_160113_005_G01      Lamp5
    ## F2S4_160113_005_H01    L2/3 IT
    ## F2S4_160113_006_A01    L2/3 IT
    ## F2S4_160113_006_B01    L2/3 IT
    ## F2S4_160113_006_C01    L2/3 IT
    ## F2S4_160113_006_D01    L2/3 IT
    ## F2S4_160113_006_E01    L2/3 IT
    ## F2S4_160113_006_F01    L2/3 IT
    ## F2S4_160113_006_G01    L2/3 IT
    ## F2S4_160113_006_H01    L2/3 IT
    ## F2S4_160113_007_B01         L4
    ## F2S4_160113_007_C01      L5 IT
    ## F2S4_160113_007_D01         L4
    ## F2S4_160113_007_E01         L4
    ## F2S4_160113_007_F01         L4
    ## F2S4_160113_007_G01         L4
    ## F2S4_160113_007_H01         L4
    ## F2S4_160113_008_A01         L4
    ## F2S4_160113_008_B01         L4
    ## F2S4_160113_008_C01         L4
    ## F2S4_160113_008_D01         L4
    ## F2S4_160113_008_E01         L4
    ## F2S4_160113_008_F01         L4
    ## F2S4_160113_008_G01         L4
    ## F2S4_160113_008_H01         L4
    ## F2S4_160113_009_A01         L4
    ## F2S4_160113_009_B01         L4
    ## F2S4_160113_009_C01         L4
    ## F2S4_160113_009_D01      L5 IT
    ## F2S4_160113_009_E01         L4
    ## F2S4_160113_009_F01         L4
    ## F2S4_160113_009_G01         L4
    ## F2S4_160113_009_H01         L4
    ## F2S4_160113_010_A01         L4
    ## F2S4_160113_010_B01         L4
    ## F2S4_160113_010_C01      L5 IT
    ## F2S4_160113_010_D01         L4
    ## F2S4_160113_010_E01         L4
    ## F2S4_160113_010_F01         L4
    ## F2S4_160113_010_G01         L4
    ## F2S4_160113_010_H01         L4
    ## F2S4_160113_011_A01         L4
    ## F2S4_160113_011_B01         L4
    ## F2S4_160113_011_C01         L4
    ## F2S4_160113_011_D01         L4
    ## F2S4_160113_011_E01         L4
    ## F2S4_160113_011_F01         L4
    ## F2S4_160113_011_G01         L4
    ## F2S4_160113_011_H01         L4
    ## F2S4_160113_012_A01      L5 IT
    ## F2S4_160113_012_B01         NP
    ## F2S4_160113_012_C01         NP
    ## F2S4_160113_012_E01        Sst
    ## F2S4_160113_012_F01        Sst
    ## F2S4_160113_012_G01         NP
    ## F2S4_160113_012_H01         NP
    ## F2S4_160113_013_A01        Vip
    ## F2S4_160113_013_B01        Vip
    ## F2S4_160113_013_C01      L5 IT
    ## F2S4_160113_013_D01        Vip
    ## F2S4_160113_013_E01        Vip
    ## F2S4_160113_013_F01        Sst
    ## F2S4_160113_013_G01         NP
    ## F2S4_160113_013_H01        Sst
    ## F2S4_160113_014_A01         NP
    ## F2S4_160113_014_B01      L5 PT
    ## F2S4_160113_014_C01        Vip
    ## F2S4_160113_014_D01         NP
    ## F2S4_160113_014_E01      L6 CT
    ## F2S4_160113_014_F01      L6 IT
    ## F2S4_160113_014_G01         NP
    ## F2S4_160113_015_A01      L5 PT
    ## F2S4_160113_015_B01         NP
    ## F2S4_160113_015_C01      L5 IT
    ## F2S4_160113_015_D01      L5 PT
    ## F2S4_160113_015_E01      L6 IT
    ## F2S4_160113_015_F01      L5 IT
    ## F2S4_160113_015_G01      L5 IT
    ## F2S4_160113_015_H01        Sst
    ## F2S4_160113_016_A01      L5 IT
    ## F2S4_160113_016_B01         NP
    ## F2S4_160113_016_C01      L5 IT
    ## F2S4_160113_016_D01        Vip
    ## F2S4_160113_016_E01         NP
    ## F2S4_160113_016_F01         NP
    ## F2S4_160113_016_G01      Lamp5
    ## F2S4_160113_016_H01      L5 IT
    ## F2S4_160113_017_A01        Sst
    ## F2S4_160113_017_C01      L5 PT
    ## F2S4_160113_017_D01         NP
    ## F2S4_160113_017_E01        Vip
    ## F2S4_160113_017_F01        Sst
    ## F2S4_160113_017_G01      L5 IT
    ## F2S4_160113_017_H01        Vip
    ## F2S4_160113_018_A01      L6 CT
    ## F2S4_160113_018_B01      L6 CT
    ## F2S4_160113_018_C01      L6 CT
    ## F2S4_160113_018_D01      L6 CT
    ## F2S4_160113_018_E01      L6 IT
    ## F2S4_160113_018_F01      L6 CT
    ## F2S4_160113_018_G01      L6 IT
    ## F2S4_160113_018_H01      L6 IT
    ## F2S4_160113_019_A01      L6 CT
    ## F2S4_160113_019_B01      L6 CT
    ## F2S4_160113_019_C01      L6 IT
    ## F2S4_160113_019_D01      L6 CT
    ## F2S4_160113_019_E01      L6 CT
    ## F2S4_160113_019_F01      L6 IT
    ## F2S4_160113_019_G01      L6 IT
    ## F2S4_160113_019_H01        Vip
    ## F2S4_160113_020_A01      L6 IT
    ## F2S4_160113_020_B01      L6 IT
    ## F2S4_160113_020_C01      L6 CT
    ## F2S4_160113_020_D01      L6 CT
    ## F2S4_160113_020_E01      L6 CT
    ## F2S4_160113_020_F01      L6 IT
    ## F2S4_160113_020_G01      L6 CT
    ## F2S4_160113_020_H01      L6 IT
    ## F2S4_160113_021_A01      L6 IT
    ## F2S4_160113_021_B01      L6 CT
    ## F2S4_160113_021_C01        Sst
    ## F2S4_160113_021_D01      L6 IT
    ## F2S4_160113_021_E01      L6 IT
    ## F2S4_160113_021_F01      L6 IT
    ## F2S4_160113_021_G01      L6 CT
    ## F2S4_160113_021_H01      L6 IT
    ## F2S4_160113_022_A01      L6 CT
    ## F2S4_160113_022_B01      L6 CT
    ## F2S4_160113_022_C01      L6 CT
    ## F2S4_160113_022_D01      L6 IT
    ## F2S4_160113_022_E01      L6 IT
    ## F2S4_160113_022_F01      L6 CT
    ## F2S4_160113_022_G01        Sst
    ## F2S4_160113_022_H01      L6 IT
    ## F2S4_160113_023_A01      L6 IT
    ## F2S4_160113_023_B01      L6 IT
    ## F2S4_160113_023_C01      L6 IT
    ## F2S4_160113_023_D01      L6 CT
    ## F2S4_160113_023_E01      L6 IT
    ## F2S4_160113_023_F01      L6 IT
    ## F2S4_160113_023_G01      L6 IT
    ## F2S4_160113_023_H01        L6b
    ## F2S4_160113_024_A01      L6 CT
    ## F2S4_160113_024_B01      L6 IT
    ## F2S4_160113_024_C01      L6 IT
    ## F2S4_160113_024_D01      L6 IT
    ## F2S4_160113_024_E01        Vip
    ## F2S4_160113_024_F01      L6 CT
    ## F2S4_160113_024_G01      L6 IT
    ## F2S4_160113_024_H01      L6 IT
    ## F2S4_160113_025_A01      L6 CT
    ## F2S4_160113_025_B01      L6 IT
    ## F2S4_160113_025_C01      L6 CT
    ## F2S4_160113_025_D01      L6 IT
    ## F2S4_160113_025_E01      L6 CT
    ## F2S4_160113_025_F01      L6 CT
    ## F2S4_160113_025_G01      L6 IT
    ## F2S4_160113_025_H01      L6 CT
    ## F2S4_160113_026_A01      L6 IT
    ## F2S4_160113_026_B01        L6b
    ## F2S4_160113_026_C01      L6 IT
    ## F2S4_160113_026_D01      L6 CT
    ## F2S4_160113_026_E01      L6 IT
    ## F2S4_160113_026_F01      L6 IT
    ## F2S4_160113_026_G01      L6 CT
    ## F2S4_160113_026_H01      L6 IT
    ## F2S4_160114_001_A01        Sst
    ## F2S4_160114_001_B01        Sst
    ## F2S4_160114_001_C01        Sst
    ## F2S4_160114_001_D01        Sst
    ## F2S4_160114_001_E01        Sst
    ## F2S4_160114_001_F01        Sst
    ## F2S4_160114_001_G01        Vip
    ## F2S4_160114_001_H01       Sncg
    ## F2S4_160114_002_A01        Sst
    ## F2S4_160114_002_B01      Pvalb
    ## F2S4_160114_002_C01        Sst
    ## F2S4_160114_002_D01        Sst
    ## F2S4_160114_002_E01        Sst
    ## F2S4_160114_002_F01        Sst
    ## F2S4_160114_002_G01        Sst
    ## F2S4_160114_002_H01        Sst
    ## F2S4_160114_003_A01        Sst
    ## F2S4_160114_003_B01      Pvalb
    ## F2S4_160114_003_C01        Sst
    ## F2S4_160114_003_D01        Sst
    ## F2S4_160114_003_E01        Vip
    ## F2S4_160114_003_F01        Sst
    ## F2S4_160114_003_G01      Lamp5
    ## F2S4_160114_003_H01      Pvalb
    ## F2S4_160114_004_A01      Pvalb
    ## F2S4_160114_004_B01      Pvalb
    ## F2S4_160114_004_C01      Pvalb
    ## F2S4_160114_004_D01        Sst
    ## F2S4_160114_004_E01      Pvalb
    ## F2S4_160114_004_F01        Vip
    ## F2S4_160114_004_G01        Sst
    ## F2S4_160114_004_H01        Sst
    ## F2S4_160114_005_A01        Sst
    ## F2S4_160114_005_B01        Sst
    ## F2S4_160114_005_C01        Vip
    ## F2S4_160114_005_D01        Vip
    ## F2S4_160114_005_E01        Vip
    ## F2S4_160114_005_F01        Vip
    ## F2S4_160114_005_G01        Sst
    ## F2S4_160114_006_A01        Vip
    ## F2S4_160114_006_B01      Pvalb
    ## F2S4_160114_006_C01      Lamp5
    ## F2S4_160114_006_D01      Pvalb
    ## F2S4_160114_006_E01        Sst
    ## F2S4_160114_006_F01        Sst
    ## F2S4_160114_006_G01        Sst
    ## F2S4_160114_006_H01        Vip
    ## F2S4_160114_007_A01      Pvalb
    ## F2S4_160114_007_B01        Sst
    ## F2S4_160114_007_E01      Pvalb
    ## F2S4_160114_007_F01        Sst
    ## F2S4_160114_007_G01      Lamp5
    ## F2S4_160114_007_H01      Pvalb
    ## F2S4_160114_008_A01        Vip
    ## F2S4_160114_008_B01      Lamp5
    ## F2S4_160114_008_C01      Pvalb
    ## F2S4_160114_008_D01        Sst
    ## F2S4_160114_008_E01        Sst
    ## F2S4_160114_008_F01        Vip
    ## F2S4_160114_008_G01      Pvalb
    ## F2S4_160114_008_H01        Sst
    ## F2S4_160114_009_A01        Sst
    ## F2S4_160114_009_B01      Pvalb
    ## F2S4_160114_009_C01        Sst
    ## F2S4_160114_009_D01      Pvalb
    ## F2S4_160114_009_E01        Sst
    ## F2S4_160114_009_F01      Meis2
    ## F2S4_160114_009_G01        Sst
    ## F2S4_160114_009_H01      Pvalb
    ## F2S4_160114_010_A01      Lamp5
    ## F2S4_160114_010_B01      Lamp5
    ## F2S4_160114_010_C01        Vip
    ## F2S4_160114_010_D01        Vip
    ## F2S4_160114_010_E01        Vip
    ## F2S4_160114_010_F01      Lamp5
    ## F2S4_160114_010_G01      Lamp5
    ## F2S4_160114_010_H01      Lamp5
    ## F2S4_160114_011_A01        Vip
    ## F2S4_160114_011_B01      Lamp5
    ## F2S4_160114_011_C01      Lamp5
    ## F2S4_160114_011_D01        Vip
    ## F2S4_160114_011_E01      Lamp5
    ## F2S4_160114_011_F01        Vip
    ## F2S4_160114_011_G01      Lamp5
    ## F2S4_160114_011_H01        Vip
    ## F2S4_160114_012_A01       Sncg
    ## F2S4_160114_012_B01      Lamp5
    ## F2S4_160114_012_C01      Lamp5
    ## F2S4_160114_012_D01      Lamp5
    ## F2S4_160114_012_E01      Lamp5
    ## F2S4_160114_012_F01      Lamp5
    ## F2S4_160114_012_G01      Lamp5
    ## F2S4_160114_012_H01        Vip
    ## F2S4_160114_013_A01        Vip
    ## F2S4_160114_013_B01        Vip
    ## F2S4_160114_013_C01        Vip
    ## F2S4_160114_013_D01      Lamp5
    ## F2S4_160114_013_E01        Vip
    ## F2S4_160114_013_F01        Vip
    ## F2S4_160114_013_G01        Vip
    ## F2S4_160114_013_H01       Sncg
    ## F2S4_160114_014_A01        Vip
    ## F2S4_160114_014_B01      Lamp5
    ## F2S4_160114_014_C01      Lamp5
    ## F2S4_160114_014_D01        Vip
    ## F2S4_160114_014_E01       Sncg
    ## F2S4_160114_014_F01      Lamp5
    ## F2S4_160114_014_G01      Lamp5
    ## F2S4_160114_014_H01      Lamp5
    ## F2S4_160114_015_A01        Sst
    ## F2S4_160114_015_C01        Sst
    ## F2S4_160114_015_D01      Pvalb
    ## F2S4_160114_015_E01      Pvalb
    ## F2S4_160114_015_F01      Pvalb
    ## F2S4_160114_015_G01        Sst
    ## F2S4_160114_015_H01      Pvalb
    ## F2S4_160114_016_A01        Vip
    ## F2S4_160114_016_B01      Pvalb
    ## F2S4_160114_016_C01        Sst
    ## F2S4_160114_016_D01        Sst
    ## F2S4_160114_016_E01      Lamp5
    ## F2S4_160114_016_F01        Sst
    ## F2S4_160114_016_G01        Sst
    ## F2S4_160114_016_H01      Pvalb
    ## F2S4_160114_017_A01      Lamp5
    ## F2S4_160114_017_B01        Sst
    ## F2S4_160114_017_C01        Vip
    ## F2S4_160114_017_D01        Sst
    ## F2S4_160114_017_E01        Vip
    ## F2S4_160114_017_F01      Lamp5
    ## F2S4_160114_017_G01        Vip
    ## F2S4_160114_017_H01      Pvalb
    ## F2S4_160114_018_A01        Vip
    ## F2S4_160114_018_B01        Vip
    ## F2S4_160114_018_C01        Sst
    ## F2S4_160114_018_D01      Lamp5
    ## F2S4_160114_018_E01        Sst
    ## F2S4_160114_018_F01        Vip
    ## F2S4_160114_018_G01        Vip
    ## F2S4_160114_018_H01      Lamp5
    ## F2S4_160114_019_A01        Vip
    ## F2S4_160114_019_B01        Vip
    ## F2S4_160114_019_C01        Vip
    ## F2S4_160114_019_D01        Vip
    ## F2S4_160114_019_E01        Sst
    ## F2S4_160114_019_F01      Lamp5
    ## F2S4_160114_019_G01      Lamp5
    ## F2S4_160114_019_H01        Vip
    ## F2S4_160114_020_A01        Vip
    ## F2S4_160114_020_B01      Lamp5
    ## F2S4_160114_020_C01        Vip
    ## F2S4_160114_020_D01      Lamp5
    ## F2S4_160114_020_E01        Vip
    ## F2S4_160114_020_F01      Lamp5
    ## F2S4_160114_020_G01        Vip
    ## F2S4_160114_020_H01        Vip
    ## F2S4_160114_021_A01      Lamp5
    ## F2S4_160114_021_B01        Vip
    ## F2S4_160114_021_C01        Vip
    ## F2S4_160114_021_D01      Lamp5
    ## F2S4_160114_021_E01        Vip
    ## F2S4_160114_021_F01        Vip
    ## F2S4_160114_021_G01      Lamp5
    ## F2S4_160114_021_H01        Vip
    ## F2S4_160115_007_A01        Vip
    ## F2S4_160115_007_B01        Vip
    ## F2S4_160115_007_C01      Lamp5
    ## F2S4_160115_007_D01        Vip
    ## F2S4_160115_007_E01        Vip
    ## F2S4_160115_007_F01        Vip
    ## F2S4_160115_007_G01        Vip
    ## F2S4_160115_007_H01      Lamp5
    ## F2S4_160115_008_A01      Astro
    ## F2S4_160115_008_B01      Astro
    ## F2S4_160115_008_C01      Lamp5
    ## F2S4_160115_008_D01        Vip
    ## F2S4_160115_008_E01        Vip
    ## F2S4_160115_008_F01      Lamp5
    ## F2S4_160115_008_G01        Vip
    ## F2S4_160115_008_H01      Lamp5
    ## F2S4_160115_009_B01        Vip
    ## F2S4_160115_009_C01      Astro
    ## F2S4_160115_009_D01       Sncg
    ## F2S4_160115_009_E01      Lamp5
    ## F2S4_160115_009_F01        Vip
    ## F2S4_160115_009_G01        Vip
    ## F2S4_160115_010_A01        Vip
    ## F2S4_160115_010_B01        Vip
    ## F2S4_160115_010_C01      Lamp5
    ## F2S4_160115_010_D01        Vip
    ## F2S4_160115_010_E01      Lamp5
    ## F2S4_160115_010_F01        Vip
    ## F2S4_160115_010_G01      Lamp5
    ## F2S4_160115_010_H01        Vip
    ## F2S4_160115_011_A01        Sst
    ## F2S4_160115_011_B01      Pvalb
    ## F2S4_160115_011_C01        Vip
    ## F2S4_160115_011_D01        Sst
    ## F2S4_160115_011_E01        Vip
    ## F2S4_160115_011_F01        Vip
    ## F2S4_160115_011_G01        Sst
    ## F2S4_160115_011_H01        Sst
    ## F2S4_160115_012_A01        Sst
    ## F2S4_160115_012_B01      Astro
    ## F2S4_160115_012_C01        Vip
    ## F2S4_160115_012_D01        Vip
    ## F2S4_160115_012_E01        Vip
    ## F2S4_160115_012_F01        Sst
    ## F2S4_160115_012_G01      Astro
    ## F2S4_160115_012_H01        Vip
    ## F2S4_160115_013_A01      Astro
    ## F2S4_160115_013_B01      Astro
    ## F2S4_160115_013_D01        Vip
    ## F2S4_160115_013_E01        Sst
    ## F2S4_160115_013_F01        Sst
    ## F2S4_160115_013_G01      Astro
    ## F2S4_160115_013_H01        Vip
    ## F2S4_160115_014_B01        Vip
    ## F2S4_160115_014_C01      Astro
    ## F2S4_160115_014_D01        Sst
    ## F2S4_160115_014_E01      Lamp5
    ## F2S4_160115_014_F01        Vip
    ## F2S4_160115_014_G01      Astro
    ## F2S4_160115_014_H01        Vip
    ## F2S4_160115_015_C01        Vip
    ## F2S4_160115_015_D01        Sst
    ## F2S4_160115_015_E01      Lamp5
    ## F2S4_160115_015_F01        Sst
    ## F2S4_160115_015_G01        Sst
    ## F2S4_160115_015_H01        Sst
    ## F2S4_160115_016_A01      Pvalb
    ## F2S4_160115_016_B01        Sst
    ## F2S4_160115_016_C01        Sst
    ## F2S4_160115_016_D01        Sst
    ## F2S4_160115_016_E01        Sst
    ## F2S4_160115_016_F01        Sst
    ## F2S4_160115_016_G01      Pvalb
    ## F2S4_160115_017_A01        Sst
    ## F2S4_160115_017_B01        Sst
    ## F2S4_160115_017_C01      Lamp5
    ## F2S4_160115_017_D01        Sst
    ## F2S4_160115_017_E01        Vip
    ## F2S4_160115_017_F01   Serpinf1
    ## F2S4_160115_017_G01        Sst
    ## F2S4_160115_017_H01        Sst
    ## F2S4_160115_018_A01        Vip
    ## F2S4_160115_018_B01        Vip
    ## F2S4_160115_018_C01        Vip
    ## F2S4_160115_018_D01      Lamp5
    ## F2S4_160115_018_E01        Sst
    ## F2S4_160115_018_F01        Sst
    ## F2S4_160115_018_G01        Sst
    ## F2S4_160115_018_H01      Pvalb
    ## F2S4_160115_019_A01      Pvalb
    ## F2S4_160115_019_B01        Sst
    ## F2S4_160115_019_C01        Sst
    ## F2S4_160115_019_D01      Lamp5
    ## F2S4_160115_019_E01        Vip
    ## F2S4_160115_019_F01       Sncg
    ## F2S4_160115_019_G01        Sst
    ## F2S4_160115_020_A01      Lamp5
    ## F2S4_160115_020_B01        Vip
    ## F2S4_160115_020_C01        Sst
    ## F2S4_160115_020_D01        Sst
    ## F2S4_160115_020_E01        Sst
    ## F2S4_160115_020_F01        Vip
    ## F2S4_160115_020_G01        Vip
    ## F2S4_160115_020_H01        Vip
    ## F2S4_160115_021_A01        Sst
    ## F2S4_160115_021_B01        Sst
    ## F2S4_160115_021_C01        Sst
    ## F2S4_160115_021_D01        Sst
    ## F2S4_160115_021_E01      Pvalb
    ## F2S4_160115_021_F01        Sst
    ## F2S4_160115_021_G01      Lamp5
    ## F2S4_160115_021_H01        Sst
    ## F2S4_160115_022_A01        Vip
    ## F2S4_160115_022_B01        Vip
    ## F2S4_160115_022_C01        Sst
    ## F2S4_160115_022_D01        Sst
    ## F2S4_160115_022_E01        Sst
    ## F2S4_160115_022_F01        Vip
    ## F2S4_160115_022_G01        Sst
    ## F2S4_160115_022_H01        Sst
    ## F2S4_160115_023_A01      Lamp5
    ## F2S4_160115_023_B01        Vip
    ## F2S4_160115_023_C01        Sst
    ## F2S4_160115_023_D01      Lamp5
    ## F2S4_160115_023_E01        Vip
    ## F2S4_160115_023_F01        Vip
    ## F2S4_160115_023_G01        Sst
    ## F2S4_160115_023_H01        Sst
    ## F2S4_160115_024_A01      Lamp5
    ## F2S4_160115_024_B01      Lamp5
    ## F2S4_160115_024_C01      Astro
    ## F2S4_160115_024_D01      Astro
    ## F2S4_160115_024_E01        Vip
    ## F2S4_160115_024_F01      Lamp5
    ## F2S4_160115_024_G01      Pvalb
    ## F2S4_160115_024_H01        Vip
    ## F2S4_160115_025_A01        Vip
    ## F2S4_160115_025_B01      Lamp5
    ## F2S4_160115_025_C01      Astro
    ## F2S4_160115_025_D01        Vip
    ## F2S4_160115_025_E01      Lamp5
    ## F2S4_160115_025_F01      Lamp5
    ## F2S4_160115_025_G01       Sncg
    ## F2S4_160115_025_H01        Vip
    ## F2S4_160115_026_A01      Lamp5
    ## F2S4_160115_026_B01      Lamp5
    ## F2S4_160115_026_C01      Lamp5
    ## F2S4_160115_026_D01        Vip
    ## F2S4_160115_026_E01      Lamp5
    ## F2S4_160115_026_F01      Astro
    ## F2S4_160115_026_G01      Lamp5
    ## F2S4_160115_026_H01        Vip
    ## F2S4_160115_027_A01      Lamp5
    ## F2S4_160115_027_B01        Vip
    ## F2S4_160115_027_C01        Vip
    ## F2S4_160115_027_D01      Lamp5
    ## F2S4_160115_027_E01       Sncg
    ## F2S4_160115_027_F01        Sst
    ## F2S4_160115_027_G01      Lamp5
    ## F2S4_160115_027_H01       Sncg
    ## F2S4_160115_028_B01      Lamp5
    ## F2S4_160115_028_C01      Lamp5
    ## F2S4_160115_028_F01      Lamp5
    ## F2S4_160115_028_G01        Vip
    ## F2S4_160115_028_H01      Lamp5
    ## F2S4_160115_029_A01        Vip
    ## F2S4_160115_029_B01      Lamp5
    ## F2S4_160115_029_C01      Lamp5
    ## F2S4_160115_029_D01      Lamp5
    ## F2S4_160115_029_E01        Vip
    ## F2S4_160115_029_F01      Lamp5
    ## F2S4_160115_029_G01        Vip
    ## F2S4_160119_001_A01        Vip
    ## F2S4_160119_001_B01        Vip
    ## F2S4_160119_001_C01        Vip
    ## F2S4_160119_001_D01        Vip
    ## F2S4_160119_001_E01        Vip
    ## F2S4_160119_001_F01    L2/3 IT
    ## F2S4_160119_001_G01        Vip
    ## F2S4_160119_001_H01      Pvalb
    ## F2S4_160119_002_A01 Macrophage
    ## F2S4_160119_002_B01      Pvalb
    ## F2S4_160119_002_C01        Vip
    ## F2S4_160119_002_D01        Vip
    ## F2S4_160119_002_E01        Vip
    ## F2S4_160119_002_F01      Lamp5
    ## F2S4_160119_002_G01      Lamp5
    ## F2S4_160119_002_H01      Lamp5
    ## F2S4_160119_003_A01      Lamp5
    ## F2S4_160119_003_B01        Vip
    ## F2S4_160119_003_C01      Lamp5
    ## F2S4_160119_003_D01        Vip
    ## F2S4_160119_003_E01        Vip
    ## F2S4_160119_003_F01        Vip
    ## F2S4_160119_003_G01       Sncg
    ## F2S4_160119_004_A01        Vip
    ## F2S4_160119_004_B01        Vip
    ## F2S4_160119_004_C01      Pvalb
    ## F2S4_160119_004_D01        Vip
    ## F2S4_160119_004_E01      Lamp5
    ## F2S4_160119_004_F01      Lamp5
    ## F2S4_160119_004_G01        Vip
    ## F2S4_160119_004_H01      Lamp5
    ## F2S4_160119_005_A01        Vip
    ## F2S4_160119_005_B01        Vip
    ## F2S4_160119_005_D01      Lamp5
    ## F2S4_160119_005_E01        Vip
    ## F2S4_160119_005_F01        Vip
    ## F2S4_160119_005_G01      Lamp5
    ## F2S4_160119_005_H01        Vip
    ## F2S4_160119_006_A01      Lamp5
    ## F2S4_160119_006_B01      Lamp5
    ## F2S4_160119_006_C01      Lamp5
    ## F2S4_160119_006_D01      Lamp5
    ## F2S4_160119_006_E01      Lamp5
    ## F2S4_160119_006_F01      Lamp5
    ## F2S4_160119_006_G01      Lamp5
    ## F2S4_160119_006_H01      Lamp5
    ## F2S4_160119_007_A01        Sst
    ## F2S4_160119_007_B01        Sst
    ## F2S4_160119_007_C01        Vip
    ## F2S4_160119_007_D01      Lamp5
    ## F2S4_160119_007_E01      Lamp5
    ## F2S4_160119_007_F01        Vip
    ## F2S4_160119_007_G01      Lamp5
    ## F2S4_160119_007_H01      Lamp5
    ## F2S4_160119_008_A01        Sst
    ## F2S4_160119_008_B01        Sst
    ## F2S4_160119_008_C01        Vip
    ## F2S4_160119_008_D01        Vip
    ## F2S4_160119_008_E01        Vip
    ## F2S4_160119_008_F01        Vip
    ## F2S4_160119_008_G01        Vip
    ## F2S4_160119_008_H01        Vip
    ## F2S4_160119_009_A01      Astro
    ## F2S4_160119_009_B01        Vip
    ## F2S4_160119_009_C01        Vip
    ## F2S4_160119_009_D01      Astro
    ## F2S4_160119_009_E01        Vip
    ## F2S4_160119_009_F01        Vip
    ## F2S4_160119_009_G01        Vip
    ## F2S4_160119_009_H01        Vip
    ## F2S4_160120_001_A01    L2/3 IT
    ## F2S4_160120_001_B01    L2/3 IT
    ## F2S4_160120_001_C01    L2/3 IT
    ## F2S4_160120_001_D01    L2/3 IT
    ## F2S4_160120_001_E01    L2/3 IT
    ## F2S4_160120_001_F01    L2/3 IT
    ## F2S4_160120_001_G01    L2/3 IT
    ## F2S4_160120_001_H01    L2/3 IT
    ## F2S4_160120_002_A01    L2/3 IT
    ## F2S4_160120_002_B01    L2/3 IT
    ## F2S4_160120_002_C01    L2/3 IT
    ## F2S4_160120_002_D01    L2/3 IT
    ## F2S4_160120_002_E01    L2/3 IT
    ## F2S4_160120_002_F01    L2/3 IT
    ## F2S4_160120_002_G01    L2/3 IT
    ## F2S4_160120_002_H01    L2/3 IT
    ## F2S4_160120_003_A01    L2/3 IT
    ## F2S4_160120_003_B01    L2/3 IT
    ## F2S4_160120_003_C01    L2/3 IT
    ## F2S4_160120_003_D01    L2/3 IT
    ## F2S4_160120_003_E01    L2/3 IT
    ## F2S4_160120_003_F01    L2/3 IT
    ## F2S4_160120_003_G01    L2/3 IT
    ## F2S4_160120_003_H01    L2/3 IT
    ## F2S4_160120_004_A01    L2/3 IT
    ## F2S4_160120_004_B01    L2/3 IT
    ## F2S4_160120_004_C01    L2/3 IT
    ## F2S4_160120_004_D01    L2/3 IT
    ## F2S4_160120_004_E01    L2/3 IT
    ## F2S4_160120_004_F01    L2/3 IT
    ## F2S4_160120_004_G01    L2/3 IT
    ## F2S4_160120_004_H01    L2/3 IT
    ## F2S4_160120_005_A01    L2/3 IT
    ## F2S4_160120_005_B01    L2/3 IT
    ## F2S4_160120_005_C01    L2/3 IT
    ## F2S4_160120_005_D01    L2/3 IT
    ## F2S4_160120_005_E01    L2/3 IT
    ## F2S4_160120_005_F01    L2/3 IT
    ## F2S4_160120_005_G01    L2/3 IT
    ## F2S4_160120_005_H01    L2/3 IT
    ## F2S4_160120_006_B01         L4
    ## F2S4_160120_006_C01         L4
    ## F2S4_160120_006_D01         L4
    ## F2S4_160120_006_E01         L4
    ## F2S4_160120_006_F01         L4
    ## F2S4_160120_006_G01         L4
    ## F2S4_160120_006_H01         L4
    ## F2S4_160120_007_A01         L4
    ## F2S4_160120_007_B01         L4
    ## F2S4_160120_007_C01         L4
    ## F2S4_160120_007_D01         L4
    ## F2S4_160120_007_E01         L4
    ## F2S4_160120_007_F01         L4
    ## F2S4_160120_007_G01         L4
    ## F2S4_160120_007_H01      L5 IT
    ## F2S4_160120_008_A01         L4
    ## F2S4_160120_008_B01         L4
    ## F2S4_160120_008_C01         L4
    ## F2S4_160120_008_D01         L4
    ## F2S4_160120_008_E01         L4
    ## F2S4_160120_008_G01         L4
    ## F2S4_160120_008_H01         L4
    ## F2S4_160120_009_A01         L4
    ## F2S4_160120_009_B01         L4
    ## F2S4_160120_009_C01         L4
    ## F2S4_160120_009_D01         L4
    ## F2S4_160120_009_E01         L4
    ## F2S4_160120_009_F01         L4
    ## F2S4_160120_009_G01         L4
    ## F2S4_160120_010_A01         L4
    ## F2S4_160120_010_B01         L4
    ## F2S4_160120_010_C01         L4
    ## F2S4_160120_010_D01         L4
    ## F2S4_160120_010_E01         L4
    ## F2S4_160120_010_F01         L4
    ## F2S4_160120_010_G01         L4
    ## F2S4_160120_011_A01      L6 IT
    ## F2S4_160120_011_B01      L6 IT
    ## F2S4_160120_011_C01         NP
    ## F2S4_160120_011_D01      L6 IT
    ## F2S4_160120_011_E01      L6 IT
    ## F2S4_160120_011_F01      L6 IT
    ## F2S4_160120_011_G01      L6 IT
    ## F2S4_160120_011_H01      L5 IT
    ## F2S4_160120_012_A01      L6 IT
    ## F2S4_160120_012_B01         NP
    ## F2S4_160120_012_C01      L6 IT
    ## F2S4_160120_012_D01      L6 IT
    ## F2S4_160120_012_E01      L5 PT
    ## F2S4_160120_012_F01      L6 IT
    ## F2S4_160120_012_G01      L6 IT
    ## F2S4_160120_012_H01      L5 PT
    ## F2S4_160120_013_A01      L6 IT
    ## F2S4_160120_013_B01      L5 IT
    ## F2S4_160120_013_C01      L6 IT
    ## F2S4_160120_013_D01      L5 PT
    ## F2S4_160120_013_E01      L6 IT
    ## F2S4_160120_013_F01      L6 IT
    ## F2S4_160120_013_G01      L6 IT
    ## F2S4_160120_013_H01      L6 IT
    ## F2S4_160120_014_A01         NP
    ## F2S4_160120_014_B01      L6 IT
    ## F2S4_160120_014_C01      L6 IT
    ## F2S4_160120_014_D01      L6 IT
    ## F2S4_160120_014_E01      L6 IT
    ## F2S4_160120_014_F01      L6 IT
    ## F2S4_160120_014_G01      L6 IT
    ## F2S4_160120_014_H01      L6 IT
    ## F2S4_160120_015_A01      L6 IT
    ## F2S4_160120_015_B01      L6 IT
    ## F2S4_160120_015_C01      L6 IT
    ## F2S4_160120_015_D01      L6 IT
    ## F2S4_160120_015_E01      L6 IT
    ## F2S4_160120_015_F01      L6 IT
    ## F2S4_160120_015_G01      L6 IT
    ## F2S4_160120_015_H01      L6 IT
    ## F2S4_160120_016_A01      L6 IT
    ## F2S4_160120_016_B01      L6 IT
    ## F2S4_160120_016_C01      L6 IT
    ## F2S4_160120_016_D01      L6 IT
    ## F2S4_160120_016_E01      L6 IT
    ## F2S4_160120_016_F01      L6 IT
    ## F2S4_160120_016_G01      L6 IT
    ## F2S4_160120_016_H01      L6 IT
    ## F2S4_160120_017_A01      L6 IT
    ## F2S4_160120_017_B01      L6 IT
    ## F2S4_160120_017_C01      L6 IT
    ## F2S4_160120_017_D01      L6 IT
    ## F2S4_160120_017_E01      L6 IT
    ## F2S4_160120_017_F01      L6 IT
    ## F2S4_160120_017_G01      L6 IT
    ## F2S4_160120_017_H01        L6b
    ## F2S4_160120_018_A01      L6 CT
    ## F2S4_160120_018_B01      L6 IT
    ## F2S4_160120_018_C01      L6 IT
    ## F2S4_160120_018_D01      L6 IT
    ## F2S4_160120_018_E01      L6 IT
    ## F2S4_160120_018_F01      L6 CT
    ## F2S4_160120_018_G01      L6 IT
    ## F2S4_160120_018_H01      L6 CT
    ## F2S4_160120_019_A01      L6 IT
    ## F2S4_160120_019_B01      L6 IT
    ## F2S4_160120_019_C01      L6 IT
    ## F2S4_160120_019_D01      L6 IT
    ## F2S4_160120_019_E01      L6 IT
    ## F2S4_160120_019_F01      L6 CT
    ## F2S4_160120_019_G01      L6 IT
    ## F2S4_160120_019_H01      L6 IT
    ## F2S4_160120_020_A01      L6 IT
    ## F2S4_160120_020_B01      L6 IT
    ## F2S4_160120_020_C01      L6 CT
    ## F2S4_160120_020_D01      L6 IT
    ## F2S4_160120_020_E01      L6 CT
    ## F2S4_160120_020_F01      L6 IT
    ## F2S4_160120_020_G01      L6 CT
    ## F2S4_160120_020_H01      L6 IT
    ## F2S4_160120_021_A01      L6 CT
    ## F2S4_160120_021_B01      L6 CT
    ## F2S4_160120_021_C01      L6 IT
    ## F2S4_160120_021_D01      L6 CT
    ## F2S4_160120_021_E01      L6 IT
    ## F2S4_160120_021_F01      L6 IT
    ## F2S4_160120_021_G01      L6 CT
    ## F2S4_160120_021_H01      L6 CT
    ## F2S4_160120_022_A01      L6 CT
    ## F2S4_160120_022_B01      L6 IT
    ## F2S4_160120_022_C01      L6 CT
    ## F2S4_160120_022_E01      L6 IT
    ## F2S4_160120_022_F01      L6 CT
    ## F2S4_160120_022_G01      L6 IT
    ## F2S4_160120_022_H01      L6 CT
    ## F2S4_160120_023_B01      L6 CT
    ## F2S4_160120_023_C01      L6 IT
    ## F2S4_160120_023_D01      L6 CT
    ## F2S4_160120_023_E01      L6 IT
    ## F2S4_160120_023_F01      L6 IT
    ## F2S4_160120_023_G01      L6 IT
    ## F2S4_160120_023_H01      L6 IT
    ## F2S4_160120_024_A01      L6 CT
    ## F2S4_160120_024_B01      L6 IT
    ## F2S4_160120_024_C01      L6 IT
    ## F2S4_160120_024_D01      L6 IT
    ## F2S4_160120_024_E01      L6 IT
    ## F2S4_160120_024_F01      L6 IT
    ## F2S4_160120_024_G01      L6 IT
    ## F2S4_160120_024_H01      L6 CT
    ## F2S4_160120_025_A01      L6 CT
    ## F2S4_160120_025_B01      L6 IT
    ## F2S4_160120_025_C01      L6 IT
    ## F2S4_160120_025_D01      L6 IT
    ## F2S4_160120_025_E01      L6 CT
    ## F2S4_160120_025_F01      L6 CT
    ## F2S4_160120_025_G01      L6 IT
    ## F2S4_160120_025_H01      L6 CT
    ## F2S4_160120_026_A01      L6 IT
    ## F2S4_160120_026_B01      L6 CT
    ## F2S4_160120_026_C01      L6 IT
    ## F2S4_160120_026_D01      L6 CT
    ## F2S4_160120_026_E01      L6 CT
    ## F2S4_160120_026_F01      L6 CT
    ## F2S4_160120_026_G01      L6 IT
    ## F2S4_160120_026_H01      L6 IT
    ## F2S4_160120_027_A01      L6 IT
    ## F2S4_160120_027_B01      L6 CT
    ## F2S4_160120_027_C01      L6 CT
    ## F2S4_160120_027_D01      L6 IT
    ## F2S4_160120_027_E01      L6 CT
    ## F2S4_160120_027_F01      L6 IT
    ## F2S4_160120_027_G01      L6 CT
    ## F2S4_160120_027_H01      L6 CT
    ## F2S4_160120_028_A01      L6 IT
    ## F2S4_160120_028_B01      L6 IT
    ## F2S4_160120_028_C01      L6 CT
    ## F2S4_160120_028_D01      L6 IT
    ## F2S4_160120_028_E01      L6 IT
    ## F2S4_160120_028_F01      L6 IT
    ## F2S4_160120_028_G01      L6 CT
    ## F2S4_160120_028_H01      L6 CT
    ## F2S4_160121_001_A01    L2/3 IT
    ## F2S4_160121_001_B01        Sst
    ## F2S4_160121_001_C01        Vip
    ## F2S4_160121_001_D01        Vip
    ## F2S4_160121_001_E01      Lamp5
    ## F2S4_160121_001_F01    L2/3 IT
    ## F2S4_160121_001_G01        Vip
    ## F2S4_160121_001_H01      Pvalb
    ## F2S4_160121_002_B01      Pvalb
    ## F2S4_160121_002_C01        Vip
    ## F2S4_160121_002_D01        Vip
    ## F2S4_160121_002_E01        Vip
    ## F2S4_160121_002_F01      Lamp5
    ## F2S4_160121_002_G01      Lamp5
    ## F2S4_160121_002_H01        Vip
    ## F2S4_160121_003_A01      Lamp5
    ## F2S4_160121_003_B01        Sst
    ## F2S4_160121_003_C01        Vip
    ## F2S4_160121_003_D01        Vip
    ## F2S4_160121_003_E01       Sncg
    ## F2S4_160121_003_F01   Serpinf1
    ## F2S4_160121_003_G01      Pvalb
    ## F2S4_160121_003_H01        Sst
    ## F2S4_160121_004_A01      Lamp5
    ## F2S4_160121_004_B01      Lamp5
    ## F2S4_160121_004_C01      Pvalb
    ## F2S4_160121_004_D01      Pvalb
    ## F2S4_160121_004_F01        Vip
    ## F2S4_160121_004_G01        Vip
    ## F2S4_160121_004_H01        Vip
    ## F2S4_160121_005_A01        Vip
    ## F2S4_160121_005_B01        Vip
    ## F2S4_160121_005_C01        Vip
    ## F2S4_160121_005_D01        Vip
    ## F2S4_160121_005_E01        Vip
    ## F2S4_160121_005_F01      Astro
    ## F2S4_160121_005_G01      Pvalb
    ## F2S4_160121_005_H01        Vip
    ## F2S4_160121_006_A01      Lamp5
    ## F2S4_160121_006_B01        Vip
    ## F2S4_160121_006_C01        Vip
    ## F2S4_160121_006_D01      Pvalb
    ## F2S4_160121_006_E01        Vip
    ## F2S4_160121_006_F01        Vip
    ## F2S4_160121_006_G01      Lamp5
    ## F2S4_160121_006_H01        Vip
    ## F2S4_160121_007_A01      Lamp5
    ## F2S4_160121_007_B01      Lamp5
    ## F2S4_160121_007_C01      Lamp5
    ## F2S4_160121_007_D01      Lamp5
    ## F2S4_160121_007_E01      Lamp5
    ## F2S4_160121_007_F01      Lamp5
    ## F2S4_160121_007_G01        Vip
    ## F2S4_160121_007_H01      Astro
    ## F2S4_160121_008_A01       Sncg
    ## F2S4_160121_008_B01      Pvalb
    ## F2S4_160121_008_C01      Lamp5
    ## F2S4_160121_008_D01      Pvalb
    ## F2S4_160121_008_E01      Lamp5
    ## F2S4_160121_008_F01      Lamp5
    ## F2S4_160121_008_G01      Astro
    ## F2S4_160121_008_H01      Pvalb
    ## F2S4_160121_009_A01        Sst
    ## F2S4_160121_009_B01        Vip
    ## F2S4_160121_009_C01        Sst
    ## F2S4_160121_009_D01        Vip
    ## F2S4_160121_009_E01        Vip
    ## F2S4_160121_009_F01        Vip
    ## F2S4_160121_009_G01      Pvalb
    ## F2S4_160121_009_H01      Astro
    ## F2S4_160121_010_A01      Pvalb
    ## F2S4_160121_010_B01      Pvalb
    ## F2S4_160121_010_C01        Sst
    ## F2S4_160121_010_D01        Vip
    ## F2S4_160121_010_E01      Pvalb
    ## F2S4_160121_010_F01        Vip
    ## F2S4_160121_010_G01        Vip
    ## F2S4_160121_010_H01        Vip
    ## F2S4_160121_011_A01        Vip
    ## F2S4_160121_011_B01      Astro
    ## F2S4_160121_011_C01       Sncg
    ## F2S4_160121_011_D01        Vip
    ## F2S4_160121_011_E01        Vip
    ## F2S4_160121_011_F01        Sst
    ## F2S4_160121_011_G01        Vip
    ## F2S4_160121_011_H01      Pvalb
    ## F2S4_160121_012_A01        Sst
    ## F2S4_160121_012_C01        Sst
    ## F2S4_160121_012_D01        Sst
    ## F2S4_160121_012_E01        Vip
    ## F2S4_160121_012_F01        Sst
    ## F2S4_160121_012_G01        Vip
    ## F2S4_160121_012_H01        Sst
    ## F2S4_160121_013_A01        Sst
    ## F2S4_160121_013_B01      Pvalb
    ## F2S4_160121_013_C01      Pvalb
    ## F2S4_160121_013_D01        Vip
    ## F2S4_160121_013_E01        Sst
    ## F2S4_160121_013_F01       Sncg
    ## F2S4_160121_013_G01       Sncg
    ## F2S4_160121_013_H01        Sst
    ## F2S4_160121_014_A01        Sst
    ## F2S4_160121_014_B01        Sst
    ## F2S4_160121_014_C01        Sst
    ## F2S4_160121_014_D01        Sst
    ## F2S4_160121_014_E01        Sst
    ## F2S4_160121_014_F01        Sst
    ## F2S4_160121_014_G01        Vip
    ## F2S4_160121_014_H01      Astro
    ## F2S4_160121_015_A01      Lamp5
    ## F2S4_160121_015_B01      Pvalb
    ## F2S4_160121_015_C01      Lamp5
    ## F2S4_160121_015_D01        Sst
    ## F2S4_160121_015_E01        Sst
    ## F2S4_160121_015_F01        Sst
    ## F2S4_160121_015_G01      Pvalb
    ## F2S4_160121_015_H01      Astro
    ## F2S4_160121_016_A01        Sst
    ## F2S4_160121_016_B01        Sst
    ## F2S4_160121_016_C01      Pvalb
    ## F2S4_160121_016_D01        Sst
    ## F2S4_160121_016_E01      Astro
    ## F2S4_160121_016_F01      Pvalb
    ## F2S4_160121_016_G01      Lamp5
    ## F2S4_160121_016_H01        Sst
    ## F2S4_160121_017_A01      Astro
    ## F2S4_160121_017_B01      Lamp5
    ## F2S4_160121_017_C01      Pvalb
    ## F2S4_160121_017_D01      Pvalb
    ## F2S4_160121_017_E01        Sst
    ## F2S4_160121_017_F01      Lamp5
    ## F2S4_160121_017_G01        Sst
    ## F2S4_160121_017_H01        Sst
    ## F2S4_160121_018_A01        Sst
    ## F2S4_160121_018_B01        Sst
    ## F2S4_160121_018_C01        Sst
    ## F2S4_160121_018_D01        Vip
    ## F2S4_160121_018_E01        Sst
    ## F2S4_160121_018_F01        Sst
    ## F2S4_160121_018_G01        Sst
    ## F2S4_160121_018_H01        Sst
    ## F2S4_160121_019_A01      Lamp5
    ## F2S4_160121_019_B01      Lamp5
    ## F2S4_160121_019_C01        Sst
    ## F2S4_160121_019_D01        Sst
    ## F2S4_160121_019_E01      Pvalb
    ## F2S4_160121_019_F01        Vip
    ## F2S4_160121_019_G01      Lamp5
    ## F2S4_160121_019_H01        Sst
    ## F2S4_160126_005_B01      L5 IT
    ## F2S4_160126_005_C01      L5 IT
    ## F2S4_160126_005_D01         NP
    ## F2S4_160126_005_E01         NP
    ## F2S4_160126_005_F01         NP
    ## F2S4_160126_005_G01      L5 IT
    ## F2S4_160126_005_H01      L6 IT
    ## F2S4_160126_006_A01      L6 IT
    ## F2S4_160126_006_B01         NP
    ## F2S4_160126_006_C01         NP
    ## F2S4_160126_006_D01      L5 PT
    ## F2S4_160126_006_E01      L6 IT
    ## F2S4_160126_006_F01      L6 IT
    ## F2S4_160126_006_G01         NP
    ## F2S4_160126_006_H01      L6 IT
    ## F2S4_160126_007_A01        Sst
    ## F2S4_160126_007_B01      L5 IT
    ## F2S4_160126_007_C01         NP
    ## F2S4_160126_007_D01         NP
    ## F2S4_160126_007_F01      L6 CT
    ## F2S4_160126_007_G01         NP
    ## F2S4_160126_007_H01        Sst
    ## F2S4_160126_008_A01         NP
    ## F2S4_160126_008_B01      L6 IT
    ## F2S4_160126_008_C01         NP
    ## F2S4_160126_008_D01         NP
    ## F2S4_160126_008_E01      L5 IT
    ## F2S4_160126_008_F01         NP
    ## F2S4_160126_008_G01        Sst
    ## F2S4_160126_008_H01         NP
    ## F2S4_160126_009_A01        Vip
    ## F2S4_160126_009_B01      L6 IT
    ## F2S4_160126_009_C01         NP
    ## F2S4_160126_009_D01      Oligo
    ## F2S4_160126_009_E01      L6 CT
    ## F2S4_160126_009_F01         NP
    ## F2S4_160126_009_G01      L5 PT
    ## F2S4_160126_009_H01      L6 IT
    ## F2S4_160126_010_A01      L6 IT
    ## F2S4_160126_010_B01      L6 CT
    ## F2S4_160126_010_C01         NP
    ## F2S4_160126_010_D01         NP
    ## F2S4_160126_010_E01      L6 IT
    ## F2S4_160126_010_F01      L6 IT
    ## F2S4_160126_010_G01      L6 IT
    ## F2S4_160126_010_H01      L5 PT
    ## F2S4_160126_011_A01        Vip
    ## F2S4_160126_011_B01      L6 IT
    ## F2S4_160126_011_C01      L6 CT
    ## F2S4_160126_011_D01      L6 CT
    ## F2S4_160126_011_E01      L6 IT
    ## F2S4_160126_011_F01      L6 IT
    ## F2S4_160126_011_G01      L6 CT
    ## F2S4_160126_011_H01      L6 CT
    ## F2S4_160126_012_A01      L6 IT
    ## F2S4_160126_012_B01      L6 IT
    ## F2S4_160126_012_C01      L6 CT
    ## F2S4_160126_012_D01      L6 IT
    ## F2S4_160126_012_E01        Vip
    ## F2S4_160126_012_F01      L6 CT
    ## F2S4_160126_012_G01      L6 IT
    ## F2S4_160126_012_H01      L6 IT
    ## F2S4_160126_013_A01      L6 IT
    ## F2S4_160126_013_B01      L6 CT
    ## F2S4_160126_013_C01      L6 CT
    ## F2S4_160126_013_D01      L6 IT
    ## F2S4_160126_013_E01      L6 CT
    ## F2S4_160126_013_F01      L6 IT
    ## F2S4_160126_013_G01      L6 CT
    ## F2S4_160126_013_H01      L6 IT
    ## F2S4_160126_014_A01      L6 CT
    ## F2S4_160126_014_B01      L6 IT
    ## F2S4_160126_014_C01      L6 IT
    ## F2S4_160126_014_D01      L6 CT
    ## F2S4_160126_014_E01      L6 IT
    ## F2S4_160126_014_F01      L6 CT
    ## F2S4_160126_014_G01      L6 IT
    ## F2S4_160126_014_H01      L6 IT
    ## F2S4_160126_015_A01        Vip
    ## F2S4_160126_015_B01      L6 CT
    ## F2S4_160126_015_C01      L6 IT
    ## F2S4_160126_015_E01      L6 CT
    ## F2S4_160126_015_F01      L6 IT
    ## F2S4_160126_015_G01      L6 IT
    ## F2S4_160126_015_H01      L6 IT
    ## F2S4_160126_016_A01      L6 CT
    ## F2S4_160126_016_B01        L6b
    ## F2S4_160126_016_C01      L6 CT
    ## F2S4_160126_016_D01      L6 CT
    ## F2S4_160126_016_E01      L6 CT
    ## F2S4_160126_016_F01      L6 IT
    ## F2S4_160126_016_G01      L6 CT
    ## F2S4_160126_016_H01      L6 IT
    ## F2S4_160126_017_A01    L2/3 IT
    ## F2S4_160126_017_B01         L4
    ## F2S4_160126_017_C01        Vip
    ## F2S4_160126_017_D01         L4
    ## F2S4_160126_017_E01    L2/3 IT
    ## F2S4_160126_017_F01    L2/3 IT
    ## F2S4_160126_017_G01         L4
    ## F2S4_160126_017_H01        Vip
    ## F2S4_160126_018_A01         L4
    ## F2S4_160126_018_B01    L2/3 IT
    ## F2S4_160126_018_C01    L2/3 IT
    ## F2S4_160126_018_D01    L2/3 IT
    ## F2S4_160126_018_E01    L2/3 IT
    ## F2S4_160126_018_F01    L2/3 IT
    ## F2S4_160126_018_G01    L2/3 IT
    ## F2S4_160126_018_H01         L4
    ## F2S4_160126_019_A01    L2/3 IT
    ## F2S4_160126_019_B01         L4
    ## F2S4_160126_019_C01    L2/3 IT
    ## F2S4_160126_019_D01    L2/3 IT
    ## F2S4_160126_019_E01        Vip
    ## F2S4_160126_019_F01        Vip
    ## F2S4_160126_019_G01        Vip
    ## F2S4_160126_019_H01         L4
    ## F2S4_160126_020_A01         L4
    ## F2S4_160126_020_B01    L2/3 IT
    ## F2S4_160126_020_C01    L2/3 IT
    ## F2S4_160126_020_D01        Vip
    ## F2S4_160126_020_E01    L2/3 IT
    ## F2S4_160126_020_F01    L2/3 IT
    ## F2S4_160126_020_G01         L4
    ## F2S4_160126_020_H01    L2/3 IT
    ## F2S4_160126_021_A01    L2/3 IT
    ## F2S4_160126_021_B01    L2/3 IT
    ## F2S4_160126_021_C01        Vip
    ## F2S4_160126_021_D01    L2/3 IT
    ## F2S4_160126_021_E01    L2/3 IT
    ## F2S4_160126_021_F01    L2/3 IT
    ## F2S4_160126_021_G01    L2/3 IT
    ## F2S4_160126_021_H01    L2/3 IT
    ## F2S4_160126_022_A01    L2/3 IT
    ## F2S4_160126_022_B01    L2/3 IT
    ## F2S4_160126_022_C01    L2/3 IT
    ## F2S4_160126_022_D01    L2/3 IT
    ## F2S4_160126_022_E01    L2/3 IT
    ## F2S4_160126_022_F01    L2/3 IT
    ## F2S4_160126_022_G01         L4
    ## F2S4_160126_022_H01    L2/3 IT
    ## F2S4_160126_023_A01         L4
    ## F2S4_160126_023_B01         L4
    ## F2S4_160126_023_C01         L4
    ## F2S4_160126_023_D01         L4
    ## F2S4_160126_023_E01         L4
    ## F2S4_160126_023_F01      L5 IT
    ## F2S4_160126_023_H01         L4
    ## F2S4_160126_024_A01         L4
    ## F2S4_160126_024_B01         L4
    ## F2S4_160126_024_C01         L4
    ## F2S4_160126_024_D01         L4
    ## F2S4_160126_024_E01      L5 IT
    ## F2S4_160126_024_F01         L4
    ## F2S4_160126_024_G01         L4
    ## F2S4_160126_024_H01         L4
    ## F2S4_160126_025_A01         L4
    ## F2S4_160126_025_B01         L4
    ## F2S4_160126_025_C01         L4
    ## F2S4_160126_025_D01         L4
    ## F2S4_160126_025_F01         L4
    ## F2S4_160126_025_G01         L4
    ## F2S4_160126_025_H01         L4
    ## F2S4_160126_026_A01         L4
    ## F2S4_160126_026_B01         L4
    ## F2S4_160126_026_D01         L4
    ## F2S4_160126_026_E01        Vip
    ## F2S4_160126_026_F01         L4
    ## F2S4_160126_026_G01         L4
    ## F2S4_160126_027_A01         L4
    ## F2S4_160126_027_B01         L4
    ## F2S4_160126_027_C01      L5 IT
    ## F2S4_160126_027_D01         L4
    ## F2S4_160126_027_E01         L4
    ## F2S4_160126_027_H01         L4
    ## F2S4_160126_028_A01         L4
    ## F2S4_160126_028_C01         L4
    ## F2S4_160126_028_E01         L4
    ## F2S4_160126_028_F01         L4
    ## F2S4_160126_028_G01         L4
    ## F2S4_160126_028_H01         L4
    ## F2S4_160126_029_A01      L6 CT
    ## F2S4_160126_029_B01      L6 IT
    ## F2S4_160126_029_C01      L6 CT
    ## F2S4_160126_029_D01      L6 IT
    ## F2S4_160126_029_E01      L6 IT
    ## F2S4_160126_029_F01      L6 IT
    ## F2S4_160126_029_G01      L6 IT
    ## F2S4_160126_029_H01      L6 CT
    ## F2S4_160126_030_A01      L6 CT
    ## F2S4_160126_030_B01      L6 IT
    ## F2S4_160126_030_C01      L6 CT
    ## F2S4_160126_030_D01      L6 IT
    ## F2S4_160126_030_E01      L6 IT
    ## F2S4_160126_030_F01      L6 IT
    ## F2S4_160126_030_G01      L6 IT
    ## F2S4_160126_030_H01      L6 IT
    ## F2S4_160126_031_A01      L6 IT
    ## F2S4_160126_031_B01      L6 CT
    ## F2S4_160126_031_C01      L6 CT
    ## F2S4_160126_031_D01      L6 IT
    ## F2S4_160126_031_E01      L6 CT
    ## F2S4_160126_031_F01      L6 IT
    ## F2S4_160126_031_G01      L6 CT
    ## F2S4_160126_031_H01      L6 IT
    ## F2S4_160126_032_A01      L6 IT
    ## F2S4_160126_032_B01        Sst
    ## F2S4_160126_032_C01      L6 CT
    ## F2S4_160126_032_D01      L6 CT
    ## F2S4_160126_032_E01      L6 IT
    ## F2S4_160126_032_G01      L6 IT
    ## F2S4_160126_032_H01      L6 CT
    ## F2S4_160126_033_B01      L6 IT
    ## F2S4_160126_033_C01      L6 IT
    ## F2S4_160126_033_D01      Lamp5
    ## F2S4_160126_033_E01      L6 IT
    ## F2S4_160126_033_F01      L6 CT
    ## F2S4_160126_033_G01      L6 IT
    ## F2S4_160126_033_H01      L6 CT
    ## F2S4_160126_034_A01        Sst
    ## F2S4_160126_034_B01      L6 IT
    ## F2S4_160126_034_C01      L6 IT
    ## F2S4_160126_034_D01      L6 IT
    ## F2S4_160126_034_E01      L6 IT
    ## F2S4_160126_034_F01      L6 IT
    ## F2S4_160126_034_G01        Sst
    ## F2S4_160126_034_H01      L6 IT
    ## F2S4_160126_035_A01      L6 IT
    ## F2S4_160126_035_B01      L6 CT
    ## F2S4_160126_035_C01      L6 CT
    ## F2S4_160126_035_D01      L6 CT
    ## F2S4_160126_035_E01      L6 IT
    ## F2S4_160126_035_F01      L6 IT
    ## F2S4_160126_035_G01      L6 CT
    ## F2S4_160126_035_H01      L6 IT
    ## F2S4_160126_036_A01      L6 CT
    ## F2S4_160126_036_B01      L6 IT
    ## F2S4_160126_036_C01      L6 IT
    ## F2S4_160126_036_D01      L6 IT
    ## F2S4_160126_036_E01      L6 CT
    ## F2S4_160126_036_F01      L6 IT
    ## F2S4_160126_036_G01      L6 IT
    ## F2S4_160126_036_H01      L6 IT
    ## F2S4_160127_001_A01    L2/3 IT
    ## F2S4_160127_001_B01        Vip
    ## F2S4_160127_001_C01    L2/3 IT
    ## F2S4_160127_001_D01    L2/3 IT
    ## F2S4_160127_001_E01      Lamp5
    ## F2S4_160127_001_F01        Vip
    ## F2S4_160127_001_G01    L2/3 IT
    ## F2S4_160127_001_H01    L2/3 IT
    ## F2S4_160127_002_A01        Vip
    ## F2S4_160127_002_B01    L2/3 IT
    ## F2S4_160127_002_C01    L2/3 IT
    ## F2S4_160127_002_D01        Vip
    ## F2S4_160127_002_E01    L2/3 IT
    ## F2S4_160127_002_F01    L2/3 IT
    ## F2S4_160127_002_G01      Lamp5
    ## F2S4_160127_002_H01    L2/3 IT
    ## F2S4_160127_003_A01        Vip
    ## F2S4_160127_003_B01    L2/3 IT
    ## F2S4_160127_003_C01       Sncg
    ## F2S4_160127_003_D01      Lamp5
    ## F2S4_160127_003_E01        Vip
    ## F2S4_160127_003_F01    L2/3 IT
    ## F2S4_160127_003_G01    L2/3 IT
    ## F2S4_160127_003_H01    L2/3 IT
    ## F2S4_160127_004_A01    L2/3 IT
    ## F2S4_160127_004_B01        Sst
    ## F2S4_160127_004_C01      Lamp5
    ## F2S4_160127_004_D01    L2/3 IT
    ## F2S4_160127_004_E01    L2/3 IT
    ## F2S4_160127_004_F01    L2/3 IT
    ## F2S4_160127_004_G01    L2/3 IT
    ## F2S4_160127_004_H01      Lamp5
    ## F2S4_160127_005_A01    L2/3 IT
    ## F2S4_160127_005_B01        Vip
    ## F2S4_160127_005_C01    L2/3 IT
    ## F2S4_160127_005_D01    L2/3 IT
    ## F2S4_160127_005_E01    L2/3 IT
    ## F2S4_160127_005_F01    L2/3 IT
    ## F2S4_160127_005_G01       Sncg
    ## F2S4_160127_005_H01      Lamp5
    ## F2S4_160127_006_A01    L2/3 IT
    ## F2S4_160127_006_B01      Lamp5
    ## F2S4_160127_006_C01    L2/3 IT
    ## F2S4_160127_006_D01    L2/3 IT
    ## F2S4_160127_006_E01        Vip
    ## F2S4_160127_006_F01        Vip
    ## F2S4_160127_006_G01    L2/3 IT
    ## F2S4_160127_006_H01    L2/3 IT
    ## F2S4_160127_007_A01         L4
    ## F2S4_160127_007_B01         L4
    ## F2S4_160127_007_C01         L4
    ## F2S4_160127_007_D01        Sst
    ## F2S4_160127_007_E01         L4
    ## F2S4_160127_007_F01        Vip
    ## F2S4_160127_007_G01         L4
    ## F2S4_160127_007_H01         L4
    ## F2S4_160127_008_A01         L4
    ## F2S4_160127_008_B01         L4
    ## F2S4_160127_008_C01         L4
    ## F2S4_160127_008_D01         L4
    ## F2S4_160127_008_E01         L4
    ## F2S4_160127_008_F01         L4
    ## F2S4_160127_008_G01         L4
    ## F2S4_160127_008_H01         L4
    ## F2S4_160127_009_A01         L4
    ## F2S4_160127_009_B01         L4
    ## F2S4_160127_009_C01         L4
    ## F2S4_160127_009_D01         L4
    ## F2S4_160127_009_E01         L4
    ## F2S4_160127_009_F01         L4
    ## F2S4_160127_009_G01         L4
    ## F2S4_160127_009_H01         L4
    ## F2S4_160127_010_A01         L4
    ## F2S4_160127_010_B01         L4
    ## F2S4_160127_010_C01         L4
    ## F2S4_160127_010_D01         L4
    ## F2S4_160127_010_E01         L4
    ## F2S4_160127_010_F01         L4
    ## F2S4_160127_010_G01         L4
    ## F2S4_160127_010_H01         L4
    ## F2S4_160127_011_A01         L4
    ## F2S4_160127_011_B01         L4
    ## F2S4_160127_011_C01         L4
    ## F2S4_160127_011_D01         L4
    ## F2S4_160127_011_E01    L2/3 IT
    ## F2S4_160127_011_F01         L4
    ## F2S4_160127_011_G01         L4
    ## F2S4_160127_011_H01         L4
    ## F2S4_160127_012_A01         L4
    ## F2S4_160127_012_B01         L4
    ## F2S4_160127_012_C01         L4
    ## F2S4_160127_012_D01         L4
    ## F2S4_160127_012_E01         L4
    ## F2S4_160127_012_F01         L4
    ## F2S4_160127_012_G01         L4
    ## F2S4_160127_012_H01         L4
    ## F2S4_160127_013_A01      L5 IT
    ## F2S4_160127_013_B01      L6 IT
    ## F2S4_160127_013_C01         NP
    ## F2S4_160127_013_D01      L5 PT
    ## F2S4_160127_013_E01      L6 IT
    ## F2S4_160127_013_F01      L6 IT
    ## F2S4_160127_013_G01         NP
    ## F2S4_160127_013_H01        Sst
    ## F2S4_160127_014_A01         NP
    ## F2S4_160127_014_B01         NP
    ## F2S4_160127_014_D01        Sst
    ## F2S4_160127_014_E01      L5 IT
    ## F2S4_160127_014_F01        Vip
    ## F2S4_160127_014_G01        Sst
    ## F2S4_160127_014_H01      L6 IT
    ## F2S4_160127_015_A01      L6 IT
    ## F2S4_160127_015_B01      L6 IT
    ## F2S4_160127_015_C01      L6 IT
    ## F2S4_160127_015_D01      L6 IT
    ## F2S4_160127_015_E01         NP
    ## F2S4_160127_015_F01         NP
    ## F2S4_160127_015_G01         NP
    ## F2S4_160127_015_H01         NP
    ## F2S4_160127_016_A01        Sst
    ## F2S4_160127_016_B01      L5 PT
    ## F2S4_160127_016_C01      L5 PT
    ## F2S4_160127_016_D01      L5 PT
    ## F2S4_160127_016_E01      L5 IT
    ## F2S4_160127_016_F01      Pvalb
    ## F2S4_160127_016_G01         NP
    ## F2S4_160127_016_H01        Vip
    ## F2S4_160127_017_A01         NP
    ## F2S4_160127_017_B01         L4
    ## F2S4_160127_017_C01        Sst
    ## F2S4_160127_017_D01        Sst
    ## F2S4_160127_017_E01        Vip
    ## F2S4_160127_017_F01         NP
    ## F2S4_160127_017_G01      Lamp5
    ## F2S4_160127_017_H01      L5 IT
    ## F2S4_160127_018_A01      L5 PT
    ## F2S4_160127_018_B01      L5 IT
    ## F2S4_160127_018_C01      L5 PT
    ## F2S4_160127_018_D01        Vip
    ## F2S4_160127_018_E01      Pvalb
    ## F2S4_160127_018_F01      L5 PT
    ## F2S4_160127_018_G01      L6 IT
    ## F2S4_160127_018_H01        Sst
    ## F2S4_160127_019_A01      L6 IT
    ## F2S4_160127_019_B01      L6 CT
    ## F2S4_160127_019_C01      L6 IT
    ## F2S4_160127_019_D01      L6 IT
    ## F2S4_160127_019_E01      L6 CT
    ## F2S4_160127_019_F01      L6 CT
    ## F2S4_160127_019_G01      L6 IT
    ## F2S4_160127_019_H01        Sst
    ## F2S4_160127_020_E01      L6 IT
    ## F2S4_160127_020_G01      L6 IT
    ## F2S4_160127_020_H01      L6 IT
    ## F2S4_160127_021_A01      L6 IT
    ## F2S4_160127_021_B01      L6 IT
    ## F2S4_160127_021_C01      L6 CT
    ## F2S4_160127_021_D01      L6 IT
    ## F2S4_160127_021_E01      L6 IT
    ## F2S4_160127_021_F01      L6 IT
    ## F2S4_160127_021_G01      L6 CT
    ## F2S4_160127_021_H01      L6 CT
    ## F2S4_160127_022_A01      L6 IT
    ## F2S4_160127_022_B01      L6 IT
    ## F2S4_160127_022_C01      L6 IT
    ## F2S4_160127_022_D01      L6 IT
    ## F2S4_160127_022_E01      L6 IT
    ## F2S4_160127_022_F01      L6 IT
    ## F2S4_160127_022_G01      L6 IT
    ## F2S4_160127_022_H01      L6 IT
    ## F2S4_160127_023_A01      L6 IT
    ## F2S4_160127_023_B01      L6 IT
    ## F2S4_160127_023_C01      L6 CT
    ## F2S4_160127_023_D01      L6 IT
    ## F2S4_160127_023_E01      Lamp5
    ## F2S4_160127_023_F01      L6 IT
    ## F2S4_160127_023_G01      L6 CT
    ## F2S4_160127_023_H01      L6 IT
    ## F2S4_160127_024_A01      L6 IT
    ## F2S4_160127_024_B01      L6 IT
    ## F2S4_160127_024_C01      L6 IT
    ## F2S4_160127_024_D01      L6 IT
    ## F2S4_160127_024_E01      L6 IT
    ## F2S4_160127_024_F01      L6 CT
    ## F2S4_160127_024_G01      L6 IT
    ## F2S4_160127_024_H01      L6 IT
    ## F2S4_160127_025_A01      L6 CT
    ## F2S4_160127_025_B01      L6 IT
    ## F2S4_160127_025_C01      L6 CT
    ## F2S4_160127_025_D01      L6 CT
    ## F2S4_160127_025_E01      L6 CT
    ## F2S4_160127_025_G01      L6 IT
    ## F2S4_160127_025_H01      L6 CT
    ## F2S4_160127_026_A01      L6 CT
    ## F2S4_160127_026_B01      L6 IT
    ## F2S4_160127_026_C01      L6 IT
    ## F2S4_160127_026_D01      L6 IT
    ## F2S4_160127_026_F01      L6 IT
    ## F2S4_160127_026_G01      L6 IT
    ## F2S4_160127_026_H01      L6 IT
    ## F2S4_160127_027_A01      L6 IT
    ## F2S4_160127_027_B01      L6 IT
    ## F2S4_160127_027_C01      L6 CT
    ## F2S4_160127_027_D01      L6 IT
    ## F2S4_160127_027_E01      L6 IT
    ## F2S4_160127_027_F01      L6 CT
    ## F2S4_160127_027_G01      L6 IT
    ## F2S4_160127_027_H01      L6 IT
    ## F2S4_160127_028_A01      L6 CT
    ## F2S4_160127_028_B01      L6 IT
    ## F2S4_160127_028_C01      L6 CT
    ## F2S4_160127_028_D01      L6 IT
    ## F2S4_160127_028_E01      L6 IT
    ## F2S4_160127_028_F01      L6 IT
    ## F2S4_160127_028_G01      L6 CT
    ## F2S4_160127_028_H01      L6 IT
    ## F2S4_160127_029_A01      L6 IT
    ## F2S4_160127_029_B01      L6 CT
    ## F2S4_160127_029_C01      L6 IT
    ## F2S4_160127_029_D01      L6 CT
    ## F2S4_160127_029_E01      L6 CT
    ## F2S4_160127_029_F01      L6 IT
    ## F2S4_160127_029_G01      L6 IT
    ## F2S4_160127_029_H01      L5 IT
    ## F2S4_160127_030_A01      L6 CT
    ## F2S4_160127_030_B01      L6 IT
    ## F2S4_160127_030_C01      L6 IT
    ## F2S4_160127_030_D01      L6 IT
    ## F2S4_160127_030_E01      L6 CT
    ## F2S4_160127_030_F01      L6 IT
    ## F2S4_160127_030_G01      L6 CT
    ## F2S4_160127_030_H01        L6b
    ## F2S4_160128_001_A01    L2/3 IT
    ## F2S4_160128_001_B01    L2/3 IT
    ## F2S4_160128_001_C01        Vip
    ## F2S4_160128_001_D01    L2/3 IT
    ## F2S4_160128_001_E01        Vip
    ## F2S4_160128_001_F01    L2/3 IT
    ## F2S4_160128_001_G01        Vip
    ## F2S4_160128_001_H01      Lamp5
    ## F2S4_160128_002_A01    L2/3 IT
    ## F2S4_160128_002_B01    L2/3 IT
    ## F2S4_160128_002_C01      Lamp5
    ## F2S4_160128_002_D01    L2/3 IT
    ## F2S4_160128_002_E01        Vip
    ## F2S4_160128_002_F01    L2/3 IT
    ## F2S4_160128_002_G01        Sst
    ## F2S4_160128_002_H01    L2/3 IT
    ## F2S4_160128_003_A01      Lamp5
    ## F2S4_160128_003_B01        Vip
    ## F2S4_160128_003_C01    L2/3 IT
    ## F2S4_160128_003_D01      Lamp5
    ## F2S4_160128_003_E01    L2/3 IT
    ## F2S4_160128_003_F01    L2/3 IT
    ## F2S4_160128_003_G01        Vip
    ## F2S4_160128_003_H01      Lamp5
    ## F2S4_160128_004_A01        Vip
    ## F2S4_160128_004_B01         L4
    ## F2S4_160128_004_C01    L2/3 IT
    ## F2S4_160128_004_D01        Vip
    ## F2S4_160128_004_E01    L2/3 IT
    ## F2S4_160128_004_F01    L2/3 IT
    ## F2S4_160128_004_G01         L4
    ## F2S4_160128_004_H01    L2/3 IT
    ## F2S4_160128_005_A01    L2/3 IT
    ## F2S4_160128_005_B01        Vip
    ## F2S4_160128_005_C01        Vip
    ## F2S4_160128_005_D01    L2/3 IT
    ## F2S4_160128_005_E01      Lamp5
    ## F2S4_160128_005_F01         L4
    ## F2S4_160128_005_G01    L2/3 IT
    ## F2S4_160128_005_H01    L2/3 IT
    ## F2S4_160128_006_A01        Vip
    ## F2S4_160128_006_B01    L2/3 IT
    ## F2S4_160128_006_C01    L2/3 IT
    ## F2S4_160128_006_D01    L2/3 IT
    ## F2S4_160128_006_E01        Vip
    ## F2S4_160128_006_F01    L2/3 IT
    ## F2S4_160128_006_G01    L2/3 IT
    ## F2S4_160128_006_H01    L2/3 IT
    ## F2S4_160128_007_A01         L4
    ## F2S4_160128_007_B01         L4
    ## F2S4_160128_007_C01         L4
    ## F2S4_160128_007_D01         L4
    ## F2S4_160128_007_E01         L4
    ## F2S4_160128_007_F01         L4
    ## F2S4_160128_007_G01      L5 PT
    ## F2S4_160128_007_H01         L4
    ## F2S4_160128_008_B01         L4
    ## F2S4_160128_008_C01         L4
    ## F2S4_160128_008_E01      L5 IT
    ## F2S4_160128_008_F01      L5 IT
    ## F2S4_160128_008_G01         L4
    ## F2S4_160128_008_H01        Vip
    ## F2S4_160128_009_A01         L4
    ## F2S4_160128_009_B01         L4
    ## F2S4_160128_009_C01      L5 IT
    ## F2S4_160128_009_D01         L4
    ## F2S4_160128_009_F01      L5 IT
    ## F2S4_160128_009_G01        Vip
    ## F2S4_160128_009_H01         L4
    ## F2S4_160128_010_A01         L4
    ## F2S4_160128_010_B01         L4
    ## F2S4_160128_010_C01        Sst
    ## F2S4_160128_010_D01         L4
    ## F2S4_160128_010_E01         L4
    ## F2S4_160128_010_F01         L4
    ## F2S4_160128_010_H01         NP
    ## F2S4_160128_011_A01         L4
    ## F2S4_160128_011_B01         L4
    ## F2S4_160128_011_C01         L4
    ## F2S4_160128_011_D01         L4
    ## F2S4_160128_011_E01         L4
    ## F2S4_160128_011_F01         L4
    ## F2S4_160128_011_G01         L4
    ## F2S4_160128_011_H01         L4
    ## F2S4_160128_012_A01         L4
    ## F2S4_160128_012_B01         L4
    ## F2S4_160128_012_C01         L4
    ## F2S4_160128_012_D01         L4
    ## F2S4_160128_012_E01         L4
    ## F2S4_160128_012_F01         L4
    ## F2S4_160128_012_G01         NP
    ## F2S4_160128_012_H01         L4
    ## F2S4_160128_013_A01      L5 IT
    ## F2S4_160128_013_B01      L5 IT
    ## F2S4_160128_013_C01      L6 IT
    ## F2S4_160128_013_D01      L6 IT
    ## F2S4_160128_013_E01      L6 IT
    ## F2S4_160128_013_F01        Sst
    ## F2S4_160128_013_G01         NP
    ## F2S4_160128_013_H01        Sst
    ## F2S4_160128_014_A01      L5 IT
    ## F2S4_160128_014_B01      L6 IT
    ## F2S4_160128_014_C01      L6 IT
    ## F2S4_160128_014_D01      Pvalb
    ## F2S4_160128_014_E01      L6 IT
    ## F2S4_160128_014_G01      L6 IT
    ## F2S4_160128_014_H01         NP
    ## F2S4_160128_015_A01       Sncg
    ## F2S4_160128_015_B01      L6 IT
    ## F2S4_160128_015_C01      L6 IT
    ## F2S4_160128_015_D01      L6 IT
    ## F2S4_160128_015_F01         NP
    ## F2S4_160128_015_G01      L6 IT
    ## F2S4_160128_015_H01         NP
    ## F2S4_160128_016_A01      L6 IT
    ## F2S4_160128_016_B01      L5 IT
    ## F2S4_160128_016_D01      L6 IT
    ## F2S4_160128_016_E01      L6 IT
    ## F2S4_160128_016_F01        Sst
    ## F2S4_160128_016_G01         NP
    ## F2S4_160128_016_H01      L6 IT
    ## F2S4_160128_017_A01      L6 IT
    ## F2S4_160128_017_B01      L6 IT
    ## F2S4_160128_017_C01      L5 IT
    ## F2S4_160128_017_D01      L5 IT
    ## F2S4_160128_017_E01      L6 IT
    ## F2S4_160128_017_F01      L6 IT
    ## F2S4_160128_017_G01      L6 IT
    ## F2S4_160128_017_H01      L6 IT
    ## F2S4_160128_018_A01      L6 IT
    ## F2S4_160128_018_B01      L6 CT
    ## F2S4_160128_018_C01        Sst
    ## F2S4_160128_018_D01        Sst
    ## F2S4_160128_018_E01      L6 IT
    ## F2S4_160128_018_F01      L5 IT
    ## F2S4_160128_018_G01      Lamp5
    ## F2S4_160128_018_H01         NP
    ## F2S4_160128_019_A01      L6 CT
    ## F2S4_160128_019_B01      L6 IT
    ## F2S4_160128_019_C01      L6 CT
    ## F2S4_160128_019_D01      L6 IT
    ## F2S4_160128_019_E01      L6 IT
    ## F2S4_160128_019_F01      L6 IT
    ## F2S4_160128_019_G01      L6 IT
    ## F2S4_160128_019_H01      L6 IT
    ## F2S4_160128_020_A01      L6 IT
    ## F2S4_160128_020_B01      L6 CT
    ## F2S4_160128_020_C01      Pvalb
    ## F2S4_160128_020_D01      L6 CT
    ## F2S4_160128_020_E01        Sst
    ## F2S4_160128_020_F01      L6 IT
    ## F2S4_160128_020_G01      L6 CT
    ## F2S4_160128_020_H01        Sst
    ## F2S4_160128_021_A01        Sst
    ## F2S4_160128_021_B01      L6 IT
    ## F2S4_160128_021_C01      L6 IT
    ## F2S4_160128_021_D01      L6 CT
    ## F2S4_160128_021_E01      Pvalb
    ## F2S4_160128_021_F01        Sst
    ## F2S4_160128_021_G01      L6 CT
    ## F2S4_160128_021_H01      L6 IT
    ## F2S4_160128_022_A01      L6 CT
    ## F2S4_160128_022_B01      L6 CT
    ## F2S4_160128_022_C01      L6 CT
    ## F2S4_160128_022_D01      L6 IT
    ## F2S4_160128_022_E01      L6 CT
    ## F2S4_160128_022_F01      L6 IT
    ## F2S4_160128_022_G01      L6 CT
    ## F2S4_160128_022_H01        Vip
    ## F2S4_160128_023_A01      L6 CT
    ## F2S4_160128_023_B01      L6 CT
    ## F2S4_160128_023_C01      L6 CT
    ## F2S4_160128_023_D01      L6 CT
    ## F2S4_160128_023_E01      L6 CT
    ## F2S4_160128_023_F01      L6 CT
    ## F2S4_160128_023_G01      L6 IT
    ## F2S4_160128_023_H01      L6 IT
    ## F2S4_160128_024_A01        Vip
    ## F2S4_160128_024_B01      L6 IT
    ## F2S4_160128_024_C01      L6 IT
    ## F2S4_160128_024_D01      L6 CT
    ## F2S4_160128_024_E01      L6 CT
    ## F2S4_160128_024_F01      L6 IT
    ## F2S4_160128_024_G01      L6 IT
    ## F2S4_160128_024_H01      L6 CT
    ## F2S4_160128_025_A01      L6 CT
    ## F2S4_160128_025_B01      L6 CT
    ## F2S4_160128_025_C01        Vip
    ## F2S4_160128_025_D01        Vip
    ## F2S4_160128_025_E01      L6 CT
    ## F2S4_160128_025_F01      L6 IT
    ## F2S4_160128_025_G01      L6 CT
    ## F2S4_160128_025_H01      L6 CT
    ## F2S4_160128_026_A01        Sst
    ## F2S4_160128_026_B01      L6 IT
    ## F2S4_160128_026_C01      Pvalb
    ## F2S4_160128_026_D01      Pvalb
    ## F2S4_160128_026_E01      L6 IT
    ## F2S4_160128_026_F01        Vip
    ## F2S4_160128_026_G01      L6 IT
    ## F2S4_160128_026_H01      L6 IT
    ## F2S4_160128_027_A01      L6 IT
    ## F2S4_160128_027_B01      L6 CT
    ## F2S4_160128_027_C01      L6 CT
    ## F2S4_160128_027_D01      L6 CT
    ## F2S4_160128_027_E01      L6 IT
    ## F2S4_160128_027_F01      L6 CT
    ## F2S4_160128_027_G01      L6 IT
    ## F2S4_160128_027_H01      L6 CT
    ## F2S4_160128_028_A01      L6 IT
    ## F2S4_160128_028_D01      Lamp5
    ## F2S4_160128_028_E01      L6 IT
    ## F2S4_160128_028_F01      L6 IT
    ## F2S4_160128_028_G01      L6 IT
    ## F2S4_160128_028_H01      L6 IT
    ## F2S4_160129_001_A01      L6 IT
    ## F2S4_160129_001_B01      L6 IT
    ## F2S4_160129_001_D01      L6 IT
    ## F2S4_160129_001_E01      L6 IT
    ## F2S4_160129_001_F01      L6 IT
    ## F2S4_160129_001_G01      L6 CT
    ## F2S4_160129_001_H01      L6 IT
    ## F2S4_160129_002_A01      L6 IT
    ## F2S4_160129_002_B01      L6 IT
    ## F2S4_160129_002_D01      L6 IT
    ## F2S4_160129_002_E01      L6 IT
    ## F2S4_160129_002_F01      L6 CT
    ## F2S4_160129_002_H01      L6 IT
    ## F2S4_160129_003_A01      L6 CT
    ## F2S4_160129_003_B01      L6 IT
    ## F2S4_160129_003_C01      L6 IT
    ## F2S4_160129_003_D01      L6 CT
    ## F2S4_160129_003_E01      L6 IT
    ## F2S4_160129_003_F01      L6 IT
    ## F2S4_160129_003_G01      L6 CT
    ## F2S4_160129_003_H01      L6 CT
    ## F2S4_160129_004_A01      L6 CT
    ## F2S4_160129_004_B01      L6 IT
    ## F2S4_160129_004_C01      L6 CT
    ## F2S4_160129_004_D01      L6 IT
    ## F2S4_160129_004_E01      L6 IT
    ## F2S4_160129_004_F01      L6 IT
    ## F2S4_160129_004_G01      L6 IT
    ## F2S4_160129_004_H01      L6 IT
    ## F2S4_160129_005_A01      L6 IT
    ## F2S4_160129_005_B01      L6 IT
    ## F2S4_160129_005_C01      L6 CT
    ## F2S4_160129_005_D01      L6 CT
    ## F2S4_160129_005_E01      L6 IT
    ## F2S4_160129_005_F01      L6 CT
    ## F2S4_160129_005_G01      L6 IT
    ## F2S4_160129_005_H01      L6 IT
    ## F2S4_160129_006_A01      L6 IT
    ## F2S4_160129_006_B01      L6 IT
    ## F2S4_160129_006_C01      L6 IT
    ## F2S4_160129_006_D01      L6 IT
    ## F2S4_160129_006_E01      L6 IT
    ## F2S4_160129_006_F01      L6 IT
    ## F2S4_160129_006_G01      L6 CT
    ## F2S4_160129_006_H01      L6 IT
    ## F2S4_160129_007_A01      L6 CT
    ## F2S4_160129_007_B01      L6 IT
    ## F2S4_160129_007_C01      L6 CT
    ## F2S4_160129_007_D01      L6 IT
    ## F2S4_160129_007_E01      L6 IT
    ## F2S4_160129_007_F01      L6 IT
    ## F2S4_160129_007_G01      L6 IT
    ## F2S4_160129_007_H01      L6 IT
    ## F2S4_160129_008_A01      L6 CT
    ## F2S4_160129_008_B01      L6 IT
    ## F2S4_160129_008_C01      L6 CT
    ## F2S4_160129_008_D01      L6 CT
    ## F2S4_160129_008_E01      L6 IT
    ## F2S4_160129_008_F01      L6 CT
    ## F2S4_160129_008_G01      L6 IT
    ## F2S4_160129_008_H01      L6 IT
    ## F2S4_160129_009_A01      L6 IT
    ## F2S4_160129_009_B01      L6 CT
    ## F2S4_160129_009_C01      L6 CT
    ## F2S4_160129_009_D01      L6 IT
    ## F2S4_160129_009_E01      L6 IT
    ## F2S4_160129_009_F01      L6 CT
    ## F2S4_160129_009_G01      L6 IT
    ## F2S4_160129_009_H01      L6 IT
    ## F2S4_160129_010_A01      L6 IT
    ## F2S4_160129_010_B01      L6 CT
    ## F2S4_160129_010_C01      L6 IT
    ## F2S4_160129_010_D01      L6 CT
    ## F2S4_160129_010_E01      L6 IT
    ## F2S4_160129_010_F01      L6 CT
    ## F2S4_160129_010_G01      L6 CT
    ## F2S4_160129_010_H01      L6 IT
    ## F2S4_160129_011_A01      L6 IT
    ## F2S4_160129_011_B01        Sst
    ## F2S4_160129_011_C01      L6 IT
    ## F2S4_160129_011_D01      L6 IT
    ## F2S4_160129_011_E01      L6 IT
    ## F2S4_160129_011_F01      L6 CT
    ## F2S4_160129_011_G01      L6 IT
    ## F2S4_160129_011_H01      L6 IT
    ## F2S4_160129_012_A01        Vip
    ## F2S4_160129_012_B01      L6 IT
    ## F2S4_160129_012_D01      L6 IT
    ## F2S4_160129_012_E01      L6 CT
    ## F2S4_160129_012_F01      L6 IT
    ## F2S4_160129_012_H01      L6 IT
    ## F2S4_160129_013_A01         NP
    ## F2S4_160129_013_B01         NP
    ## F2S4_160129_013_C01         NP
    ## F2S4_160129_013_D01      L5 IT
    ## F2S4_160129_013_E01         NP
    ## F2S4_160129_013_F01      Pvalb
    ## F2S4_160129_013_G01         NP
    ## F2S4_160129_013_H01         NP
    ## F2S4_160129_014_A01      Pvalb
    ## F2S4_160129_014_C01         NP
    ## F2S4_160129_014_D01        Sst
    ## F2S4_160129_014_E01         NP
    ## F2S4_160129_014_F01         NP
    ## F2S4_160129_014_G01         NP
    ## F2S4_160129_014_H01      L5 IT
    ## F2S4_160129_015_A01      L5 IT
    ## F2S4_160129_015_B01      Oligo
    ## F2S4_160129_015_C01      L5 IT
    ## F2S4_160129_015_D01      L5 PT
    ## F2S4_160129_015_E01      L5 IT
    ## F2S4_160129_015_F01        Sst
    ## F2S4_160129_015_H01         NP
    ## F2S4_160129_016_A01        Sst
    ## F2S4_160129_016_B01      L5 IT
    ## F2S4_160129_016_C01        Sst
    ## F2S4_160129_016_D01      Pvalb
    ## F2S4_160129_016_E01         NP
    ## F2S4_160129_016_F01         NP
    ## F2S4_160129_016_G01         NP
    ## F2S4_160129_016_H01         NP
    ## F2S4_160129_017_C01      Pvalb
    ## F2S4_160129_017_D01      Pvalb
    ## F2S4_160129_017_E01      Pvalb
    ## F2S4_160129_017_F01      L5 PT
    ## F2S4_160129_017_G01        Sst
    ## F2S4_160129_017_H01         NP
    ## F2S4_160129_018_A01      L6 IT
    ## F2S4_160129_018_B01         NP
    ## F2S4_160129_018_C01      L5 IT
    ## F2S4_160129_018_D01      Pvalb
    ## F2S4_160129_018_E01         NP
    ## F2S4_160129_018_F01      L5 IT
    ## F2S4_160129_018_G01        Vip
    ## F2S4_160129_018_H01         NP
    ## F2S4_160129_019_A01         L4
    ## F2S4_160129_019_B01         L4
    ## F2S4_160129_019_C01         L4
    ## F2S4_160129_019_D01         L4
    ## F2S4_160129_019_E01         L4
    ## F2S4_160129_019_F01         L4
    ## F2S4_160129_019_G01         L4
    ## F2S4_160129_019_H01         L4
    ## F2S4_160129_020_A01         L4
    ## F2S4_160129_020_B01         L4
    ## F2S4_160129_020_D01         L4
    ## F2S4_160129_020_E01         L4
    ## F2S4_160129_020_F01         L4
    ## F2S4_160129_020_G01         L4
    ## F2S4_160129_020_H01         L4
    ## F2S4_160129_021_A01         L4
    ## F2S4_160129_021_B01         L4
    ## F2S4_160129_021_C01         L4
    ## F2S4_160129_021_D01         L4
    ## F2S4_160129_021_E01         L4
    ## F2S4_160129_021_F01         L4
    ## F2S4_160129_021_G01         L4
    ## F2S4_160129_021_H01         L4
    ## F2S4_160129_022_A01         L4
    ## F2S4_160129_022_B01         L4
    ## F2S4_160129_022_C01         L4
    ## F2S4_160129_022_E01         L4
    ## F2S4_160129_022_F01      Lamp5
    ## F2S4_160129_022_G01         L4
    ## F2S4_160129_022_H01         L4
    ## F2S4_160129_023_A01         L4
    ## F2S4_160129_023_B01         L4
    ## F2S4_160129_023_C01         L4
    ## F2S4_160129_023_E01         L4
    ## F2S4_160129_023_F01         L4
    ## F2S4_160129_023_G01         L4
    ## F2S4_160129_023_H01         L4
    ## F2S4_160129_024_A01         L4
    ## F2S4_160129_024_B01         L4
    ## F2S4_160129_024_C01         L4
    ## F2S4_160129_024_D01         L4
    ## F2S4_160129_024_E01         L4
    ## F2S4_160129_024_F01         L4
    ## F2S4_160129_024_G01         L4
    ## F2S4_160129_024_H01         L4
    ## F2S4_160129_025_A01    L2/3 IT
    ## F2S4_160129_025_B01         L4
    ## F2S4_160129_025_C01      Lamp5
    ## F2S4_160129_025_D01    L2/3 IT
    ## F2S4_160129_025_E01    L2/3 IT
    ## F2S4_160129_025_F01    L2/3 IT
    ## F2S4_160129_025_G01    L2/3 IT
    ## F2S4_160129_025_H01        Vip
    ## F2S4_160129_026_A01    L2/3 IT
    ## F2S4_160129_026_B01    L2/3 IT
    ## F2S4_160129_026_C01    L2/3 IT
    ## F2S4_160129_026_D01       Sncg
    ## F2S4_160129_026_E01    L2/3 IT
    ## F2S4_160129_026_F01    L2/3 IT
    ## F2S4_160129_026_G01    L2/3 IT
    ## F2S4_160129_026_H01         L4
    ## F2S4_160129_027_A01        Vip
    ## F2S4_160129_027_B01      Lamp5
    ## F2S4_160129_027_C01        Vip
    ## F2S4_160129_027_D01      Lamp5
    ## F2S4_160129_027_E01    L2/3 IT
    ## F2S4_160129_027_F01    L2/3 IT
    ## F2S4_160129_027_G01    L2/3 IT
    ## F2S4_160129_027_H01        Vip
    ## F2S4_160129_028_A01    L2/3 IT
    ## F2S4_160129_028_B01      Lamp5
    ## F2S4_160129_028_C01        Vip
    ## F2S4_160129_028_D01    L2/3 IT
    ## F2S4_160129_028_E01        Vip
    ## F2S4_160129_028_F01    L2/3 IT
    ## F2S4_160129_028_G01    L2/3 IT
    ## F2S4_160129_028_H01         L4
    ## F2S4_160129_029_A01      Lamp5
    ## F2S4_160129_029_B01    L2/3 IT
    ## F2S4_160129_029_C01    L2/3 IT
    ## F2S4_160129_029_D01      Lamp5
    ## F2S4_160129_029_E01      Lamp5
    ## F2S4_160129_029_F01    L2/3 IT
    ## F2S4_160129_029_G01    L2/3 IT
    ## F2S4_160129_029_H01    L2/3 IT
    ## F2S4_160129_030_A01    L2/3 IT
    ## F2S4_160129_030_B01    L2/3 IT
    ## F2S4_160129_030_C01    L2/3 IT
    ## F2S4_160129_030_D01        Sst
    ## F2S4_160129_030_E01      Lamp5
    ## F2S4_160129_030_F01    L2/3 IT
    ## F2S4_160129_030_G01      Lamp5
    ## F2S4_160129_030_H01    L2/3 IT
    ## F2S4_160303_013_A01        Sst
    ## F2S4_160303_013_B01      Lamp5
    ## F2S4_160303_013_C01        Sst
    ## F2S4_160303_013_D01        Sst
    ## F2S4_160303_013_E01        Sst
    ## F2S4_160303_013_F01        Sst
    ## F2S4_160303_013_G01      Pvalb
    ## F2S4_160303_013_H01        Sst
    ## F2S4_160303_014_A01      Pvalb
    ## F2S4_160303_014_B01        Sst
    ## F2S4_160303_014_C01        Sst
    ## F2S4_160303_014_D01        Sst
    ## F2S4_160303_014_E01        Sst
    ## F2S4_160303_014_F01        Vip
    ## F2S4_160303_014_G01      Pvalb
    ## F2S4_160303_014_H01        Sst
    ## F2S4_160303_015_A01      Lamp5
    ## F2S4_160303_015_C01        Vip
    ## F2S4_160303_015_D01        Sst
    ## F2S4_160303_015_F01        Vip
    ## F2S4_160303_015_G01        Sst
    ## F2S4_160303_015_H01      Pvalb
    ## F2S4_160303_016_A01      Lamp5
    ## F2S4_160303_016_B01      Lamp5
    ## F2S4_160303_016_C01        Sst
    ## F2S4_160303_016_D01        Sst
    ## F2S4_160303_016_E01       Sncg
    ## F2S4_160303_016_F01        Vip
    ## F2S4_160303_016_G01        Sst
    ## F2S4_160303_016_H01      Lamp5
    ## F2S4_160303_017_A01        Vip
    ## F2S4_160303_017_C01      Lamp5
    ## F2S4_160303_017_D01        Vip
    ## F2S4_160303_017_E01      Lamp5
    ## F2S4_160303_017_F01        Vip
    ## F2S4_160303_017_G01      Lamp5
    ## F2S4_160303_017_H01      Lamp5
    ## F2S4_160303_018_A01        Vip
    ## F2S4_160303_018_B01        Vip
    ## F2S4_160303_018_C01      Lamp5
    ## F2S4_160303_018_D01      Pvalb
    ## F2S4_160303_018_E01      Lamp5
    ## F2S4_160303_018_F01        Vip
    ## F2S4_160303_018_H01        Vip
    ## F2S4_160303_019_A01        Vip
    ## F2S4_160303_019_B01        Vip
    ## F2S4_160303_019_C01        Vip
    ## F2S4_160303_019_D01        Vip
    ## F2S4_160303_019_E01      Lamp5
    ## F2S4_160303_019_F01      Lamp5
    ## F2S4_160303_019_G01      Pvalb
    ## F2S4_160303_019_H01        Vip
    ## F2S4_160303_020_A01        Sst
    ## F2S4_160303_020_B01        Vip
    ## F2S4_160303_020_C01       Sncg
    ## F2S4_160303_020_D01        Vip
    ## F2S4_160303_020_E01        Sst
    ## F2S4_160303_020_F01      Lamp5
    ## F2S4_160303_020_G01      Lamp5
    ## F2S4_160303_020_H01        Vip
    ## F2S4_160303_021_A01        Vip
    ## F2S4_160303_021_B01       Sncg
    ## F2S4_160303_021_C01        Vip
    ## F2S4_160303_021_D01        Vip
    ## F2S4_160303_021_E01        Vip
    ## F2S4_160303_021_F01      Lamp5
    ## F2S4_160303_021_G01        Vip
    ## F2S4_160303_021_H01        Vip
    ## F2S4_160303_022_A01        Vip
    ## F2S4_160303_022_B01        Vip
    ## F2S4_160303_022_C01        Vip
    ## F2S4_160303_022_D01        Vip
    ## F2S4_160303_022_E01        Vip
    ## F2S4_160303_022_F01        Vip
    ## F2S4_160303_022_G01       Sncg
    ## F2S4_160303_022_H01        Vip
    ## F2S4_160303_023_A01        Vip
    ## F2S4_160303_023_C01        Vip
    ## F2S4_160303_023_D01        Vip
    ## F2S4_160303_023_E01        Vip
    ## F2S4_160303_023_F01        Vip
    ## F2S4_160303_023_G01        Vip
    ## F2S4_160303_023_H01        Vip
    ## F2S4_160304_001_A01      Lamp5
    ## F2S4_160304_001_C01      Lamp5
    ## F2S4_160304_001_D01      Lamp5
    ## F2S4_160304_001_F01      Lamp5
    ## F2S4_160304_001_G01        Vip
    ## F2S4_160304_001_H01        Vip
    ## F2S4_160304_002_A01      Lamp5
    ## F2S4_160304_002_C01      Lamp5
    ## F2S4_160304_002_D01      Astro
    ## F2S4_160304_002_E01        Vip
    ## F2S4_160304_002_F01      Lamp5
    ## F2S4_160304_002_G01      Lamp5
    ## F2S4_160304_002_H01        Vip
    ## F2S4_160304_003_A01      Lamp5
    ## F2S4_160304_003_B01        Vip
    ## F2S4_160304_003_C01      Lamp5
    ## F2S4_160304_003_D01        Vip
    ## F2S4_160304_003_E01      Lamp5
    ## F2S4_160304_003_F01      Lamp5
    ## F2S4_160304_003_G01      Lamp5
    ## F2S4_160304_003_H01      Lamp5
    ## F2S4_160304_004_A01        Vip
    ## F2S4_160304_004_B01      Lamp5
    ## F2S4_160304_004_C01        Vip
    ## F2S4_160304_004_D01      Lamp5
    ## F2S4_160304_004_E01      Lamp5
    ## F2S4_160304_004_F01        Vip
    ## F2S4_160304_004_G01      Lamp5
    ## F2S4_160304_004_H01        Vip
    ## F2S4_160304_005_A01      Pvalb
    ## F2S4_160304_005_B01        Sst
    ## F2S4_160304_005_C01        Sst
    ## F2S4_160304_005_D01        Sst
    ## F2S4_160304_005_E01      Astro
    ## F2S4_160304_005_F01        Vip
    ## F2S4_160304_005_G01        Sst
    ## F2S4_160304_005_H01       Sncg
    ## F2S4_160304_006_A01        Vip
    ## F2S4_160304_006_B01      Astro
    ## F2S4_160304_006_C01        Vip
    ## F2S4_160304_006_D01        Vip
    ## F2S4_160304_006_E01        Vip
    ## F2S4_160304_006_F01        Sst
    ## F2S4_160304_006_G01      Lamp5
    ## F2S4_160304_006_H01      Pvalb
    ## F2S4_160304_007_A01        Sst
    ## F2S4_160304_007_B01      Astro
    ## F2S4_160304_007_C01   Serpinf1
    ## F2S4_160304_007_D01        Sst
    ## F2S4_160304_007_E01      Pvalb
    ## F2S4_160304_007_F01        Vip
    ## F2S4_160304_007_G01        Vip
    ## F2S4_160304_007_H01      Astro
    ## F2S4_160304_008_A01        Sst
    ## F2S4_160304_008_B01        Vip
    ## F2S4_160304_008_C01   Serpinf1
    ## F2S4_160304_008_D01       Sncg
    ## F2S4_160304_008_E01      Pvalb
    ## F2S4_160304_008_F01        Sst
    ## F2S4_160304_008_G01      Pvalb
    ## F2S4_160304_008_H01        Sst
    ## F2S4_160304_009_A01        Sst
    ## F2S4_160304_009_B01        Vip
    ## F2S4_160304_009_C01        Sst
    ## F2S4_160304_009_D01        Vip
    ## F2S4_160304_009_F01       Sncg
    ## F2S4_160304_009_G01        Sst
    ## F2S4_160304_009_H01        Sst
    ## F2S4_160304_010_A01        Vip
    ## F2S4_160304_010_B01        Sst
    ## F2S4_160304_010_C01        Vip
    ## F2S4_160304_010_D01        Sst
    ## F2S4_160304_010_E01      Pvalb
    ## F2S4_160304_010_F01      Lamp5
    ## F2S4_160304_010_G01        Vip
    ## F2S4_160304_010_H01        Sst
    ## F2S4_160304_011_A01        Vip
    ## F2S4_160304_011_B01        Vip
    ## F2S4_160304_011_C01        Sst
    ## F2S4_160304_011_D01      Lamp5
    ## F2S4_160304_011_E01      Astro
    ## F2S4_160304_011_F01      Pvalb
    ## F2S4_160304_011_G01        Vip
    ## F2S4_160304_012_A01      Lamp5
    ## F2S4_160304_012_B01        Vip
    ## F2S4_160304_012_C01        Vip
    ## F2S4_160304_012_D01        Vip
    ## F2S4_160304_012_E01        Vip
    ## F2S4_160304_012_F01        Vip
    ## F2S4_160304_012_G01        Vip
    ## F2S4_160304_012_H01        Vip
    ## F2S4_160304_013_A01      Pvalb
    ## F2S4_160304_013_B01        Vip
    ## F2S4_160304_013_C01      Lamp5
    ## F2S4_160304_013_D01        Vip
    ## F2S4_160304_013_E01        Vip
    ## F2S4_160304_013_F01        Vip
    ## F2S4_160304_013_G01        Sst
    ## F2S4_160304_013_H01        Vip
    ## F2S4_160304_014_A01        Vip
    ## F2S4_160304_014_C01      Lamp5
    ## F2S4_160304_014_D01      Lamp5
    ## F2S4_160304_014_F01        Vip
    ## F2S4_160304_014_G01        Vip
    ## F2S4_160304_014_H01        Vip
    ## F2S4_160304_015_A01        Vip
    ## F2S4_160304_015_B01        Vip
    ## F2S4_160304_015_C01        Vip
    ## F2S4_160304_015_D01        Vip
    ## F2S4_160304_015_E01        Vip
    ## F2S4_160304_015_G01        Vip
    ## F2S4_160304_015_H01      Lamp5
    ## F2S4_160304_016_A01      Lamp5
    ## F2S4_160304_016_B01      Lamp5
    ## F2S4_160304_016_C01      Pvalb
    ## F2S4_160304_016_D01    L2/3 IT
    ## F2S4_160304_016_E01        Vip
    ## F2S4_160304_016_F01      Lamp5
    ## F2S4_160304_016_G01        Vip
    ## F2S4_160304_016_H01        Vip
    ## F2S4_160304_017_A01        Vip
    ## F2S4_160304_017_B01        Vip
    ## F2S4_160304_017_C01      Lamp5
    ## F2S4_160304_017_D01        Vip
    ## F2S4_160304_017_E01    L2/3 IT
    ## F2S4_160304_017_F01        Vip
    ## F2S4_160304_017_G01      Lamp5
    ## F2S4_160304_017_H01       Sncg
    ## F2S4_160304_018_A01      Pvalb
    ## F2S4_160304_018_B01      Astro
    ## F2S4_160304_018_C01      Lamp5
    ## F2S4_160304_018_D01       Sncg
    ## F2S4_160304_018_E01      Lamp5
    ## F2S4_160304_018_F01      Lamp5
    ## F2S4_160304_018_G01        Vip
    ## F2S4_160304_018_H01        Vip
    ## F2S4_160304_019_A01        Vip
    ## F2S4_160304_019_B01      Lamp5
    ## F2S4_160304_019_C01        Vip
    ## F2S4_160304_019_D01        Vip
    ## F2S4_160304_019_E01      Lamp5
    ## F2S4_160304_019_F01      Lamp5
    ## F2S4_160304_019_G01        Vip
    ## F2S4_160304_019_H01      Lamp5
    ## F2S4_160304_020_A01      Astro
    ## F2S4_160304_020_C01        Vip
    ## F2S4_160304_020_D01        Vip
    ## F2S4_160304_020_E01      Lamp5
    ## F2S4_160304_020_F01      Lamp5
    ## F2S4_160304_020_G01        Vip
    ## F2S4_160304_020_H01        Sst
    ## F2S4_160304_021_A01      Lamp5
    ## F2S4_160304_021_B01      Lamp5
    ## F2S4_160304_021_C01      Lamp5
    ## F2S4_160304_021_D01      Lamp5
    ## F2S4_160304_021_F01      Pvalb
    ## F2S4_160304_021_G01      Lamp5
    ## F2S4_160304_021_H01      Lamp5
    ## F2S4_160317_002_A01      L6 IT
    ## F2S4_160317_002_B01      L6 CT
    ## F2S4_160317_002_C01      L5 IT
    ## F2S4_160317_002_D01        Vip
    ## F2S4_160317_002_E01      L6 IT
    ## F2S4_160317_002_F01      L6 CT
    ## F2S4_160317_002_G01      L6 IT
    ## F2S4_160317_002_H01      L6 IT
    ## F2S4_160317_003_A01      L6 IT
    ## F2S4_160317_003_B01      L6 CT
    ## F2S4_160317_003_C01      L6 CT
    ## F2S4_160317_003_D01      L6 IT
    ## F2S4_160317_003_E01      L6 IT
    ## F2S4_160317_003_F01      L6 CT
    ## F2S4_160317_003_G01        Sst
    ## F2S4_160317_003_H01      L6 IT
    ## F2S4_160317_004_A01      L6 IT
    ## F2S4_160317_004_B01      L6 IT
    ## F2S4_160317_004_C01      L6 IT
    ## F2S4_160317_004_D01      L6 CT
    ## F2S4_160317_004_E01      L6 IT
    ## F2S4_160317_004_F01      L6 IT
    ## F2S4_160317_004_G01      L6 CT
    ## F2S4_160317_004_H01      L6 IT
    ## F2S4_160317_005_A01      L6 CT
    ## F2S4_160317_005_B01      L6 IT
    ## F2S4_160317_005_C01      L6 IT
    ## F2S4_160317_005_D01      L6 CT
    ## F2S4_160317_005_E01      L6 IT
    ## F2S4_160317_005_F01      L6 IT
    ## F2S4_160317_005_G01      L6 CT
    ## F2S4_160317_005_H01      L6 IT
    ## F2S4_160317_006_A01      L6 IT
    ## F2S4_160317_006_B01      L6 IT
    ## F2S4_160317_006_C01      L6 CT
    ## F2S4_160317_006_D01      L6 IT
    ## F2S4_160317_006_E01      L6 IT
    ## F2S4_160317_006_F01      L6 CT
    ## F2S4_160317_006_G01      L6 IT
    ## F2S4_160317_006_H01      L6 IT
    ## F2S4_160317_007_A01        Sst
    ## F2S4_160317_007_B01         NP
    ## F2S4_160317_007_C01      L5 PT
    ## F2S4_160317_007_D01         NP
    ## F2S4_160317_007_E01         NP
    ## F2S4_160317_007_F01         NP
    ## F2S4_160317_007_G01      L5 IT
    ## F2S4_160317_007_H01         NP
    ## F2S4_160317_008_A01         NP
    ## F2S4_160317_008_B01      L5 PT
    ## F2S4_160317_008_C01         NP
    ## F2S4_160317_008_D01      L5 IT
    ## F2S4_160317_008_E01         NP
    ## F2S4_160317_008_F01      L6 IT
    ## F2S4_160317_008_G01      L6 IT
    ## F2S4_160317_008_H01         NP
    ## F2S4_160317_009_A01      L5 IT
    ## F2S4_160317_009_B01      L6 IT
    ## F2S4_160317_009_C01         NP
    ## F2S4_160317_009_D01      L5 IT
    ## F2S4_160317_009_E01      L5 IT
    ## F2S4_160317_009_F01      Lamp5
    ## F2S4_160317_009_G01      L5 IT
    ## F2S4_160317_009_H01        Vip
    ## F2S4_160317_010_B01         NP
    ## F2S4_160317_010_C01      L5 IT
    ## F2S4_160317_010_D01         NP
    ## F2S4_160317_010_E01         NP
    ## F2S4_160317_010_F01         NP
    ## F2S4_160317_010_G01         NP
    ## F2S4_160317_011_A01      L5 IT
    ## F2S4_160317_011_B01      L5 IT
    ## F2S4_160317_011_C01         NP
    ## F2S4_160317_011_D01      L5 IT
    ## F2S4_160317_011_E01         NP
    ## F2S4_160317_011_F01      L5 IT
    ## F2S4_160317_011_G01        Vip
    ## F2S4_160317_011_H01         NP
    ## F2S4_160317_013_A01         L4
    ## F2S4_160317_013_B01         L4
    ## F2S4_160317_013_C01         L4
    ## F2S4_160317_013_D01         L4
    ## F2S4_160317_013_E01         L4
    ## F2S4_160317_013_F01         L4
    ## F2S4_160317_013_G01         L4
    ## F2S4_160317_013_H01         L4
    ## F2S4_160317_014_B01         L4
    ## F2S4_160317_014_C01         L4
    ## F2S4_160317_014_D01         L4
    ## F2S4_160317_014_E01         L4
    ## F2S4_160317_014_F01        Vip
    ## F2S4_160317_014_H01         L4
    ## F2S4_160317_015_A01         L4
    ## F2S4_160317_015_B01         L4
    ## F2S4_160317_015_C01        Vip
    ## F2S4_160317_015_D01         L4
    ## F2S4_160317_015_E01        Vip
    ## F2S4_160317_015_F01         L4
    ## F2S4_160317_015_G01         L4
    ## F2S4_160317_015_H01         L4
    ## F2S4_160317_016_A01         L4
    ## F2S4_160317_016_B01         L4
    ## F2S4_160317_016_C01         L4
    ## F2S4_160317_016_D01         L4
    ## F2S4_160317_016_E01         L4
    ## F2S4_160317_016_F01         L4
    ## F2S4_160317_016_G01         L4
    ## F2S4_160317_016_H01         L4
    ## F2S4_160317_017_A01    L2/3 IT
    ## F2S4_160317_017_B01        Vip
    ## F2S4_160317_017_C01    L2/3 IT
    ## F2S4_160317_017_D01    L2/3 IT
    ## F2S4_160317_017_E01    L2/3 IT
    ## F2S4_160317_017_F01    L2/3 IT
    ## F2S4_160317_017_G01    L2/3 IT
    ## F2S4_160317_017_H01      Lamp5
    ## F2S4_160317_018_A01    L2/3 IT
    ## F2S4_160317_018_B01    L2/3 IT
    ## F2S4_160317_018_C01    L2/3 IT
    ## F2S4_160317_018_D01         L4
    ## F2S4_160317_018_E01         L4
    ## F2S4_160317_018_F01    L2/3 IT
    ## F2S4_160317_018_G01    L2/3 IT
    ## F2S4_160317_018_H01    L2/3 IT
    ## F2S4_160317_019_A01    L2/3 IT
    ## F2S4_160317_019_B01    L2/3 IT
    ## F2S4_160317_019_C01         L4
    ## F2S4_160317_019_D01    L2/3 IT
    ## F2S4_160317_019_E01    L2/3 IT
    ## F2S4_160317_019_F01    L2/3 IT
    ## F2S4_160317_019_G01    L2/3 IT
    ## F2S4_160317_019_H01    L2/3 IT
    ## F2S4_160317_020_A01    L2/3 IT
    ## F2S4_160317_020_B01    L2/3 IT
    ## F2S4_160317_020_C01    L2/3 IT
    ## F2S4_160317_020_D01    L2/3 IT
    ## F2S4_160317_020_E01    L2/3 IT
    ## F2S4_160317_020_F01    L2/3 IT
    ## F2S4_160317_020_G01    L2/3 IT
    ## F2S4_160317_020_H01        Vip
    ## F2S4_160317_021_A01    L2/3 IT
    ## F2S4_160317_021_B01        Vip
    ## F2S4_160317_021_C01         L4
    ## F2S4_160317_021_D01    L2/3 IT
    ## F2S4_160317_021_E01    L2/3 IT
    ## F2S4_160317_021_F01      Lamp5
    ## F2S4_160317_021_G01    L2/3 IT
    ## F2S4_160317_021_H01    L2/3 IT
    ## F2S4_160317_022_A01      L6 IT
    ## F2S4_160317_022_B01      L6 IT
    ## F2S4_160317_022_C01        Vip
    ## F2S4_160317_022_D01      L6 IT
    ## F2S4_160317_022_E01      L6 IT
    ## F2S4_160317_022_F01      L6 IT
    ## F2S4_160317_022_G01        Sst
    ## F2S4_160317_022_H01      L6 IT
    ## F2S4_160317_023_A01         NP
    ## F2S4_160317_023_B01         NP
    ## F2S4_160317_023_C01        Vip
    ## F2S4_160317_023_D01      L6 IT
    ## F2S4_160317_023_E01      L6 IT
    ## F2S4_160317_023_F01      L5 IT
    ## F2S4_160317_023_G01        Sst
    ## F2S4_160317_023_H01      L5 IT
    ## F2S4_160317_024_A01      L6 IT
    ## F2S4_160317_024_B01      L5 IT
    ## F2S4_160317_024_C01      L6 IT
    ## F2S4_160317_024_D01      L6 IT
    ## F2S4_160317_024_E01       Sncg
    ## F2S4_160317_024_F01        Sst
    ## F2S4_160317_024_G01      L6 CT
    ## F2S4_160317_024_H01        Sst
    ## F2S4_160317_025_A01      L6 IT
    ## F2S4_160317_025_B01        Sst
    ## F2S4_160317_025_C01         NP
    ## F2S4_160317_025_D01         NP
    ## F2S4_160317_025_E01      L6 IT
    ## F2S4_160317_025_F01      L5 IT
    ## F2S4_160317_025_G01      L5 IT
    ## F2S4_160317_025_H01      L6 IT
    ## F2S4_160317_026_A01      L6 IT
    ## F2S4_160317_026_B01        Vip
    ## F2S4_160317_026_D01      L5 PT
    ## F2S4_160317_026_E01      L6 IT
    ## F2S4_160317_026_F01      L6 IT
    ## F2S4_160317_026_G01      L6 IT
    ## F2S4_160317_026_H01        Vip
    ## F2S4_160317_027_A01      L6 CT
    ## F2S4_160317_027_B01      L6 CT
    ## F2S4_160317_027_C01      L6 CT
    ## F2S4_160317_027_D01      L6 IT
    ## F2S4_160317_027_E01      L6 IT
    ## F2S4_160317_027_F01      L6 CT
    ## F2S4_160317_027_G01      L6 CT
    ## F2S4_160317_027_H01      L6 CT
    ## F2S4_160317_028_A01      L6 IT
    ## F2S4_160317_028_B01      L6 IT
    ## F2S4_160317_028_C01      L6 IT
    ## F2S4_160317_028_D01      L6 IT
    ## F2S4_160317_028_F01      L6 CT
    ## F2S4_160317_028_G01      L6 CT
    ## F2S4_160317_028_H01      L6 CT
    ## F2S4_160317_029_A01      L6 CT
    ## F2S4_160317_029_B01      L6 CT
    ## F2S4_160317_029_C01      L6 CT
    ## F2S4_160317_029_D01      L6 IT
    ## F2S4_160317_029_E01      L6 IT
    ## F2S4_160317_029_F01      L6 CT
    ## F2S4_160317_029_G01      L6 CT
    ## F2S4_160317_029_H01      L6 IT
    ## F2S4_160317_030_A01      L6 IT
    ## F2S4_160317_030_B01      L6 IT
    ## F2S4_160317_030_C01      L6 IT
    ## F2S4_160317_030_D01      L6 IT
    ## F2S4_160317_030_E01      L6 IT
    ## F2S4_160317_030_F01      L6 IT
    ## F2S4_160317_030_G01      L6 IT
    ## F2S4_160317_030_H01      L6 IT
    ## F2S4_160317_031_A01      L6 IT
    ## F2S4_160317_031_B01      L6 IT
    ## F2S4_160317_031_C01      L6 CT
    ## F2S4_160317_031_D01      L6 CT
    ## F2S4_160317_031_E01      L6 CT
    ## F2S4_160317_031_F01      L6 CT
    ## F2S4_160317_031_G01      L6 CT
    ## F2S4_160317_031_H01      L6 IT
    ## F2S4_160321_001_A01      Lamp5
    ## F2S4_160321_001_B01        Vip
    ## F2S4_160321_001_C01         L4
    ## F2S4_160321_001_D01        Vip
    ## F2S4_160321_001_E01      Lamp5
    ## F2S4_160321_001_F01      Lamp5
    ## F2S4_160321_001_G01    L2/3 IT
    ## F2S4_160321_001_H01        Vip
    ## F2S4_160321_002_A01        Vip
    ## F2S4_160321_002_B01    L2/3 IT
    ## F2S4_160321_002_C01    L2/3 IT
    ## F2S4_160321_002_D01        Vip
    ## F2S4_160321_002_E01        Vip
    ## F2S4_160321_002_F01      Lamp5
    ## F2S4_160321_002_G01    L2/3 IT
    ## F2S4_160321_002_H01      Lamp5
    ## F2S4_160321_003_B01         L4
    ## F2S4_160321_003_C01        Vip
    ## F2S4_160321_003_D01    L2/3 IT
    ## F2S4_160321_003_E01    L2/3 IT
    ## F2S4_160321_003_F01    L2/3 IT
    ## F2S4_160321_003_G01    L2/3 IT
    ## F2S4_160321_003_H01      Lamp5
    ## F2S4_160321_006_A01         L4
    ## F2S4_160321_006_C01         L4
    ## F2S4_160321_006_D01         L4
    ## F2S4_160321_006_E01         L4
    ## F2S4_160321_006_F01      Lamp5
    ## F2S4_160321_006_G01         L4
    ## F2S4_160321_006_H01      L5 IT
    ## F2S4_160321_007_A01      Pvalb
    ## F2S4_160321_007_B01         L4
    ## F2S4_160321_007_C01         L4
    ## F2S4_160321_007_D01         L4
    ## F2S4_160321_007_F01        Vip
    ## F2S4_160321_007_G01         L4
    ## F2S4_160321_007_H01         L4
    ## F2S4_160321_008_A01         L4
    ## F2S4_160321_008_B01      L5 IT
    ## F2S4_160321_008_C01         L4
    ## F2S4_160321_008_E01      L5 IT
    ## F2S4_160321_008_F01         NP
    ## F2S4_160321_011_A01         NP
    ## F2S4_160321_011_B01      L6 IT
    ## F2S4_160321_011_C01         L4
    ## F2S4_160321_011_D01      L5 IT
    ## F2S4_160321_011_E01      L6 IT
    ## F2S4_160321_011_F01      L5 IT
    ## F2S4_160321_011_G01        Sst
    ## F2S4_160321_011_H01         NP
    ## F2S4_160321_012_A01      L5 PT
    ## F2S4_160321_012_B01      L6 IT
    ## F2S4_160321_012_C01      L5 IT
    ## F2S4_160321_012_D01      Pvalb
    ## F2S4_160321_012_E01      L5 PT
    ## F2S4_160321_012_F01        Sst
    ## F2S4_160321_012_G01         NP
    ## F2S4_160321_012_H01        Sst
    ## F2S4_160321_016_A01      L6 IT
    ## F2S4_160321_016_B01      L6 IT
    ## F2S4_160321_016_C01      L6 IT
    ## F2S4_160321_016_D01      L6 IT
    ## F2S4_160321_016_E01      L6 CT
    ## F2S4_160321_016_F01      L6 CT
    ## F2S4_160321_016_G01      L6 CT
    ## F2S4_160321_016_H01      L6 IT
    ## F2S4_160321_017_A01      L6 IT
    ## F2S4_160321_017_B01      L6 IT
    ## F2S4_160321_017_C01      L6 CT
    ## F2S4_160321_017_D01      L6 CT
    ## F2S4_160321_017_E01      L6 IT
    ## F2S4_160321_017_F01      L6 IT
    ## F2S4_160321_017_G01      L6 IT
    ## F2S4_160321_017_H01      L6 IT
    ## F2S4_160321_022_A01      L5 PT
    ## F2S4_160321_022_B01      L6 IT
    ## F2S4_160321_022_C01      L6 IT
    ## F2S4_160321_022_D01      L5 IT
    ## F2S4_160321_022_E01      L6 IT
    ## F2S4_160321_022_F01      Lamp5
    ## F2S4_160321_022_G01      L5 IT
    ## F2S4_160321_022_H01      L5 IT
    ## F2S4_160321_023_A01      L5 IT
    ## F2S4_160321_023_B01      L5 IT
    ## F2S4_160321_023_C01      L5 IT
    ## F2S4_160321_023_D01      L5 PT
    ## F2S4_160321_023_E01      L5 IT
    ## F2S4_160321_023_F01      L6 IT
    ## F2S4_160321_023_G01      L6 IT
    ## F2S4_160321_023_H01      L6 IT
    ## F2S4_160321_024_A01      L5 PT
    ## F2S4_160321_024_B01        Vip
    ## F2S4_160321_024_C01      L6 IT
    ## F2S4_160321_024_D01        Sst
    ## F2S4_160321_024_E01         NP
    ## F2S4_160321_024_F01      L5 IT
    ## F2S4_160321_024_G01      L5 IT
    ## F2S4_160321_024_H01      L5 IT
    ## F2S4_160321_025_A01        Sst
    ## F2S4_160321_025_B01         NP
    ## F2S4_160321_025_C01      L5 IT
    ## F2S4_160321_025_D01        Sst
    ## F2S4_160321_025_E01         NP
    ## F2S4_160321_025_F01      Lamp5
    ## F2S4_160321_025_G01         NP
    ## F2S4_160321_025_H01      L5 PT
    ## F2S4_160321_026_A01      L6 CT
    ## F2S4_160321_026_B01      L6 CT
    ## F2S4_160321_026_C01      L6 IT
    ## F2S4_160321_026_D01        Vip
    ## F2S4_160321_026_E01      L6 IT
    ## F2S4_160321_026_F01      L6 IT
    ## F2S4_160321_026_H01      L6 IT
    ## F2S4_160321_027_A01      L6 IT
    ## F2S4_160321_027_B01      L6 IT
    ## F2S4_160321_027_C01      L6 CT
    ## F2S4_160321_027_D01      L6 IT
    ## F2S4_160321_027_E01      L6 IT
    ## F2S4_160321_027_F01      L6 IT
    ## F2S4_160321_027_G01      L6 IT
    ## F2S4_160321_027_H01      L6 IT
    ## F2S4_160321_028_A01      L6 IT
    ## F2S4_160321_028_B01      L6 CT
    ## F2S4_160321_028_C01      Lamp5
    ## F2S4_160321_028_D01      L6 IT
    ## F2S4_160321_028_E01      L6 CT
    ## F2S4_160321_028_F01      L6 CT
    ## F2S4_160321_028_G01      L6 IT
    ## F2S4_160321_028_H01      L6 IT
    ## F2S4_160321_029_A01      L6 IT
    ## F2S4_160321_029_B01      L6 IT
    ## F2S4_160321_029_C01      L6 IT
    ## F2S4_160321_029_D01      L6 IT
    ## F2S4_160321_029_E01      L6 IT
    ## F2S4_160321_029_F01      L6 IT
    ## F2S4_160321_029_G01      L6 IT
    ## F2S4_160321_029_H01      L6 IT
    ## F2S4_160321_030_A01      L6 IT
    ## F2S4_160321_030_B01      L6 IT
    ## F2S4_160321_030_C01      L6 CT
    ## F2S4_160321_030_D01      L6 CT
    ## F2S4_160321_030_E01      L6 IT
    ## F2S4_160321_030_F01      L6 IT
    ## F2S4_160321_030_G01       Sncg
    ## F2S4_160321_030_H01        Sst
    ## F2S4_160328_005_A01      L5 PT
    ## F2S4_160328_005_B01      L5 PT
    ## F2S4_160328_005_C01      L5 IT
    ## F2S4_160328_005_D01      L5 IT
    ## F2S4_160328_005_E01      L5 PT
    ## F2S4_160328_005_F01      L5 IT
    ## F2S4_160328_005_G01      L5 IT
    ## F2S4_160328_005_H01      L5 IT
    ## F2S4_160328_006_A01      L5 IT
    ## F2S4_160328_006_B01      L5 PT
    ## F2S4_160328_006_C01      L5 IT
    ## F2S4_160328_006_D01      L5 PT
    ## F2S4_160328_006_E01      L5 IT
    ## F2S4_160328_006_F01      L5 PT
    ## F2S4_160328_006_G01      L5 IT
    ## F2S4_160328_006_H01      L5 PT
    ## F2S4_160328_007_A01      L5 IT
    ## F2S4_160328_007_B01      L5 IT
    ## F2S4_160328_007_D01      L5 IT
    ## F2S4_160328_007_E01      L5 IT
    ## F2S4_160328_007_F01      L5 PT
    ## F2S4_160328_007_G01      L5 IT
    ## F2S4_160328_008_A01      L5 PT
    ## F2S4_160328_008_B01      L5 PT
    ## F2S4_160328_008_C01      L5 IT
    ## F2S4_160328_008_D01      L5 PT
    ## F2S4_160328_008_E01      L5 IT
    ## F2S4_160328_008_F01      L5 IT
    ## F2S4_160328_009_A01      L6 IT
    ## F2S4_160328_009_B01        Sst
    ## F2S4_160328_009_C01      L6 IT
    ## F2S4_160328_009_D01      L6 CT
    ## F2S4_160328_009_E01      L6 IT
    ## F2S4_160328_009_F01      L6 CT
    ## F2S4_160328_009_G01      L6 IT
    ## F2S4_160328_009_H01      L6 IT
    ## F2S4_160328_010_D01      L6 CT
    ## F2S4_160328_010_E01        Sst
    ## F2S4_160328_010_F01      L6 IT
    ## F2S4_160328_010_G01      L6 IT
    ## F2S4_160328_010_H01      L6 IT
    ## F2S4_160328_011_A01      L6 CT
    ## F2S4_160328_011_B01      L5 IT
    ## F2S4_160328_011_C01      L5 IT
    ## F2S4_160328_011_D01      L5 IT
    ## F2S4_160328_011_E01      L5 IT
    ## F2S4_160328_011_F01      L5 PT
    ## F2S4_160328_011_G01      L5 IT
    ## F2S4_160328_011_H01      L5 IT
    ## F2S4_160328_012_A01      L5 PT
    ## F2S4_160328_012_D01        Sst
    ## F2S4_160328_012_F01      L5 PT
    ## F2S4_160328_012_G01      L5 IT
    ## F2S4_160328_013_A01      L5 IT
    ## F2S4_160328_013_B01      L5 IT
    ## F2S4_160328_013_D01      L5 PT
    ## F2S4_160328_013_E01      L5 PT
    ## F2S4_160328_013_F01      L5 IT
    ## F2S4_160328_013_G01      L5 IT
    ## F2S4_160328_013_H01      L5 PT
    ## F2S4_160328_014_A01      L5 PT
    ## F2S4_160328_014_B01      L5 IT
    ## F2S4_160328_014_C01      L5 IT
    ## F2S4_160328_014_D01      L6 IT
    ## F2S4_160329_001_A01    L2/3 IT
    ## F2S4_160329_001_B01        Vip
    ## F2S4_160329_001_C01    L2/3 IT
    ## F2S4_160329_001_D01    L2/3 IT
    ## F2S4_160329_001_E01    L2/3 IT
    ## F2S4_160329_001_F01    L2/3 IT
    ## F2S4_160329_001_G01    L2/3 IT
    ## F2S4_160329_001_H01        Vip
    ## F2S4_160329_002_A01    L2/3 IT
    ## F2S4_160329_002_B01    L2/3 IT
    ## F2S4_160329_002_C01    L2/3 IT
    ## F2S4_160329_002_D01    L2/3 IT
    ## F2S4_160329_002_E01    L2/3 IT
    ## F2S4_160329_002_F01        Vip
    ## F2S4_160329_002_G01    L2/3 IT
    ## F2S4_160329_002_H01    L2/3 IT
    ## F2S4_160329_003_A01    L2/3 IT
    ## F2S4_160329_003_B01    L2/3 IT
    ## F2S4_160329_003_C01    L2/3 IT
    ## F2S4_160329_003_D01    L2/3 IT
    ## F2S4_160329_003_E01    L2/3 IT
    ## F2S4_160329_003_F01    L2/3 IT
    ## F2S4_160329_003_G01    L2/3 IT
    ## F2S4_160329_003_H01    L2/3 IT
    ## F2S4_160329_004_A01    L2/3 IT
    ## F2S4_160329_004_B01    L2/3 IT
    ## F2S4_160329_004_C01    L2/3 IT
    ## F2S4_160329_004_D01    L2/3 IT
    ## F2S4_160329_004_E01    L2/3 IT
    ## F2S4_160329_004_F01    L2/3 IT
    ## F2S4_160329_004_G01    L2/3 IT
    ## F2S4_160329_004_H01    L2/3 IT
    ## F2S4_160329_005_A01    L2/3 IT
    ## F2S4_160329_005_B01    L2/3 IT
    ## F2S4_160329_005_C01    L2/3 IT
    ## F2S4_160329_005_D01    L2/3 IT
    ## F2S4_160329_005_E01    L2/3 IT
    ## F2S4_160329_005_F01    L2/3 IT
    ## F2S4_160329_005_G01    L2/3 IT
    ## F2S4_160329_005_H01    L2/3 IT
    ## F2S4_160329_006_A01         L4
    ## F2S4_160329_006_B01         L4
    ## F2S4_160329_006_C01         L4
    ## F2S4_160329_006_D01      L5 IT
    ## F2S4_160329_006_E01         L4
    ## F2S4_160329_006_F01         L4
    ## F2S4_160329_006_G01         L4
    ## F2S4_160329_006_H01         L4
    ## F2S4_160329_007_A01         L4
    ## F2S4_160329_007_B01         L4
    ## F2S4_160329_007_C01         L4
    ## F2S4_160329_007_D01         L4
    ## F2S4_160329_007_E01         L4
    ## F2S4_160329_007_F01         L4
    ## F2S4_160329_007_G01         L4
    ## F2S4_160329_007_H01         L4
    ## F2S4_160329_008_A01         L4
    ## F2S4_160329_008_B01         L4
    ## F2S4_160329_008_C01         L4
    ## F2S4_160329_008_D01         L4
    ## F2S4_160329_008_E01         L4
    ## F2S4_160329_008_F01        Vip
    ## F2S4_160329_008_G01        Vip
    ## F2S4_160329_009_A01         L4
    ## F2S4_160329_009_B01         L4
    ## F2S4_160329_009_C01         L4
    ## F2S4_160329_009_E01        Sst
    ## F2S4_160329_009_F01         L4
    ## F2S4_160329_009_G01      L5 IT
    ## F2S4_160329_009_H01         L4
    ## F2S4_160329_010_A01         L4
    ## F2S4_160329_010_B01         L4
    ## F2S4_160329_010_C01         L4
    ## F2S4_160329_010_D01         L4
    ## F2S4_160329_010_E01         L4
    ## F2S4_160329_010_F01      L5 IT
    ## F2S4_160329_010_G01         L4
    ## F2S4_160329_010_H01         L4
    ## F2S4_160329_011_A01        Sst
    ## F2S4_160329_011_B01      L6 IT
    ## F2S4_160329_011_C01        Vip
    ## F2S4_160329_011_D01         NP
    ## F2S4_160329_011_E01         NP
    ## F2S4_160329_011_F01      L6 IT
    ## F2S4_160329_011_G01         NP
    ## F2S4_160329_011_H01      L5 PT
    ## F2S4_160329_012_A01         NP
    ## F2S4_160329_012_C01         NP
    ## F2S4_160329_012_D01         NP
    ## F2S4_160329_012_E01      L5 IT
    ## F2S4_160329_012_G01      L6 IT
    ## F2S4_160329_012_H01      L5 IT
    ## F2S4_160329_013_A01         NP
    ## F2S4_160329_013_B01      L6 IT
    ## F2S4_160329_013_C01      L5 IT
    ## F2S4_160329_013_D01        Vip
    ## F2S4_160329_013_E01         NP
    ## F2S4_160329_013_F01         NP
    ## F2S4_160329_013_G01         NP
    ## F2S4_160329_013_H01        Sst
    ## F2S4_160329_014_A01         NP
    ## F2S4_160329_014_B01         NP
    ## F2S4_160329_014_C01        Vip
    ## F2S4_160329_014_D01      L5 IT
    ## F2S4_160329_014_E01         NP
    ## F2S4_160329_014_F01         NP
    ## F2S4_160329_014_G01        Vip
    ## F2S4_160329_014_H01         NP
    ## F2S4_160329_015_A01      L6 IT
    ## F2S4_160329_015_B01         L4
    ## F2S4_160329_015_C01         NP
    ## F2S4_160329_015_D01        Vip
    ## F2S4_160329_015_E01         NP
    ## F2S4_160329_015_F01         NP
    ## F2S4_160329_015_G01      Oligo
    ## F2S4_160329_015_H01         NP
    ## F2S4_160329_016_A01        Vip
    ## F2S4_160329_016_B01        Vip
    ## F2S4_160329_016_C01        Sst
    ## F2S4_160329_016_D01      L6 IT
    ## F2S4_160329_016_E01      L6 IT
    ## F2S4_160329_016_F01      L6 CT
    ## F2S4_160329_016_G01      L6 IT
    ## F2S4_160329_016_H01      L6 IT
    ## F2S4_160329_017_A01      L6 IT
    ## F2S4_160329_017_C01      L6 IT
    ## F2S4_160329_017_D01      L6 CT
    ## F2S4_160329_017_E01      L6 IT
    ## F2S4_160329_017_F01      L6 IT
    ## F2S4_160329_017_G01      L6 IT
    ## F2S4_160329_018_A01      L6 IT
    ## F2S4_160329_018_C01        Vip
    ## F2S4_160329_018_D01      L6 CT
    ## F2S4_160329_018_E01      L6 IT
    ## F2S4_160329_018_F01      L6 CT
    ## F2S4_160329_018_G01      L6 IT
    ## F2S4_160329_018_H01      L6 IT
    ## F2S4_160329_019_A01        Vip
    ## F2S4_160329_019_B01      L6 IT
    ## F2S4_160329_019_C01      L6 IT
    ## F2S4_160329_019_D01      L6 IT
    ## F2S4_160329_019_E01      L6 CT
    ## F2S4_160329_019_F01      L6 CT
    ## F2S4_160329_019_G01      L6 CT
    ## F2S4_160329_019_H01      L6 CT
    ## F2S4_160329_020_A01      L6 IT
    ## F2S4_160329_020_B01      L6 CT
    ## F2S4_160329_020_C01      L6 IT
    ## F2S4_160329_020_D01      L6 IT
    ## F2S4_160329_020_E01      L6 IT
    ## F2S4_160329_020_F01      L6 IT
    ## F2S4_160329_020_G01      L6 CT
    ## F2S4_160329_020_H01        Sst
    ## F2S4_160329_021_A01      L6 IT
    ## F2S4_160329_021_B01      L6 CT
    ## F2S4_160329_021_C01      L6 IT
    ## F2S4_160329_021_D01      L6 IT
    ## F2S4_160329_021_E01      L6 IT
    ## F2S4_160329_021_F01      L6 CT
    ## F2S4_160329_021_G01      L6 IT
    ## F2S4_160329_021_H01      L6 IT
    ## F2S4_160329_022_A01      L6 IT
    ## F2S4_160329_022_B01      L6 IT
    ## F2S4_160329_022_C01      Meis2
    ## F2S4_160329_022_D01      L6 IT
    ## F2S4_160329_022_E01      L6 IT
    ## F2S4_160329_022_F01      L6 IT
    ## F2S4_160329_022_G01      L6 IT
    ## F2S4_160329_023_A01      L6 IT
    ## F2S4_160329_023_B01      L6 IT
    ## F2S4_160329_023_C01      L6 IT
    ## F2S4_160329_023_D01      L6 IT
    ## F2S4_160329_023_E01      L6 IT
    ## F2S4_160329_023_F01      L6 IT
    ## F2S4_160329_023_G01      L6 CT
    ## F2S4_160329_023_H01      L6 IT
    ## F2S4_160329_024_A01      L6 IT
    ## F2S4_160329_024_B01      L6 IT
    ## F2S4_160329_024_C01      L6 CT
    ## F2S4_160329_024_D01      L6 IT
    ## F2S4_160329_024_E01      L6 IT
    ## F2S4_160329_024_F01      L6 IT
    ## F2S4_160329_024_G01        Sst
    ## F2S4_160329_024_H01      L6 CT
    ## F2S4_160329_025_A01      L6 IT
    ## F2S4_160329_025_B01      L6 IT
    ## F2S4_160329_025_C01      L6 IT
    ## F2S4_160329_025_E01      L6 IT
    ## F2S4_160329_025_F01      L6 IT
    ## F2S4_160329_025_G01      L6 CT
    ## F2S4_160329_025_H01      L6 IT
    ## F2S4_160329_026_B01      L5 IT
    ## F2S4_160329_026_C01         NP
    ## F2S4_160329_026_D01         NP
    ## F2S4_160329_026_E01         NP
    ## F2S4_160329_026_F01        Vip
    ## F2S4_160329_026_G01      L5 IT
    ## F2S4_160329_026_H01      L6 IT
    ## F2S4_160329_027_A01      L5 IT
    ## F2S4_160329_027_B01      L5 IT
    ## F2S4_160329_027_C01      L6 IT
    ## F2S4_160329_027_D01      L5 IT
    ## F2S4_160329_027_E01      L6 IT
    ## F2S4_160329_027_F01      L5 IT
    ## F2S4_160329_027_G01      L6 IT
    ## F2S4_160329_027_H01      L5 IT
    ## F2S4_160329_028_A01        Sst
    ## F2S4_160329_028_B01         NP
    ## F2S4_160329_028_C01      L5 IT
    ## F2S4_160329_028_D01        Vip
    ## F2S4_160329_028_E01      L6 IT
    ## F2S4_160329_028_F01        Vip
    ## F2S4_160329_028_G01      L5 IT
    ## F2S4_160329_028_H01      L5 IT
    ## F2S4_160329_029_A01      L6 IT
    ## F2S4_160329_029_B01      L5 IT
    ## F2S4_160329_029_C01      L6 IT
    ## F2S4_160329_029_D01      L6 IT
    ## F2S4_160329_029_E01         NP
    ## F2S4_160329_029_F01      L6 IT
    ## F2S4_160329_029_G01      L6 IT
    ## F2S4_160329_029_H01      L5 PT
    ## F2S4_160329_030_B01      L6 IT
    ## F2S4_160329_030_C01        Vip
    ## F2S4_160329_030_D01      L6 IT
    ## F2S4_160329_030_E01      L6 IT
    ## F2S4_160329_030_F01         NP
    ## F2S4_160329_030_G01      L5 IT
    ## F2S4_160329_030_H01      L5 PT
    ## F2S4_160330_001_A01      Lamp5
    ## F2S4_160330_001_B01      Lamp5
    ## F2S4_160330_001_C01        Vip
    ## F2S4_160330_001_D01        Vip
    ## F2S4_160330_001_E01        Vip
    ## F2S4_160330_001_F01        Vip
    ## F2S4_160330_001_G01        Vip
    ## F2S4_160330_001_H01        Vip
    ## F2S4_160330_002_A01    L2/3 IT
    ## F2S4_160330_002_B01    L2/3 IT
    ## F2S4_160330_002_C01      Lamp5
    ## F2S4_160330_002_D01      Lamp5
    ## F2S4_160330_002_E01        Sst
    ## F2S4_160330_002_F01      Lamp5
    ## F2S4_160330_002_G01        Vip
    ## F2S4_160330_002_H01      Lamp5
    ## F2S4_160330_003_A01         L4
    ## F2S4_160330_003_B01    L2/3 IT
    ## F2S4_160330_003_C01    L2/3 IT
    ## F2S4_160330_003_E01      Lamp5
    ## F2S4_160330_003_F01      Lamp5
    ## F2S4_160330_003_G01      Lamp5
    ## F2S4_160330_003_H01      Lamp5
    ## F2S4_160330_004_A01        Vip
    ## F2S4_160330_004_B01        Vip
    ## F2S4_160330_004_C01        Vip
    ## F2S4_160330_004_D01    L2/3 IT
    ## F2S4_160330_004_E01         L4
    ## F2S4_160330_004_F01        Vip
    ## F2S4_160330_004_G01      Lamp5
    ## F2S4_160330_004_H01         L4
    ## F2S4_160330_005_A01        Vip
    ## F2S4_160330_005_B01        Vip
    ## F2S4_160330_005_C01      Lamp5
    ## F2S4_160330_005_D01        Vip
    ## F2S4_160330_005_E01      Lamp5
    ## F2S4_160330_005_F01        Vip
    ## F2S4_160330_005_G01        Vip
    ## F2S4_160330_005_H01        Vip
    ## F2S4_160330_006_A01      Pvalb
    ## F2S4_160330_006_B01      L6 CT
    ## F2S4_160330_006_C01      L6 CT
    ## F2S4_160330_006_D01      L6 IT
    ## F2S4_160330_006_E01      L6 IT
    ## F2S4_160330_006_F01        Sst
    ## F2S4_160330_006_G01      L6 CT
    ## F2S4_160330_006_H01      L6 CT
    ## F2S4_160330_007_A01      L6 CT
    ## F2S4_160330_007_B01      L6 IT
    ## F2S4_160330_007_C01      L6 CT
    ## F2S4_160330_007_D01      L6 IT
    ## F2S4_160330_007_E01      L6 CT
    ## F2S4_160330_007_F01      L6 CT
    ## F2S4_160330_007_G01         CR
    ## F2S4_160330_007_H01      L6 IT
    ## F2S4_160330_008_A01      L6 CT
    ## F2S4_160330_008_B01      L6 CT
    ## F2S4_160330_008_C01      L6 CT
    ## F2S4_160330_008_D01      L6 CT
    ## F2S4_160330_008_E01        Sst
    ## F2S4_160330_008_F01      L6 CT
    ## F2S4_160330_008_G01      L6 IT
    ## F2S4_160330_008_H01        Vip
    ## F2S4_160330_009_A01      L6 CT
    ## F2S4_160330_009_B01      L6 CT
    ## F2S4_160330_009_C01        Sst
    ## F2S4_160330_009_D01      L6 CT
    ## F2S4_160330_009_E01      L6 CT
    ## F2S4_160330_009_F01      L6 CT
    ## F2S4_160330_009_G01      L6 CT
    ## F2S4_160330_009_H01      L6 CT
    ## F2S4_160330_010_A01      L6 IT
    ## F2S4_160330_010_B01      L6 IT
    ## F2S4_160330_010_C01      L6 CT
    ## F2S4_160330_010_E01      L6 IT
    ## F2S4_160330_010_F01      L6 CT
    ## F2S4_160330_010_G01      L6 CT
    ## F2S4_160330_010_H01      L6 CT
    ## F2S4_160330_011_A01         L4
    ## F2S4_160330_011_B01         L4
    ## F2S4_160330_011_C01         L4
    ## F2S4_160330_011_D01         L4
    ## F2S4_160330_011_F01         L4
    ## F2S4_160330_011_G01         L4
    ## F2S4_160330_011_H01         L4
    ## F2S4_160330_012_A01         L4
    ## F2S4_160330_012_B01         L4
    ## F2S4_160330_012_C01         L4
    ## F2S4_160330_012_D01         L4
    ## F2S4_160330_012_E01      Lamp5
    ## F2S4_160330_012_F01         L4
    ## F2S4_160330_012_G01         L4
    ## F2S4_160330_012_H01         L4
    ## F2S4_160330_013_A01         L4
    ## F2S4_160330_013_B01         L4
    ## F2S4_160330_013_C01         L4
    ## F2S4_160330_013_D01        Vip
    ## F2S4_160330_013_E01         L4
    ## F2S4_160330_013_F01         L4
    ## F2S4_160330_013_G01         L4
    ## F2S4_160330_013_H01         L4
    ## F2S4_160330_014_A01         L4
    ## F2S4_160330_014_B01         L4
    ## F2S4_160330_014_C01         L4
    ## F2S4_160330_014_D01        Vip
    ## F2S4_160330_014_E01         L4
    ## F2S4_160330_014_F01        Vip
    ## F2S4_160330_014_G01         L4
    ## F2S4_160330_014_H01         L4
    ## F2S4_160330_015_A01         L4
    ## F2S4_160330_015_B01         L4
    ## F2S4_160330_015_D01         L4
    ## F2S4_160330_015_E01         L4
    ## F2S4_160330_015_F01        Vip
    ## F2S4_160330_015_G01         L4
    ## F2S4_160330_015_H01         L4
    ## F2S4_160330_016_A01    L2/3 IT
    ## F2S4_160330_016_B01    L2/3 IT
    ## F2S4_160330_016_C01    L2/3 IT
    ## F2S4_160330_016_D01    L2/3 IT
    ## F2S4_160330_016_E01    L2/3 IT
    ## F2S4_160330_016_F01    L2/3 IT
    ## F2S4_160330_016_H01    L2/3 IT
    ## F2S4_160330_017_A01       Sncg
    ## F2S4_160330_017_B01    L2/3 IT
    ## F2S4_160330_017_C01    L2/3 IT
    ## F2S4_160330_017_D01        Vip
    ## F2S4_160330_017_E01    L2/3 IT
    ## F2S4_160330_017_F01    L2/3 IT
    ## F2S4_160330_017_G01    L2/3 IT
    ## F2S4_160330_017_H01    L2/3 IT
    ## F2S4_160330_018_A01         L4
    ## F2S4_160330_018_B01    L2/3 IT
    ## F2S4_160330_018_C01    L2/3 IT
    ## F2S4_160330_018_D01    L2/3 IT
    ## F2S4_160330_018_E01         L4
    ## F2S4_160330_018_F01    L2/3 IT
    ## F2S4_160330_018_G01        Vip
    ## F2S4_160330_018_H01    L2/3 IT
    ## F2S4_160330_019_A01    L2/3 IT
    ## F2S4_160330_019_B01      Lamp5
    ## F2S4_160330_019_C01      Lamp5
    ## F2S4_160330_019_D01        Vip
    ## F2S4_160330_019_E01    L2/3 IT
    ## F2S4_160330_019_F01    L2/3 IT
    ## F2S4_160330_019_G01    L2/3 IT
    ## F2S4_160330_019_H01        Vip
    ## F2S4_160330_020_A01        Vip
    ## F2S4_160330_020_B01      Lamp5
    ## F2S4_160330_020_C01      Lamp5
    ## F2S4_160330_020_D01    L2/3 IT
    ## F2S4_160330_020_E01    L2/3 IT
    ## F2S4_160330_020_F01    L2/3 IT
    ## F2S4_160330_020_G01    L2/3 IT
    ## F2S4_160330_020_H01         L4
    ## F2S4_160330_021_A01      Lamp5
    ## F2S4_160330_021_B01   Serpinf1
    ## F2S4_160330_021_C01         L4
    ## F2S4_160330_021_D01      Lamp5
    ## F2S4_160330_021_E01         L4
    ## F2S4_160330_021_F01         L4
    ## F2S4_160330_021_G01         L4
    ## F2S4_160330_021_H01         L4
    ## F2S4_160330_022_A01         L4
    ## F2S4_160330_022_B01        Sst
    ## F2S4_160330_022_C01        Sst
    ## F2S4_160330_022_D01         L4
    ## F2S4_160330_022_E01        Vip
    ## F2S4_160330_022_F01         L4
    ## F2S4_160330_022_G01      L5 PT
    ## F2S4_160330_022_H01         L4
    ## F2S4_160330_023_A01         L4
    ## F2S4_160330_023_B01         L4
    ## F2S4_160330_023_C01         L4
    ## F2S4_160330_023_D01      L5 IT
    ## F2S4_160330_023_E01         L4
    ## F2S4_160330_023_F01      L5 IT
    ## F2S4_160330_023_G01         L4
    ## F2S4_160330_023_H01         L4
    ## F2S4_160330_024_A01         L4
    ## F2S4_160330_024_B01         L4
    ## F2S4_160330_024_C01         L4
    ## F2S4_160330_024_D01        Sst
    ## F2S4_160330_024_E01      L5 IT
    ## F2S4_160330_024_F01         L4
    ## F2S4_160330_024_G01         L4
    ## F2S4_160330_024_H01         NP
    ## F2S4_160330_025_B01        Sst
    ## F2S4_160330_025_C01         L4
    ## F2S4_160330_025_D01         L4
    ## F2S4_160330_025_E01         L4
    ## F2S4_160330_025_F01         L4
    ## F2S4_160330_025_G01         L4
    ## F2S4_160330_025_H01         NP
    ## F2S4_160330_026_A01      L6 CT
    ## F2S4_160330_026_B01      L6 IT
    ## F2S4_160330_026_C01      L6 IT
    ## F2S4_160330_026_D01      L6 CT
    ## F2S4_160330_026_E01      L6 CT
    ## F2S4_160330_026_F01      L6 IT
    ## F2S4_160330_026_G01      L6 IT
    ## F2S4_160330_026_H01      L6 IT
    ## F2S4_160330_027_A01      L6 CT
    ## F2S4_160330_027_C01      L6 CT
    ## F2S4_160330_027_D01      L6 IT
    ## F2S4_160330_027_E01      L6 IT
    ## F2S4_160330_027_F01      L6 CT
    ## F2S4_160330_027_G01      L6 IT
    ## F2S4_160330_027_H01      L6 CT
    ## F2S4_160330_028_A01      L6 CT
    ## F2S4_160330_028_B01      L6 CT
    ## F2S4_160330_028_C01      L6 CT
    ## F2S4_160330_028_D01      L6 IT
    ## F2S4_160330_028_E01      L6 IT
    ## F2S4_160330_028_F01      L6 CT
    ## F2S4_160330_028_G01      L6 IT
    ## F2S4_160330_028_H01      L6 IT
    ## F2S4_160330_029_A01      L6 IT
    ## F2S4_160330_029_B01      L6 CT
    ## F2S4_160330_029_C01      L6 IT
    ## F2S4_160330_029_D01      L6 IT
    ## F2S4_160330_029_E01      Meis2
    ## F2S4_160330_029_F01      L6 CT
    ## F2S4_160330_029_G01      L6 CT
    ## F2S4_160330_029_H01      L6 IT
    ## F2S4_160330_030_A01      L6 CT
    ## F2S4_160330_030_B01      L6 IT
    ## F2S4_160330_030_C01      L6 IT
    ## F2S4_160330_030_D01      L6 CT
    ## F2S4_160330_030_E01      L6 CT
    ## F2S4_160330_030_F01      L6 CT
    ## F2S4_160330_030_G01      L6 CT
    ## F2S4_160330_030_H01      L6 CT
    ## F2S4_160331_001_A01      Pvalb
    ## F2S4_160331_001_B01      Pvalb
    ## F2S4_160331_001_C01        Sst
    ## F2S4_160331_001_D01      Pvalb
    ## F2S4_160331_001_E01      Pvalb
    ## F2S4_160331_001_F01      Pvalb
    ## F2S4_160331_001_G01      Pvalb
    ## F2S4_160331_001_H01      Pvalb
    ## F2S4_160331_002_A01      Pvalb
    ## F2S4_160331_002_B01      Pvalb
    ## F2S4_160331_002_C01      Pvalb
    ## F2S4_160331_002_D01      Pvalb
    ## F2S4_160331_002_E01      Pvalb
    ## F2S4_160331_002_F01      Pvalb
    ## F2S4_160401_001_A01      Pvalb
    ## F2S4_160401_001_B01      Pvalb
    ## F2S4_160401_001_C01      Pvalb
    ## F2S4_160401_001_D01      Pvalb
    ## F2S4_160401_001_E01      Pvalb
    ## F2S4_160401_001_F01      Pvalb
    ## F2S4_160401_001_G01      Pvalb
    ## F2S4_160401_001_H01      Pvalb
    ## F2S4_160401_002_A01      Pvalb
    ## F2S4_160401_002_B01      Pvalb
    ## F2S4_160401_002_C01      Pvalb
    ## F2S4_160401_002_D01      Pvalb
    ## F2S4_160401_002_E01      Pvalb
    ## F2S4_160401_002_F01      Pvalb
    ## F2S4_160401_002_G01      Pvalb
    ## F2S4_160401_002_H01      Pvalb
    ## F2S4_160401_003_A01      Pvalb
    ## F2S4_160401_003_B01      Pvalb
    ## F2S4_160401_003_C01      Pvalb
    ## F2S4_160401_003_D01      Pvalb
    ## F2S4_160401_003_E01      Pvalb
    ## F2S4_160401_003_F01      Pvalb
    ## F2S4_160401_003_G01      Pvalb
    ## F2S4_160401_003_H01      Pvalb
    ## F2S4_160401_004_A01      Pvalb
    ## F2S4_160401_004_B01      Pvalb
    ## F2S4_160401_004_C01      Pvalb
    ## F2S4_160401_004_D01      Pvalb
    ## F2S4_160401_004_E01      Pvalb
    ## F2S4_160404_001_A01      Lamp5
    ## F2S4_160404_001_B01      Lamp5
    ## F2S4_160404_001_C01        Vip
    ## F2S4_160404_001_D01        Vip
    ## F2S4_160404_001_E01      Lamp5
    ## F2S4_160404_001_F01      Lamp5
    ## F2S4_160404_001_G01      Lamp5
    ## F2S4_160404_001_H01        Vip
    ## F2S4_160404_002_A01       Sncg
    ## F2S4_160404_002_B01        Vip
    ## F2S4_160404_002_D01      Lamp5
    ## F2S4_160404_002_E01      Lamp5
    ## F2S4_160404_002_F01      Lamp5
    ## F2S4_160404_002_G01        Vip
    ## F2S4_160404_002_H01      Lamp5
    ## F2S4_160404_003_A01        Vip
    ## F2S4_160404_003_B01      Lamp5
    ## F2S4_160404_003_C01      Lamp5
    ## F2S4_160404_003_D01      Lamp5
    ## F2S4_160404_003_E01      Lamp5
    ## F2S4_160404_003_F01      Lamp5
    ## F2S4_160404_003_H01      Lamp5
    ## F2S4_160404_004_A01      Lamp5
    ## F2S4_160404_004_B01      Lamp5
    ## F2S4_160404_004_C01      Lamp5
    ## F2S4_160404_004_D01      Lamp5
    ## F2S4_160404_004_E01      Lamp5
    ## F2S4_160404_004_F01      Lamp5
    ## F2S4_160404_004_G01      Lamp5
    ## F2S4_160404_004_H01       Sncg
    ## F2S4_160404_005_A01      Lamp5
    ## F2S4_160404_005_B01        Vip
    ## F2S4_160404_005_C01      Lamp5
    ## F2S4_160404_005_D01      Lamp5
    ## F2S4_160404_005_E01      Lamp5
    ## F2S4_160404_005_F01      Lamp5
    ## F2S4_160404_005_G01      Lamp5
    ## F2S4_160404_005_H01        Vip
    ## F2S4_160404_006_A01        Sst
    ## F2S4_160404_006_B01      Pvalb
    ## F2S4_160404_006_C01        Vip
    ## F2S4_160404_006_D01      Lamp5
    ## F2S4_160404_006_E01        Vip
    ## F2S4_160404_006_F01    L2/3 IT
    ## F2S4_160404_006_G01      Lamp5
    ## F2S4_160404_006_H01      Lamp5
    ## F2S4_160404_007_A01      Lamp5
    ## F2S4_160404_007_B01      Lamp5
    ## F2S4_160404_007_C01      Astro
    ## F2S4_160404_007_D01        Sst
    ## F2S4_160404_007_E01        Sst
    ## F2S4_160404_007_F01        Vip
    ## F2S4_160404_007_G01        Vip
    ## F2S4_160404_007_H01        Vip
    ## F2S4_160404_008_A01      Pvalb
    ## F2S4_160404_008_B01      Pvalb
    ## F2S4_160404_008_C01        Sst
    ## F2S4_160404_008_D01        Sst
    ## F2S4_160404_008_E01        Sst
    ## F2S4_160404_008_F01      Lamp5
    ## F2S4_160404_008_G01        Sst
    ## F2S4_160404_008_H01        Sst
    ## F2S4_160404_009_A01        Sst
    ## F2S4_160404_009_B01        Vip
    ## F2S4_160404_009_C01        Sst
    ## F2S4_160404_009_D01        Sst
    ## F2S4_160404_009_E01        Sst
    ## F2S4_160404_009_F01      Lamp5
    ## F2S4_160404_009_G01        Vip
    ## F2S4_160404_009_H01        Sst
    ## F2S4_160404_010_A01      Astro
    ## F2S4_160404_010_B01        Sst
    ## F2S4_160404_010_C01      Pvalb
    ## F2S4_160404_010_D01      Lamp5
    ## F2S4_160404_010_E01        Sst
    ## F2S4_160404_010_F01      Pvalb
    ## F2S4_160404_010_G01        Sst
    ## F2S4_160404_010_H01      Pvalb
    ## F2S4_160404_011_A01      Pvalb
    ## F2S4_160404_011_B01        Sst
    ## F2S4_160404_011_C01        Sst
    ## F2S4_160404_011_D01        Sst
    ## F2S4_160404_011_E01      Lamp5
    ## F2S4_160404_011_F01        Sst
    ## F2S4_160404_011_G01      Astro
    ## F2S4_160404_011_H01      Pvalb
    ## F2S4_160404_012_A01      Pvalb
    ## F2S4_160404_012_B01        Vip
    ## F2S4_160404_012_C01      Pvalb
    ## F2S4_160404_012_D01        Vip
    ## F2S4_160404_012_E01      Astro
    ## F2S4_160404_012_F01      Pvalb
    ## F2S4_160404_012_G01      Pvalb
    ## F2S4_160404_012_H01        Sst
    ## F2S4_160404_013_A01        Vip
    ## F2S4_160404_013_B01        Vip
    ## F2S4_160404_013_C01      Pvalb
    ## F2S4_160404_013_D01        Vip
    ## F2S4_160404_013_E01      Lamp5
    ## F2S4_160404_013_F01        Vip
    ## F2S4_160404_013_G01      Lamp5
    ## F2S4_160404_013_H01      Lamp5
    ## F2S4_160404_014_A01      Pvalb
    ## F2S4_160404_014_B01        Vip
    ## F2S4_160404_014_C01       Sncg
    ## F2S4_160404_014_D01        Vip
    ## F2S4_160404_014_E01    L2/3 IT
    ## F2S4_160404_014_F01        Vip
    ## F2S4_160404_014_G01        Vip
    ## F2S4_160404_014_H01        Sst
    ## F2S4_160404_015_A01    L2/3 IT
    ## F2S4_160404_015_B01      Pvalb
    ## F2S4_160404_015_C01        Vip
    ## F2S4_160404_015_D01        Vip
    ## F2S4_160404_015_E01        Vip
    ## F2S4_160404_015_F01        Vip
    ## F2S4_160404_015_G01      Pvalb
    ## F2S4_160404_015_H01        Vip
    ## F2S4_160404_016_A01        Vip
    ## F2S4_160404_016_B01        Vip
    ## F2S4_160404_016_C01        Vip
    ## F2S4_160404_016_D01        Vip
    ## F2S4_160404_016_E01        Vip
    ## F2S4_160404_016_F01        Vip
    ## F2S4_160404_016_G01      Lamp5
    ## F2S4_160404_016_H01        Vip
    ## F2S4_160404_017_A01      Lamp5
    ## F2S4_160404_017_B01    L2/3 IT
    ## F2S4_160404_017_C01    L2/3 IT
    ## F2S4_160404_017_D01      Lamp5
    ## F2S4_160404_017_E01      Lamp5
    ## F2S4_160404_017_F01        Vip
    ## F2S4_160404_017_G01        Sst
    ## F2S4_160404_017_H01      Pvalb
    ## F2S4_160404_018_A01      Pvalb
    ## F2S4_160404_018_B01        Vip
    ## F2S4_160404_018_C01        Vip
    ## F2S4_160404_018_D01        Vip
    ## F2S4_160404_018_E01        Vip
    ## F2S4_160404_018_F01      Lamp5
    ## F2S4_160404_018_G01        Vip
    ## F2S4_160404_018_H01        Vip
    ## F2S4_160404_019_A01      Lamp5
    ## F2S4_160404_019_B01      Lamp5
    ## F2S4_160404_019_C01      Pvalb
    ## F2S4_160404_019_D01        Vip
    ## F2S4_160404_019_E01       Sncg
    ## F2S4_160404_019_F01        Vip
    ## F2S4_160404_019_G01        Vip
    ## F2S4_160404_019_H01        Vip
    ## F2S4_160404_020_A01        Vip
    ## F2S4_160404_020_B01      Pvalb
    ## F2S4_160404_020_C01        Vip
    ## F2S4_160404_020_D01      Lamp5
    ## F2S4_160404_020_E01        Vip
    ## F2S4_160404_020_F01        Vip
    ## F2S4_160404_020_G01        Sst
    ## F2S4_160404_020_H01      Lamp5
    ## F2S4_160404_021_A01      Pvalb
    ## F2S4_160404_021_B01        Vip
    ## F2S4_160404_021_C01        Vip
    ## F2S4_160404_021_D01        Vip
    ## F2S4_160404_021_E01        Vip
    ## F2S4_160404_021_F01        Vip
    ## F2S4_160404_021_G01        Vip
    ## F2S4_160404_021_H01        Vip
    ## F2S4_160404_022_A01      Pvalb
    ## F2S4_160404_022_B01        Sst
    ## F2S4_160404_022_D01      Lamp5
    ## F2S4_160404_022_E01        Vip
    ## F2S4_160404_022_F01        Sst
    ## F2S4_160404_022_G01      Pvalb
    ## F2S4_160404_022_H01      Pvalb
    ## F2S4_160404_023_A01      Astro
    ## F2S4_160404_023_B01      Pvalb
    ## F2S4_160404_023_C01        Sst
    ## F2S4_160404_023_D01        Vip
    ## F2S4_160404_023_E01      Pvalb
    ## F2S4_160404_023_G01      Pvalb
    ## F2S4_160404_023_H01      Pvalb
    ## F2S4_160404_024_A01      Pvalb
    ## F2S4_160404_024_B01      Pvalb
    ## F2S4_160404_024_C01        Vip
    ## F2S4_160404_024_D01        Sst
    ## F2S4_160404_024_E01        Vip
    ## F2S4_160404_024_F01        Vip
    ## F2S4_160404_024_G01      Pvalb
    ## F2S4_160404_024_H01      Pvalb
    ## F2S4_160404_025_A01      Lamp5
    ## F2S4_160404_025_B01        Vip
    ## F2S4_160404_025_C01      Pvalb
    ## F2S4_160404_025_D01        Vip
    ## F2S4_160404_025_E01        Sst
    ## F2S4_160404_025_G01        Vip
    ## F2S4_160404_025_H01        Sst
    ## F2S4_160404_026_A01        Vip
    ## F2S4_160404_026_B01        Sst
    ## F2S4_160404_026_C01        Vip
    ## F2S4_160404_026_D01        Vip
    ## F2S4_160404_026_E01      Meis2
    ## F2S4_160404_026_F01      Lamp5
    ## F2S4_160404_026_G01        Sst
    ## F2S4_160404_026_H01        Vip
    ## F2S4_160404_027_A01        Sst
    ## F2S4_160404_027_B01      Astro
    ## F2S4_160404_027_C01        Sst
    ## F2S4_160404_027_D01      Pvalb
    ## F2S4_160404_027_E01      Pvalb
    ## F2S4_160404_027_F01        Vip
    ## F2S4_160404_027_G01      Lamp5
    ## F2S4_160404_028_C01      Lamp5
    ## F2S4_160404_028_D01        Vip
    ## F2S4_160404_028_F01      Lamp5
    ## F2S4_160404_028_G01        Sst
    ## F2S4_160404_028_H01        Sst
    ## F2S4_160404_029_C01      Lamp5
    ## F2S4_160404_029_D01       Sncg
    ## F2S4_160404_029_E01      Lamp5
    ## F2S4_160404_029_F01        Vip
    ## F2S4_160404_029_G01      Lamp5
    ## F2S4_160404_029_H01      Lamp5
    ## F2S4_160404_030_A01        Vip
    ## F2S4_160404_030_B01      Oligo
    ## F2S4_160404_030_C01        Vip
    ## F2S4_160404_030_D01      Pvalb
    ## F2S4_160404_030_E01      Lamp5
    ## F2S4_160404_030_F01      Lamp5
    ## F2S4_160404_030_G01        Sst
    ## F2S4_160404_030_H01        Vip
    ## F2S4_160404_031_A01      Lamp5
    ## F2S4_160404_031_B01      Lamp5
    ## F2S4_160404_031_C01        Vip
    ## F2S4_160404_031_D01      Lamp5
    ## F2S4_160404_031_E01      Lamp5
    ## F2S4_160404_031_F01        Vip
    ## F2S4_160404_031_G01       Sncg
    ## F2S4_160404_031_H01       Sncg
    ## F2S4_160404_032_A01      Lamp5
    ## F2S4_160404_032_B01    L2/3 IT
    ## F2S4_160404_032_C01       Sncg
    ## F2S4_160404_032_D01      Lamp5
    ## F2S4_160404_032_E01      Lamp5
    ## F2S4_160404_032_F01      Lamp5
    ## F2S4_160404_032_G01    L2/3 IT
    ## F2S4_160404_032_H01      Pvalb
    ## F2S4_160404_033_A01        Sst
    ## F2S4_160404_033_B01        Vip
    ## F2S4_160404_033_C01    L2/3 IT
    ## F2S4_160404_033_D01      Pvalb
    ## F2S4_160404_033_E01        Vip
    ## F2S4_160404_033_F01        Vip
    ## F2S4_160404_033_G01      Lamp5
    ## F2S4_160404_033_H01      Lamp5
    ## F2S4_160404_034_A01        Vip
    ## F2S4_160404_034_B01        Vip
    ## F2S4_160404_034_C01      Lamp5
    ## F2S4_160404_034_D01        Vip
    ## F2S4_160404_034_E01      Lamp5
    ## F2S4_160404_034_F01        Vip
    ## F2S4_160404_034_G01      Lamp5
    ## F2S4_160404_034_H01      Lamp5
    ## F2S4_160404_035_A01      Lamp5
    ## F2S4_160404_035_B01        Vip
    ## F2S4_160404_035_C01      Lamp5
    ## F2S4_160404_035_D01      Lamp5
    ## F2S4_160404_035_E01      Lamp5
    ## F2S4_160404_035_F01      Lamp5
    ## F2S4_160404_035_G01      Lamp5
    ## F2S4_160404_035_H01      Lamp5
    ## F2S4_160404_036_A01        Vip
    ## F2S4_160404_036_B01      Lamp5
    ## F2S4_160404_036_C01        Vip
    ## F2S4_160404_036_D01        Vip
    ## F2S4_160404_036_E01      Lamp5
    ## F2S4_160404_037_A01      Lamp5
    ## F2S4_160404_037_B01 Macrophage
    ## F2S4_160404_037_C01      Lamp5
    ## F2S4_160404_037_D01      Lamp5
    ## F2S4_160404_037_E01      Lamp5
    ## F2S4_160404_037_F01 Macrophage
    ## F2S4_160404_037_G01      Lamp5
    ## F2S4_160404_037_H01      Astro
    ## F2S4_160404_038_A01        Vip
    ## F2S4_160404_038_B01        Sst
    ## F2S4_160404_038_C01      Astro
    ## F2S4_160404_038_D01      Astro
    ## F2S4_160404_038_E01        Sst
    ## F2S4_160404_038_F01        Sst
    ## F2S4_160404_038_G01        Sst
    ## F2S4_160404_038_H01        Sst
    ## F2S4_160404_039_A01        Sst
    ## F2S4_160404_039_B01        Sst
    ## F2S4_160404_039_C01        Sst
    ## F2S4_160404_039_D01        Sst
    ## F2S4_160404_039_E01        Vip
    ## F2S4_160404_039_F01        Vip
    ## F2S4_160404_039_G01        Vip
    ## F2S4_160404_039_H01        Vip
    ## F2S4_160404_040_A01        Sst
    ## F2S4_160404_040_C01      Astro
    ## F2S4_160404_040_D01        Vip
    ## F2S4_160404_040_E01      Lamp5
    ## F2S4_160404_040_F01        Vip
    ## F2S4_160404_040_H01        Sst
    ## F2S4_160404_041_A01      Pvalb
    ## F2S4_160404_041_B01        Sst
    ## F2S4_160404_041_C01        Vip
    ## F2S4_160404_041_D01        Sst
    ## F2S4_160404_041_E01      Lamp5
    ## F2S4_160404_041_F01        Sst
    ## F2S4_160404_041_G01        Sst
    ## F2S4_160404_041_H01        Sst
    ## F2S4_160404_042_A01        Vip
    ## F2S4_160404_042_B01      Astro
    ## F2S4_160404_042_C01        Sst
    ## F2S4_160404_042_D01        Sst
    ## F2S4_160404_042_E01      Lamp5
    ## F2S4_160404_042_F01        Sst
    ## F2S4_160404_042_G01      Lamp5
    ## F2S4_160404_042_H01        Sst
    ## F2S4_160404_043_A01        Vip
    ## F2S4_160404_043_B01      Pvalb
    ## F2S4_160404_043_C01        Vip
    ## F2S4_160404_043_D01        Vip
    ## F2S4_160404_043_E01    L2/3 IT
    ## F2S4_160404_043_F01      Pvalb
    ## F2S4_160404_043_G01        Sst
    ## F2S4_160404_043_H01      Pvalb
    ## F2S4_160404_044_A01        Vip
    ## F2S4_160404_044_B01      Lamp5
    ## F2S4_160404_044_C01        Vip
    ## F2S4_160404_044_D01      Astro
    ## F2S4_160404_044_E01        Vip
    ## F2S4_160404_044_F01        Vip
    ## F2S4_160404_044_G01       Sncg
    ## F2S4_160404_044_H01        Vip
    ## F2S4_160404_045_A01      Lamp5
    ## F2S4_160404_045_B01        Vip
    ## F2S4_160404_045_C01        Vip
    ## F2S4_160404_045_D01        Vip
    ## F2S4_160404_045_E01        Vip
    ## F2S4_160404_045_F01        Vip
    ## F2S4_160404_045_G01      Lamp5
    ## F2S4_160404_045_H01        Vip
    ## F2S4_160404_046_A01        Vip
    ## F2S4_160404_046_B01        Vip
    ## F2S4_160404_046_C01        Vip
    ## F2S4_160404_046_D01        Vip
    ## F2S4_160404_046_E01      Lamp5
    ## F2S4_160404_046_F01      Lamp5
    ## F2S4_160404_046_G01        Vip
    ## F2S4_160404_046_H01    L2/3 IT
    ## F2S4_160404_047_A01      Lamp5
    ## F2S4_160404_047_B01        Sst
    ## F2S4_160404_047_C01        Vip
    ## F2S4_160404_047_D01        Vip
    ## F2S4_160404_047_E01      Lamp5
    ## F2S4_160404_047_F01        Vip
    ## F2S4_160404_047_G01        Vip
    ## F2S4_160404_047_H01        Vip
    ## F2S4_160404_048_A01      Lamp5
    ## F2S4_160404_048_B01        Vip
    ## F2S4_160404_048_C01      Lamp5
    ## F2S4_160404_048_D01      Pvalb
    ## F2S4_160404_048_E01        Vip
    ## F2S4_160404_048_F01        Sst
    ## F2S4_160404_048_G01        Vip
    ## F2S4_160404_048_H01        Vip
    ## F2S4_160404_049_A01        Sst
    ## F2S4_160404_049_B01        Vip
    ## F2S4_160404_049_C01      Lamp5
    ## F2S4_160404_049_D01       Sncg
    ## F2S4_160404_049_F01        Vip
    ## F2S4_160404_049_H01      Lamp5
    ## F2S4_160404_050_A01      Lamp5
    ## F2S4_160404_050_B01      Lamp5
    ## F2S4_160404_050_C01      Lamp5
    ## F2S4_160404_050_D01      Lamp5
    ## F2S4_160404_050_E01      Lamp5
    ## F2S4_160404_050_F01        Sst
    ## F2S4_160404_050_G01        Vip
    ## F2S4_160404_050_H01      Meis2
    ## F2S4_160404_051_A01        Sst
    ## F2S4_160404_051_B01      Lamp5
    ## F2S4_160404_051_C01      Astro
    ## F2S4_160404_051_D01      Lamp5
    ## F2S4_160404_051_E01      Astro
    ## F2S4_160404_051_F01        Sst
    ## F2S4_160404_051_G01      Astro
    ## F2S4_160404_051_H01      Lamp5
    ## F2S4_160405_001_B01       Endo
    ## F2S4_160405_001_C01      Lamp5
    ## F2S4_160405_001_D01       Endo
    ## F2S4_160405_001_E01       VLMC
    ## F2S4_160405_001_F01      Lamp5
    ## F2S4_160405_001_H01       Endo
    ## F2S4_160405_002_A01      Lamp5
    ## F2S4_160405_002_B01      Lamp5
    ## F2S4_160405_002_C01       VLMC
    ## F2S4_160405_002_D01      Lamp5
    ## F2S4_160405_002_E01      Lamp5
    ## F2S4_160405_002_F01      Lamp5
    ## F2S4_160405_002_G01       VLMC
    ## F2S4_160405_002_H01       Endo
    ## F2S4_160405_003_A01      Lamp5
    ## F2S4_160405_003_C01      Lamp5
    ## F2S4_160405_003_D01       VLMC
    ## F2S4_160405_003_F01      Lamp5
    ## F2S4_160405_003_G01      Lamp5
    ## F2S4_160405_003_H01       Endo
    ## F2S4_160405_009_A01       VLMC
    ## F2S4_160405_009_B01      Lamp5
    ## F2S4_160405_009_C01      Lamp5
    ## F2S4_160405_009_D01      Lamp5
    ## F2S4_160405_009_E01      Lamp5
    ## F2S4_160405_009_F01      Lamp5
    ## F2S4_160405_009_G01      Lamp5
    ## F2S4_160405_009_H01       Endo
    ## F2S4_160405_010_B01      Lamp5
    ## F2S4_160405_010_C01       VLMC
    ## F2S4_160405_010_D01      Lamp5
    ## F2S4_160405_010_E01       Endo
    ## F2S4_160405_010_F01       VLMC
    ## F2S4_160405_010_G01      Lamp5
    ## F2S4_160405_010_H01       Endo
    ## F2S4_160405_011_A01      Lamp5
    ## F2S4_160405_011_C01      Lamp5
    ## F2S4_160405_011_D01      Lamp5
    ## F2S4_160405_011_E01      Lamp5
    ## F2S4_160405_011_F01      Lamp5
    ## F2S4_160405_011_G01      Lamp5
    ## F2S4_160405_011_H01       Endo
    ## F2S4_160405_012_B01      Lamp5
    ## F2S4_160405_012_C01       Endo
    ## F2S4_160405_012_D01       VLMC
    ## F2S4_160405_012_E01      Lamp5
    ## F2S4_160405_012_F01       VLMC
    ## F2S4_160405_012_G01       VLMC
    ## F2S4_160405_012_H01       Endo
    ## F2S4_160405_014_A01       VLMC
    ## F2S4_160405_014_B01       VLMC
    ## F2S4_160405_014_C01      Lamp5
    ## F2S4_160405_014_D01       VLMC
    ## F2S4_160405_014_E01      Lamp5
    ## F2S4_160405_014_F01      Lamp5
    ## F2S4_160405_014_G01      Lamp5
    ## F2S4_160405_014_H01       VLMC
    ## F2S4_160405_015_A01      Lamp5
    ## F2S4_160405_015_B01       Endo
    ## F2S4_160405_015_C01       VLMC
    ## F2S4_160405_015_E01       VLMC
    ## F2S4_160405_015_F01      Lamp5
    ## F2S4_160405_015_H01      Lamp5
    ## F2S4_160405_016_A01      Lamp5
    ## F2S4_160405_016_B01       VLMC
    ## F2S4_160405_016_C01       VLMC
    ## F2S4_160405_016_D01       VLMC
    ## F2S4_160405_016_F01       VLMC
    ## F2S4_160405_016_G01      Lamp5
    ## F2S4_160405_016_H01       Endo
    ## F2S4_160405_017_A01       VLMC
    ## F2S4_160405_017_B01       VLMC
    ## F2S4_160405_017_C01      Lamp5
    ## F2S4_160405_017_D01      Lamp5
    ## F2S4_160405_017_E01      Lamp5
    ## F2S4_160405_017_F01       VLMC
    ## F2S4_160405_017_G01       VLMC
    ## F2S4_160405_017_H01       VLMC
    ## F2S4_160405_020_A01      Lamp5
    ## F2S4_160405_020_B01      Lamp5
    ## F2S4_160405_020_C01       VLMC
    ## F2S4_160405_020_D01      Lamp5
    ## F2S4_160405_020_F01      Lamp5
    ## F2S4_160405_020_G01      Lamp5
    ## F2S4_160405_020_H01       VLMC
    ## F2S4_160405_021_A01       Endo
    ## F2S4_160405_021_B01       VLMC
    ## F2S4_160405_021_C01      Lamp5
    ## F2S4_160405_021_D01      Lamp5
    ## F2S4_160405_021_E01       VLMC
    ## F2S4_160405_021_F01       Endo
    ## F2S4_160405_021_G01      Lamp5
    ## F2S4_160405_021_H01       Endo
    ## F2S4_160405_022_A01       VLMC
    ## F2S4_160405_022_B01       Endo
    ## F2S4_160405_022_C01       VLMC
    ## F2S4_160405_022_D01       Endo
    ## F2S4_160405_022_E01       Endo
    ## F2S4_160405_022_F01      Lamp5
    ## F2S4_160405_022_G01       VLMC
    ## F2S4_160405_022_H01       VLMC
    ## F2S4_160405_023_B01       VLMC
    ## F2S4_160405_023_C01       VLMC
    ## F2S4_160405_023_D01         CR
    ## F2S4_160405_023_E01      Lamp5
    ## F2S4_160405_023_F01      Lamp5
    ## F2S4_160405_023_G01      Lamp5
    ## F2S4_160405_023_H01       Endo
    ## F2S4_160405_024_A01      Lamp5
    ## F2S4_160405_024_B01        Vip
    ## F2S4_160405_024_C01        Sst
    ## F2S4_160405_024_D01        Sst
    ## F2S4_160405_024_E01        Vip
    ## F2S4_160405_024_F01        Vip
    ## F2S4_160405_024_G01        Sst
    ## F2S4_160405_024_H01      Astro
    ## F2S4_160405_025_A01        Sst
    ## F2S4_160405_025_B01        Sst
    ## F2S4_160405_025_C01        Sst
    ## F2S4_160405_025_E01        Sst
    ## F2S4_160405_025_F01        Sst
    ## F2S4_160405_025_G01        Sst
    ## F2S4_160405_025_H01        Sst
    ## F2S4_160405_026_B01      Pvalb
    ## F2S4_160405_026_C01      Pvalb
    ## F2S4_160405_026_D01        Vip
    ## F2S4_160405_026_E01      Astro
    ## F2S4_160405_026_F01      Pvalb
    ## F2S4_160405_026_G01        Sst
    ## F2S4_160405_026_H01      Pvalb
    ## F2S4_160405_027_A01        Sst
    ## F2S4_160405_027_B01      Pvalb
    ## F2S4_160405_027_C01       Sncg
    ## F2S4_160405_027_D01        Vip
    ## F2S4_160405_027_E01      Astro
    ## F2S4_160405_027_F01        Sst
    ## F2S4_160405_027_G01        Sst
    ## F2S4_160405_027_H01        Sst
    ## F2S4_160405_028_A01      Pvalb
    ## F2S4_160405_028_B01        Sst
    ## F2S4_160405_028_C01      Pvalb
    ## F2S4_160405_028_D01      Pvalb
    ## F2S4_160405_028_E01        Sst
    ## F2S4_160405_028_F01        Vip
    ## F2S4_160405_028_G01        Sst
    ## F2S4_160405_028_H01        Sst
    ## F2S4_160405_029_A01        Sst
    ## F2S4_160405_029_B01      Pvalb
    ## F2S4_160405_029_D01      Pvalb
    ## F2S4_160405_029_E01      Pvalb
    ## F2S4_160405_029_G01        Sst
    ## F2S4_160405_029_H01        Sst
    ## F2S4_160405_030_A01        Sst
    ## F2S4_160405_030_B01        Sst
    ## F2S4_160405_030_C01        Sst
    ## F2S4_160405_030_D01        Vip
    ## F2S4_160405_030_E01        Sst
    ## F2S4_160405_030_F01        Vip
    ## F2S4_160405_030_G01        Sst
    ## F2S4_160405_030_H01      Meis2
    ## F2S4_160405_031_A01      Lamp5
    ## F2S4_160405_031_B01      Lamp5
    ## F2S4_160405_031_C01        Vip
    ## F2S4_160405_031_D01        Vip
    ## F2S4_160405_031_E01        Vip
    ## F2S4_160405_031_F01      Lamp5
    ## F2S4_160405_031_G01      Lamp5
    ## F2S4_160405_031_H01      Lamp5
    ## F2S4_160405_032_A01      Lamp5
    ## F2S4_160405_032_B01       Sncg
    ## F2S4_160405_032_C01      Lamp5
    ## F2S4_160405_032_D01      Lamp5
    ## F2S4_160405_032_E01        Vip
    ## F2S4_160405_032_F01      Lamp5
    ## F2S4_160405_032_G01      Lamp5
    ## F2S4_160405_032_H01      Lamp5
    ## F2S4_160405_033_A01      Lamp5
    ## F2S4_160405_033_B01      Lamp5
    ## F2S4_160405_033_C01      Lamp5
    ## F2S4_160405_033_D01      Lamp5
    ## F2S4_160405_033_F01        Vip
    ## F2S4_160405_033_G01        Vip
    ## F2S4_160405_033_H01        Vip
    ## F2S4_160405_034_A01    L2/3 IT
    ## F2S4_160405_034_B01        Vip
    ## F2S4_160405_034_C01        Vip
    ## F2S4_160405_034_D01      Lamp5
    ## F2S4_160405_034_E01      Lamp5
    ## F2S4_160405_034_F01      Lamp5
    ## F2S4_160405_034_G01        Vip
    ## F2S4_160405_034_H01      Astro
    ## F2S4_160405_035_A01        Vip
    ## F2S4_160405_035_B01        Vip
    ## F2S4_160405_035_C01        Vip
    ## F2S4_160405_035_D01        Vip
    ## F2S4_160405_035_E01        Vip
    ## F2S4_160405_035_F01        Sst
    ## F2S4_160405_035_G01        Vip
    ## F2S4_160405_035_H01      Lamp5
    ## F2S4_160405_036_A01        Vip
    ## F2S4_160405_036_B01        Sst
    ## F2S4_160405_036_C01        Vip
    ## F2S4_160405_036_D01        Vip
    ## F2S4_160405_036_E01        Sst
    ## F2S4_160405_036_F01        Vip
    ## F2S4_160405_036_G01        Vip
    ## F2S4_160405_036_H01        Vip
    ## F2S4_160405_037_A01        Sst
    ## F2S4_160405_037_B01        Vip
    ## F2S4_160405_037_C01        Vip
    ## F2S4_160405_037_D01        Vip
    ## F2S4_160405_037_E01        Vip
    ## F2S4_160405_037_F01      Lamp5
    ## F2S4_160405_037_G01      Lamp5
    ## F2S4_160405_037_H01      Lamp5
    ## F2S4_160405_038_A01      Lamp5
    ## F2S4_160405_038_B01        Vip
    ## F2S4_160405_038_C01        Vip
    ## F2S4_160405_038_D01        Vip
    ## F2S4_160405_038_E01        Vip
    ## F2S4_160405_038_F01        Vip
    ## F2S4_160405_038_G01        Vip
    ## F2S4_160405_038_H01      Lamp5
    ## F2S4_160405_039_A01      Lamp5
    ## F2S4_160405_039_B01      Lamp5
    ## F2S4_160405_039_C01        Vip
    ## F2S4_160405_039_D01        Vip
    ## F2S4_160405_039_F01        Vip
    ## F2S4_160405_039_G01        Vip
    ## F2S4_160405_039_H01        Vip
    ## F2S4_160405_040_A01        Sst
    ## F2S4_160405_040_B01        Vip
    ## F2S4_160405_040_C01        Vip
    ## F2S4_160405_040_D01      Lamp5
    ## F2S4_160405_040_E01      Astro
    ## F2S4_160405_040_F01        Vip
    ## F2S4_160405_040_G01        Sst
    ## F2S4_160405_040_H01      Astro
    ## F2S4_160405_041_A01        Vip
    ## F2S4_160405_041_B01      Pvalb
    ## F2S4_160405_041_C01        Vip
    ## F2S4_160405_041_D01        Vip
    ## F2S4_160405_041_E01      Astro
    ## F2S4_160405_041_F01        Vip
    ## F2S4_160405_041_G01        Vip
    ## F2S4_160405_041_H01        Vip
    ## F2S4_160405_042_B01      Lamp5
    ## F2S4_160405_042_C01      Lamp5
    ## F2S4_160405_042_D01      Lamp5
    ## F2S4_160405_042_E01      Lamp5
    ## F2S4_160405_042_F01      Lamp5
    ## F2S4_160405_042_G01      Lamp5
    ## F2S4_160405_042_H01      Astro
    ## F2S4_160405_043_A01      Lamp5
    ## F2S4_160405_043_B01      Lamp5
    ## F2S4_160405_043_C01      Lamp5
    ## F2S4_160405_043_D01      Lamp5
    ## F2S4_160405_043_E01      Lamp5
    ## F2S4_160405_043_F01      Lamp5
    ## F2S4_160405_043_G01      Lamp5
    ## F2S4_160405_043_H01      Lamp5
    ## F2S4_160405_044_A01      Lamp5
    ## F2S4_160405_044_B01        Vip
    ## F2S4_160405_044_C01      Lamp5
    ## F2S4_160405_044_D01      Lamp5
    ## F2S4_160405_044_E01      Lamp5
    ## F2S4_160405_044_F01        Vip
    ## F2S4_160405_044_G01      Lamp5
    ## F2S4_160405_044_H01      Lamp5
    ## F2S4_160406_001_A01        Vip
    ## F2S4_160406_001_B01        Vip
    ## F2S4_160406_001_C01        Vip
    ## F2S4_160406_001_D01        Vip
    ## F2S4_160406_001_E01        Vip
    ## F2S4_160406_001_F01        Vip
    ## F2S4_160406_001_G01        Vip
    ## F2S4_160406_001_H01      Lamp5
    ## F2S4_160406_002_A01        Vip
    ## F2S4_160406_002_B01        Vip
    ## F2S4_160406_002_C01       Sncg
    ## F2S4_160406_002_D01        Vip
    ## F2S4_160406_002_E01        Vip
    ## F2S4_160406_002_F01        Vip
    ## F2S4_160406_002_G01        Vip
    ## F2S4_160406_002_H01        Vip
    ## F2S4_160406_003_A01        Vip
    ## F2S4_160406_003_B01        Vip
    ## F2S4_160406_003_C01        Vip
    ## F2S4_160406_003_D01        Vip
    ## F2S4_160406_003_E01        Vip
    ## F2S4_160408_001_A01      Pvalb
    ## F2S4_160408_001_B01      Pvalb
    ## F2S4_160408_001_C01      Pvalb
    ## F2S4_160408_001_D01      Pvalb
    ## F2S4_160408_001_E01      Pvalb
    ## F2S4_160408_001_F01      Pvalb
    ## F2S4_160408_001_G01      Pvalb
    ## F2S4_160408_001_H01      Pvalb
    ## F2S4_160408_002_A01      Pvalb
    ## F2S4_160408_002_B01      Pvalb
    ## F2S4_160408_002_C01      Pvalb
    ## F2S4_160408_002_D01      Pvalb
    ## F2S4_160408_002_E01      Pvalb
    ## F2S4_160408_002_F01      Pvalb
    ## F2S4_160408_002_G01      Pvalb
    ## F2S4_160408_002_H01      Pvalb
    ## F2S4_160408_003_A01      Pvalb
    ## F2S4_160408_003_B01      Pvalb
    ## F2S4_160408_003_C01      Pvalb
    ## F2S4_160408_003_D01      Pvalb
    ## F2S4_160408_003_E01      Pvalb
    ## F2S4_160408_003_F01      Pvalb
    ## F2S4_160408_003_G01      Pvalb
    ## F2S4_160408_003_H01      Pvalb
    ## F2S4_160408_004_A01      Pvalb
    ## F2S4_160411_001_B01      Pvalb
    ## F2S4_160411_001_C01      Pvalb
    ## F2S4_160411_001_D01      Pvalb
    ## F2S4_160412_007_A01        Sst
    ## F2S4_160412_007_B01        Sst
    ## F2S4_160412_007_C01      L5 IT
    ## F2S4_160412_007_D01        Sst
    ## F2S4_160412_007_E01        Sst
    ## F2S4_160412_007_F01        Sst
    ## F2S4_160412_007_G01        Sst
    ## F2S4_160412_007_H01        Sst
    ## F2S4_160412_008_A01        Sst
    ## F2S4_160412_008_B01        Sst
    ## F2S4_160412_008_C01        Sst
    ## F2S4_160412_008_D01      L5 IT
    ## F2S4_160412_008_E01        Sst
    ## F2S4_160412_008_F01        Sst
    ## F2S4_160412_008_G01        Vip
    ## F2S4_160412_008_H01        Sst
    ## F2S4_160412_009_A01        Sst
    ## F2S4_160412_009_B01        Sst
    ## F2S4_160412_009_C01        Sst
    ## F2S4_160412_009_D01        Sst
    ## F2S4_160412_009_E01      Pvalb
    ## F2S4_160412_009_F01        Sst
    ## F2S4_160412_009_G01        Sst
    ## F2S4_160412_009_H01        Vip
    ## F2S4_160412_010_A01        Sst
    ## F2S4_160412_010_B01        Sst
    ## F2S4_160412_010_C01      Pvalb
    ## F2S4_160412_010_D01      Pvalb
    ## F2S4_160412_010_E01       Sncg
    ## F2S4_160412_010_F01        Vip
    ## F2S4_160412_010_G01        Sst
    ## F2S4_160412_010_H01      Pvalb
    ## F2S4_160413_006_A01      Lamp5
    ## F2S4_160413_006_B01      Lamp5
    ## F2S4_160413_006_C01      Pvalb
    ## F2S4_160413_006_D01        Vip
    ## F2S4_160413_006_E01      Astro
    ## F2S4_160413_006_F01      Lamp5
    ## F2S4_160413_006_G01      Lamp5
    ## F2S4_160413_006_H01      Lamp5
    ## F2S4_160413_007_B01 Macrophage
    ## F2S4_160413_007_C01      Lamp5
    ## F2S4_160413_007_D01       Sncg
    ## F2S4_160413_007_E01      Lamp5
    ## F2S4_160413_007_F01 Macrophage
    ## F2S4_160413_007_G01      Oligo
    ## F2S4_160413_007_H01      Lamp5
    ## F2S4_160413_008_A01      Lamp5
    ## F2S4_160413_008_B01      Lamp5
    ## F2S4_160413_008_C01      Lamp5
    ## F2S4_160413_008_D01        Vip
    ## F2S4_160413_008_E01      Lamp5
    ## F2S4_160413_008_F01        Vip
    ## F2S4_160413_008_G01      Lamp5
    ## F2S4_160413_008_H01      Lamp5
    ## F2S4_160413_009_A01      Lamp5
    ## F2S4_160413_009_B01      Lamp5
    ## F2S4_160413_009_C01      Lamp5
    ## F2S4_160413_009_D01 Macrophage
    ## F2S4_160413_009_E01      Lamp5
    ## F2S4_160413_009_F01      Lamp5
    ## F2S4_160413_009_G01        Vip
    ## F2S4_160413_009_H01      Lamp5
    ## F2S4_160413_010_A01      Lamp5
    ## F2S4_160413_010_B01      Lamp5
    ## F2S4_160413_010_C01        Vip
    ## F2S4_160413_010_D01      Pvalb
    ## F2S4_160413_010_E01      Lamp5
    ## F2S4_160413_010_G01        Sst
    ## F2S4_160413_010_H01        Vip
    ## F2S4_160413_011_A01        Vip
    ## F2S4_160413_011_B01      Pvalb
    ## F2S4_160413_011_C01        Sst
    ## F2S4_160413_011_D01       Sncg
    ## F2S4_160413_011_E01      Pvalb
    ## F2S4_160413_011_F01      Pvalb
    ## F2S4_160413_011_G01        Sst
    ## F2S4_160413_011_H01      Pvalb
    ## F2S4_160413_012_A01        Sst
    ## F2S4_160413_012_B01      Pvalb
    ## F2S4_160413_012_C01      Lamp5
    ## F2S4_160413_012_D01      Pvalb
    ## F2S4_160413_012_E01      Pvalb
    ## F2S4_160413_012_F01        Sst
    ## F2S4_160413_012_G01      Lamp5
    ## F2S4_160413_012_H01        Sst
    ## F2S4_160413_013_A01        Sst
    ## F2S4_160413_013_B01        Sst
    ## F2S4_160413_013_C01      Lamp5
    ## F2S4_160413_013_D01       Sncg
    ## F2S4_160413_013_E01      Lamp5
    ## F2S4_160413_013_F01        Sst
    ## F2S4_160413_013_G01        Vip
    ## F2S4_160413_013_H01        Sst
    ## F2S4_160413_014_A01        Sst
    ## F2S4_160413_014_B01        Sst
    ## F2S4_160413_014_C01        Sst
    ## F2S4_160413_014_D01      Astro
    ## F2S4_160413_014_E01      Pvalb
    ## F2S4_160413_014_F01      Pvalb
    ## F2S4_160413_014_G01      Lamp5
    ## F2S4_160413_014_H01        Vip
    ## F2S4_160413_015_A01        Sst
    ## F2S4_160413_015_B01        Sst
    ## F2S4_160413_015_C01        Sst
    ## F2S4_160413_015_D01        Vip
    ## F2S4_160413_015_E01        Sst
    ## F2S4_160413_015_F01        Vip
    ## F2S4_160413_015_G01      Pvalb
    ## F2S4_160413_015_H01      Oligo
    ## F2S4_160413_016_A01        Sst
    ## F2S4_160413_016_B01        Sst
    ## F2S4_160413_016_C01        Sst
    ## F2S4_160413_016_D01        Sst
    ## F2S4_160413_016_E01        Sst
    ## F2S4_160413_016_F01        Sst
    ## F2S4_160413_016_G01        Sst
    ## F2S4_160413_016_H01       Sncg
    ## F2S4_160413_017_A01        Sst
    ## F2S4_160413_017_B01        Sst
    ## F2S4_160413_017_C01        Sst
    ## F2S4_160413_017_D01        Vip
    ## F2S4_160413_017_E01        Sst
    ## F2S4_160413_017_F01        Sst
    ## F2S4_160413_017_G01        Vip
    ## F2S4_160413_017_H01        Sst
    ## F2S4_160413_018_A01      Pvalb
    ## F2S4_160413_018_B01      Lamp5
    ## F2S4_160413_018_C01      Pvalb
    ## F2S4_160413_018_D01       Sncg
    ## F2S4_160413_018_E01      Pvalb
    ## F2S4_160413_018_F01       Sncg
    ## F2S4_160413_018_G01        Sst
    ## F2S4_160413_018_H01        Sst
    ## F2S4_160413_019_A01      Pvalb
    ## F2S4_160413_019_B01        Sst
    ## F2S4_160413_019_C01        Vip
    ## F2S4_160413_019_D01        Sst
    ## F2S4_160413_019_E01        Sst
    ## F2S4_160413_019_F01        Vip
    ## F2S4_160413_019_G01      Pvalb
    ## F2S4_160413_019_H01      Pvalb
    ## F2S4_160413_020_A01      Pvalb
    ## F2S4_160413_020_B01        Vip
    ## F2S4_160413_020_C01        Sst
    ## F2S4_160413_020_D01        Sst
    ## F2S4_160413_020_E01      Lamp5
    ## F2S4_160413_020_F01        Sst
    ## F2S4_160413_020_G01        Sst
    ## F2S4_160413_020_H01        Sst
    ## F2S4_160413_021_A01        Sst
    ## F2S4_160413_021_B01        Sst
    ## F2S4_160413_021_C01      Astro
    ## F2S4_160413_021_D01        Vip
    ## F2S4_160413_021_E01        Vip
    ## F2S4_160413_021_F01      Astro
    ## F2S4_160413_021_G01        Sst
    ## F2S4_160413_021_H01        Sst
    ## F2S4_160413_024_A01      Astro
    ## F2S4_160413_024_B01        Vip
    ## F2S4_160413_024_C01        Sst
    ## F2S4_160413_024_D01      Lamp5
    ## F2S4_160413_024_E01        Vip
    ## F2S4_160413_024_F01        Vip
    ## F2S4_160413_024_G01      Lamp5
    ## F2S4_160413_024_H01      Astro
    ## F2S4_160414_002_A01        Sst
    ## F2S4_160414_002_B01        Sst
    ## F2S4_160414_002_C01        Sst
    ## F2S4_160414_002_D01        Sst
    ## F2S4_160414_002_E01        Vip
    ## F2S4_160414_002_F01      Lamp5
    ## F2S4_160414_002_H01      Lamp5
    ## F2S4_160414_003_A01        Sst
    ## F2S4_160414_003_B01        Sst
    ## F2S4_160414_003_C01      Pvalb
    ## F2S4_160414_003_D01      Lamp5
    ## F2S4_160414_003_E01      Pvalb
    ## F2S4_160414_003_F01        Sst
    ## F2S4_160414_003_G01      Lamp5
    ## F2S4_160414_003_H01        Sst
    ## F2S4_160414_004_A01       Sncg
    ## F2S4_160414_004_B01      Pvalb
    ## F2S4_160414_004_C01        Sst
    ## F2S4_160414_004_D01      Pvalb
    ## F2S4_160414_004_E01        Sst
    ## F2S4_160414_004_F01        Sst
    ## F2S4_160414_004_G01      Pvalb
    ## F2S4_160414_004_H01        Sst
    ## F2S4_160414_005_A01        Sst
    ## F2S4_160414_005_B01        Vip
    ## F2S4_160414_005_C01      Pvalb
    ## F2S4_160414_005_D01      Pvalb
    ## F2S4_160414_005_E01        Vip
    ## F2S4_160414_005_F01        Vip
    ## F2S4_160414_005_G01        Sst
    ## F2S4_160414_005_H01      Pvalb
    ## F2S4_160414_006_A01      Lamp5
    ## F2S4_160414_006_B01        Sst
    ## F2S4_160414_006_C01      Pvalb
    ## F2S4_160414_006_D01      Lamp5
    ## F2S4_160414_006_E01        Sst
    ## F2S4_160414_006_F01      Pvalb
    ## F2S4_160414_006_G01        Sst
    ## F2S4_160414_006_H01        Sst
    ## F2S4_160414_007_A01      Lamp5
    ## F2S4_160414_007_B01      Lamp5
    ## F2S4_160414_007_C01        Sst
    ## F2S4_160414_007_D01      Lamp5
    ## F2S4_160414_007_E01        Sst
    ## F2S4_160414_007_F01        Vip
    ## F2S4_160414_007_G01        Sst
    ## F2S4_160414_007_H01        Sst
    ## F2S4_160414_008_A01      Lamp5
    ## F2S4_160414_008_C01      Lamp5
    ## F2S4_160414_008_D01        Vip
    ## F2S4_160414_008_E01        Sst
    ## F2S4_160414_008_F01        Sst
    ## F2S4_160414_008_G01        Sst
    ## F2S4_160414_008_H01        Sst
    ## F2S4_160414_009_A01      Lamp5
    ## F2S4_160414_009_B01 Macrophage
    ## F2S4_160414_009_C01      Lamp5
    ## F2S4_160414_009_D01      Lamp5
    ## F2S4_160414_009_E01      Lamp5
    ## F2S4_160414_009_F01      Lamp5
    ## F2S4_160414_009_G01      Lamp5
    ## F2S4_160414_009_H01      Lamp5
    ## F2S4_160414_010_A01      Lamp5
    ## F2S4_160414_010_B01      Lamp5
    ## F2S4_160414_010_C01        Vip
    ## F2S4_160414_010_D01      Lamp5
    ## F2S4_160414_010_E01        Vip
    ## F2S4_160414_010_F01      Lamp5
    ## F2S4_160414_010_G01      Lamp5
    ## F2S4_160414_010_H01      Lamp5
    ## F2S4_160414_011_A01        Vip
    ## F2S4_160414_011_B01        Vip
    ## F2S4_160414_011_C01        Vip
    ## F2S4_160414_011_D01      Lamp5
    ## F2S4_160414_011_E01        Vip
    ## F2S4_160414_011_F01      Lamp5
    ## F2S4_160414_011_G01        Sst
    ## F2S4_160414_011_H01        Vip
    ## F2S4_160415_001_A01      Pvalb
    ## F2S4_160415_001_C01      Pvalb
    ## F2S4_160415_001_D01        Sst
    ## F2S4_160415_001_E01        Sst
    ## F2S4_160415_001_F01        Sst
    ## F2S4_160415_001_G01      Pvalb
    ## F2S4_160415_001_H01        Vip
    ## F2S4_160415_002_A01        Sst
    ## F2S4_160415_002_B01        Vip
    ## F2S4_160415_002_C01        Sst
    ## F2S4_160415_002_D01        Sst
    ## F2S4_160415_002_E01      Pvalb
    ## F2S4_160415_002_F01        Vip
    ## F2S4_160415_002_G01        Sst
    ## F2S4_160415_002_H01      Lamp5
    ## F2S4_160415_003_A01        Vip
    ## F2S4_160415_003_B01      Pvalb
    ## F2S4_160415_003_C01        Sst
    ## F2S4_160415_003_D01        Vip
    ## F2S4_160415_003_E01        Sst
    ## F2S4_160415_003_F01        Sst
    ## F2S4_160415_003_H01        Sst
    ## F2S4_160415_004_A01        Sst
    ## F2S4_160415_004_B01      Pvalb
    ## F2S4_160415_004_C01   Serpinf1
    ## F2S4_160415_004_D01        Sst
    ## F2S4_160415_004_E01        Sst
    ## F2S4_160415_004_F01        Sst
    ## F2S4_160415_004_G01        Vip
    ## F2S4_160415_004_H01        Sst
    ## F2S4_160415_005_A01      Pvalb
    ## F2S4_160415_005_B01        Sst
    ## F2S4_160415_005_C01      Lamp5
    ## F2S4_160415_005_D01        Sst
    ## F2S4_160415_005_E01        Sst
    ## F2S4_160415_005_G01      Lamp5
    ## F2S4_160415_005_H01        Sst
    ## F2S4_160415_006_A01        Sst
    ## F2S4_160415_006_B01        Sst
    ## F2S4_160415_006_C01        Sst
    ## F2S4_160415_006_D01      Pvalb
    ## F2S4_160415_006_E01        Sst
    ## F2S4_160415_006_F01      Pvalb
    ## F2S4_160415_006_G01      Pvalb
    ## F2S4_160415_006_H01      Pvalb
    ## F2S4_160415_007_A01        Vip
    ## F2S4_160415_007_B01        Sst
    ## F2S4_160415_007_C01        Sst
    ## F2S4_160415_007_D01      Pvalb
    ## F2S4_160415_007_E01        Sst
    ## F2S4_160415_007_F01      Lamp5
    ## F2S4_160415_007_G01        Sst
    ## F2S4_160415_007_H01        Sst
    ## F2S4_160415_008_A01        Sst
    ## F2S4_160415_008_B01      Pvalb
    ## F2S4_160415_008_C01        Sst
    ## F2S4_160415_008_D01        Sst
    ## F2S4_160415_008_E01        Sst
    ## F2S4_160415_008_F01        Sst
    ## F2S4_160415_008_G01      Pvalb
    ## F2S4_160415_008_H01      Pvalb
    ## F2S4_160415_009_A01        Sst
    ## F2S4_160415_009_B01        Sst
    ## F2S4_160415_009_C01      Pvalb
    ## F2S4_160415_009_D01        Sst
    ## F2S4_160415_009_E01   Serpinf1
    ## F2S4_160415_009_F01      Lamp5
    ## F2S4_160415_009_G01      Pvalb
    ## F2S4_160415_009_H01        Sst
    ## F2S4_160415_010_A01      Pvalb
    ## F2S4_160415_010_B01        Sst
    ## F2S4_160415_010_C01        Sst
    ## F2S4_160415_010_D01        Sst
    ## F2S4_160415_010_E01        Sst
    ## F2S4_160415_010_F01        Vip
    ## F2S4_160415_010_G01        Sst
    ## F2S4_160415_010_H01        Vip
    ## F2S4_160415_011_A01      Pvalb
    ## F2S4_160415_011_B01      Meis2
    ## F2S4_160415_011_C01        Sst
    ## F2S4_160415_011_D01        Sst
    ## F2S4_160415_011_E01      Lamp5
    ## F2S4_160415_011_F01        Vip
    ## F2S4_160415_011_G01        Sst
    ## F2S4_160415_011_H01       Sncg
    ## F2S4_160415_012_A01        Sst
    ## F2S4_160415_012_C01      Lamp5
    ## F2S4_160415_012_E01   Serpinf1
    ## F2S4_160415_012_F01      Pvalb
    ## F2S4_160415_012_G01        Sst
    ## F2S4_160415_012_H01        Sst
    ## F2S4_160415_013_A01      Lamp5
    ## F2S4_160415_013_B01        Vip
    ## F2S4_160415_013_C01        Sst
    ## F2S4_160415_013_D01        Vip
    ## F2S4_160415_013_E01      Pvalb
    ## F2S4_160415_013_F01        Sst
    ## F2S4_160415_013_G01        Sst
    ## F2S4_160415_013_H01      Lamp5
    ## F2S4_160415_014_A01        Vip
    ## F2S4_160415_014_B01      Pvalb
    ## F2S4_160415_014_C01      Lamp5
    ## F2S4_160415_014_D01      Pvalb
    ## F2S4_160415_014_E01      Lamp5
    ## F2S4_160415_014_F01        Sst
    ## F2S4_160415_014_G01        Sst
    ## F2S4_160415_014_H01        Vip
    ## F2S4_160415_015_A01      Pvalb
    ## F2S4_160415_015_B01      Pvalb
    ## F2S4_160415_015_C01        Sst
    ## F2S4_160415_015_D01        Sst
    ## F2S4_160415_015_E01        Vip
    ## F2S4_160415_015_F01      Lamp5
    ## F2S4_160415_015_G01        Vip
    ## F2S4_160415_015_H01      Lamp5
    ## F2S4_160415_016_A01      Lamp5
    ## F2S4_160415_016_B01      Lamp5
    ## F2S4_160415_016_C01      Lamp5
    ## F2S4_160415_016_D01        Vip
    ## F2S4_160415_016_E01      Lamp5
    ## F2S4_160415_016_F01      Lamp5
    ## F2S4_160415_016_G01        Vip
    ## F2S4_160415_016_H01        Vip
    ## F2S4_160415_017_A01      Lamp5
    ## F2S4_160415_017_B01        Vip
    ## F2S4_160415_017_C01      Lamp5
    ## F2S4_160415_017_D01        Vip
    ## F2S4_160415_017_E01        Vip
    ## F2S4_160415_017_F01      Lamp5
    ## F2S4_160415_017_G01      Lamp5
    ## F2S4_160415_017_H01       Sncg
    ## F2S4_160415_018_A01      Lamp5
    ## F2S4_160415_018_B01      Lamp5
    ## F2S4_160415_018_C01      Lamp5
    ## F2S4_160415_018_D01       Sncg
    ## F2S4_160415_018_E01      Lamp5
    ## F2S4_160415_018_F01        Vip
    ## F2S4_160415_018_G01      Lamp5
    ## F2S4_160415_018_H01      Lamp5
    ## F2S4_160415_019_A01        Vip
    ## F2S4_160415_019_B01      Lamp5
    ## F2S4_160415_019_C01      Lamp5
    ## F2S4_160415_019_D01      Lamp5
    ## F2S4_160415_019_E01        Vip
    ## F2S4_160415_019_F01      Pvalb
    ## F2S4_160415_019_G01      Lamp5
    ## F2S4_160415_019_H01      Lamp5
    ## F2S4_160415_020_A01      Lamp5
    ## F2S4_160415_020_B01      Lamp5
    ## F2S4_160415_020_C01      Lamp5
    ## F2S4_160415_020_D01      Lamp5
    ## F2S4_160415_020_E01        Vip
    ## F2S4_160415_020_F01        Vip
    ## F2S4_160415_020_G01      Lamp5
    ## F2S4_160415_020_H01       Sncg
    ## F2S4_160415_021_A01      Lamp5
    ## F2S4_160415_021_B01        Vip
    ## F2S4_160415_021_C01      Lamp5
    ## F2S4_160415_021_D01      Lamp5
    ## F2S4_160415_021_E01      Lamp5
    ## F2S4_160415_021_F01      Lamp5
    ## F2S4_160415_021_G01      Lamp5
    ## F2S4_160415_021_H01      Lamp5
    ## F2S4_160415_022_A01      Lamp5
    ## F2S4_160415_022_B01      Lamp5
    ## F2S4_160415_022_C01      Lamp5
    ## F2S4_160415_022_D01        Vip
    ## F2S4_160415_022_E01      Lamp5
    ## F2S4_160415_022_F01      Lamp5
    ## F2S4_160415_022_G01      Lamp5
    ## F2S4_160415_022_H01      Lamp5
    ## F2S4_160415_023_A01        Vip
    ## F2S4_160415_023_B01        Vip
    ## F2S4_160415_023_C01   Serpinf1
    ## F2S4_160415_023_D01      Lamp5
    ## F2S4_160415_023_E01      Lamp5
    ## F2S4_160415_023_F01      Lamp5
    ## F2S4_160415_023_G01      Lamp5
    ## F2S4_160415_023_H01      Lamp5
    ## F2S4_160415_024_A01        Vip
    ## F2S4_160415_024_B01        Vip
    ## F2S4_160415_024_C01        Vip
    ## F2S4_160415_024_D01        Vip
    ## F2S4_160415_024_E01        Vip
    ## F2S4_160415_024_F01        Vip
    ## F2S4_160415_024_G01        Vip
    ## F2S4_160415_024_H01        Vip
    ## F2S4_160415_025_A01        Sst
    ## F2S4_160415_025_B01        Vip
    ## F2S4_160415_025_C01        Vip
    ## F2S4_160415_025_D01        Vip
    ## F2S4_160415_025_E01        Vip
    ## F2S4_160415_025_F01        Vip
    ## F2S4_160415_025_G01        Vip
    ## F2S4_160415_025_H01        Vip
    ## F2S4_160419_001_A01      Astro
    ## F2S4_160419_001_B01      Astro
    ## F2S4_160419_001_C01      Lamp5
    ## F2S4_160419_001_D01      Lamp5
    ## F2S4_160419_001_E01      Lamp5
    ## F2S4_160419_001_F01      Lamp5
    ## F2S4_160419_001_G01      Lamp5
    ## F2S4_160419_001_H01      Lamp5
    ## F2S4_160419_002_A01        Sst
    ## F2S4_160419_002_B01      Lamp5
    ## F2S4_160419_002_C01        Vip
    ## F2S4_160419_002_D01      Lamp5
    ## F2S4_160419_002_E01        Vip
    ## F2S4_160419_002_F01      Lamp5
    ## F2S4_160419_002_G01      Lamp5
    ## F2S4_160419_002_H01      Lamp5
    ## F2S4_160419_003_A01        Vip
    ## F2S4_160419_003_B01      Lamp5
    ## F2S4_160419_003_C01        Vip
    ## F2S4_160419_003_D01      Lamp5
    ## F2S4_160419_003_E01        Vip
    ## F2S4_160419_003_F01      Lamp5
    ## F2S4_160419_003_G01      Lamp5
    ## F2S4_160419_003_H01        Sst
    ## F2S4_160419_004_A01        Sst
    ## F2S4_160419_004_B01      Astro
    ## F2S4_160419_004_C01        Sst
    ## F2S4_160419_004_D01        Vip
    ## F2S4_160419_004_E01        Vip
    ## F2S4_160419_004_F01        Sst
    ## F2S4_160419_004_G01      Pvalb
    ## F2S4_160419_004_H01      Lamp5
    ## F2S4_160419_005_A01      Pvalb
    ## F2S4_160419_005_B01        Vip
    ## F2S4_160419_005_D01        Vip
    ## F2S4_160419_005_E01        Vip
    ## F2S4_160419_005_F01        Sst
    ## F2S4_160419_005_G01        Sst
    ## F2S4_160419_005_H01        Vip
    ## F2S4_160419_006_A01      Oligo
    ## F2S4_160419_006_B01        Vip
    ## F2S4_160419_006_C01        Vip
    ## F2S4_160419_006_D01        Vip
    ## F2S4_160419_006_E01        Vip
    ## F2S4_160419_006_F01      Pvalb
    ## F2S4_160419_006_G01        Vip
    ## F2S4_160419_006_H01      Astro
    ## F2S4_160419_007_A01        Vip
    ## F2S4_160419_007_C01        Vip
    ## F2S4_160419_007_D01        Vip
    ## F2S4_160419_007_E01        Sst
    ## F2S4_160419_007_F01        Sst
    ## F2S4_160419_007_G01        Vip
    ## F2S4_160419_007_H01        Vip
    ## F2S4_160419_008_A01      Pvalb
    ## F2S4_160419_008_B01        Sst
    ## F2S4_160419_008_C01       Sncg
    ## F2S4_160419_008_D01      Pvalb
    ## F2S4_160419_008_E01        Sst
    ## F2S4_160419_008_F01   Serpinf1
    ## F2S4_160419_008_G01        Sst
    ## F2S4_160419_008_H01      Meis2
    ## F2S4_160419_009_A01        Sst
    ## F2S4_160419_009_B01        Sst
    ## F2S4_160419_009_C01      Pvalb
    ## F2S4_160419_009_D01       Sncg
    ## F2S4_160419_009_E01        Vip
    ## F2S4_160419_009_F01        Sst
    ## F2S4_160419_009_G01        Sst
    ## F2S4_160419_009_H01        Sst
    ## F2S4_160419_010_A01        Sst
    ## F2S4_160419_010_B01        Sst
    ## F2S4_160419_010_C01      Lamp5
    ## F2S4_160419_010_D01        Sst
    ## F2S4_160419_010_E01        Sst
    ## F2S4_160419_010_F01        Sst
    ## F2S4_160419_010_G01        Sst
    ## F2S4_160419_010_H01        Sst
    ## F2S4_160419_011_A01        Vip
    ## F2S4_160419_011_B01        Vip
    ## F2S4_160419_011_C01        Vip
    ## F2S4_160419_011_D01        Vip
    ## F2S4_160419_011_E01        Vip
    ## F2S4_160419_011_F01      Lamp5
    ## F2S4_160419_011_G01        Vip
    ## F2S4_160419_011_H01      Lamp5
    ## F2S4_160419_012_A01    L2/3 IT
    ## F2S4_160419_012_B01      Pvalb
    ## F2S4_160419_012_C01        Vip
    ## F2S4_160419_012_D01        Vip
    ## F2S4_160419_012_E01      Lamp5
    ## F2S4_160419_012_F01        Vip
    ## F2S4_160419_012_G01        Vip
    ## F2S4_160419_012_H01        Vip
    ## F2S4_160419_013_A01    L2/3 IT
    ## F2S4_160419_013_B01        Vip
    ## F2S4_160419_013_C01       Sncg
    ## F2S4_160419_013_D01        Vip
    ## F2S4_160419_013_E01        Vip
    ## F2S4_160419_013_F01        Sst
    ## F2S4_160419_013_G01        Vip
    ## F2S4_160419_013_H01        Vip
    ## F2S4_160419_014_B01      Lamp5
    ## F2S4_160419_014_C01        Vip
    ## F2S4_160419_014_D01      Lamp5
    ## F2S4_160419_014_F01        Sst
    ## F2S4_160419_014_G01        Sst
    ## F2S4_160419_014_H01        Vip
    ## F2S4_160419_015_A01      Pvalb
    ## F2S4_160419_015_B01        Sst
    ## F2S4_160419_015_C01        Sst
    ## F2S4_160419_015_D01        Sst
    ## F2S4_160419_015_E01        Vip
    ## F2S4_160419_015_F01      Lamp5
    ## F2S4_160419_015_G01      Pvalb
    ## F2S4_160419_015_H01   Serpinf1
    ## F2S4_160419_016_A01        Vip
    ## F2S4_160419_016_B01        Sst
    ## F2S4_160419_016_C01      Pvalb
    ## F2S4_160419_016_D01        Vip
    ## F2S4_160419_016_E01        Sst
    ## F2S4_160419_016_F01        Vip
    ## F2S4_160419_016_G01      Lamp5
    ## F2S4_160419_016_H01      Pvalb
    ## F2S4_160419_017_A01      Astro
    ## F2S4_160419_017_B01      Astro
    ## F2S4_160419_017_C01      Pvalb
    ## F2S4_160419_017_D01        Vip
    ## F2S4_160419_017_E01      Astro
    ## F2S4_160419_017_F01        Vip
    ## F2S4_160419_017_G01      Pvalb
    ## F2S4_160419_017_H01        Vip
    ## F2S4_160419_018_B01        Sst
    ## F2S4_160419_018_C01        Vip
    ## F2S4_160419_018_E01        Sst
    ## F2S4_160419_018_F01        Sst
    ## F2S4_160419_018_G01      Pvalb
    ## F2S4_160419_018_H01        Vip
    ## F2S4_160419_019_A01      Pvalb
    ## F2S4_160419_019_B01        Vip
    ## F2S4_160419_019_C01        Sst
    ## F2S4_160419_019_D01        Sst
    ## F2S4_160419_019_E01        Sst
    ## F2S4_160419_019_F01        Sst
    ## F2S4_160419_019_G01        Sst
    ## F2S4_160419_019_H01      Pvalb
    ## F2S4_160419_020_A01        Sst
    ## F2S4_160419_020_B01      Pvalb
    ## F2S4_160419_020_C01        Sst
    ## F2S4_160419_020_D01        Sst
    ## F2S4_160419_020_E01      Pvalb
    ## F2S4_160419_020_F01      Pvalb
    ## F2S4_160419_020_G01        Sst
    ## F2S4_160419_020_H01      Pvalb
    ## F2S4_160419_021_A01      Astro
    ## F2S4_160419_021_B01        Vip
    ## F2S4_160419_021_C01        Sst
    ## F2S4_160419_021_D01      Astro
    ## F2S4_160419_021_E01        Vip
    ## F2S4_160419_021_F01   Serpinf1
    ## F2S4_160419_021_G01      Astro
    ## F2S4_160419_021_H01        Vip
    ## F2S4_160419_022_A01      Lamp5
    ## F2S4_160419_022_B01      Astro
    ## F2S4_160419_022_C01       Sncg
    ## F2S4_160419_022_D01        Vip
    ## F2S4_160419_022_E01        Vip
    ## F2S4_160419_022_G01      Astro
    ## F2S4_160419_022_H01        Vip
    ## F2S4_160419_023_A01 Macrophage
    ## F2S4_160419_023_B01        Vip
    ## F2S4_160419_023_C01        Vip
    ## F2S4_160419_023_D01        Vip
    ## F2S4_160419_023_E01      Lamp5
    ## F2S4_160419_023_F01        Vip
    ## F2S4_160419_023_G01        Vip
    ## F2S4_160419_023_H01        Vip
    ## F2S4_160422_001_A01      L5 IT
    ## F2S4_160422_001_B01         L4
    ## F2S4_160422_001_C01      L5 IT
    ## F2S4_160422_001_D01         L4
    ## F2S4_160422_001_E01         L4
    ## F2S4_160422_001_F01         L4
    ## F2S4_160422_001_G01         L4
    ## F2S4_160422_001_H01         L4
    ## F2S4_160422_002_B01        Sst
    ## F2S4_160422_002_C01         L4
    ## F2S4_160422_002_D01         L4
    ## F2S4_160422_002_E01         L4
    ## F2S4_160422_002_F01         L4
    ## F2S4_160422_002_G01        Vip
    ## F2S4_160422_002_H01         L4
    ## F2S4_160422_003_A01         L4
    ## F2S4_160422_003_B01        Vip
    ## F2S4_160422_003_C01         L4
    ## F2S4_160422_003_D01        Vip
    ## F2S4_160422_003_E01         L4
    ## F2S4_160422_003_F01         L4
    ## F2S4_160422_003_G01         L4
    ## F2S4_160422_003_H01         L4
    ## F2S4_160422_004_A01         L4
    ## F2S4_160422_004_B01      L5 IT
    ## F2S4_160422_004_C01         L4
    ## F2S4_160422_004_D01         L4
    ## F2S4_160422_004_E01         L4
    ## F2S4_160422_004_F01         L4
    ## F2S4_160422_004_G01         L4
    ## F2S4_160422_004_H01         L4
    ## F2S4_160422_005_A01         L4
    ## F2S4_160422_005_B01         L4
    ## F2S4_160422_005_C01         L4
    ## F2S4_160422_005_D01         L4
    ## F2S4_160422_005_E01         L4
    ## F2S4_160422_005_F01         L4
    ## F2S4_160422_005_G01         L4
    ## F2S4_160422_005_H01         L4
    ## F2S4_160422_006_A01         L4
    ## F2S4_160422_006_B01         L4
    ## F2S4_160422_006_C01         L4
    ## F2S4_160422_006_D01         L4
    ## F2S4_160422_006_F01         L4
    ## F2S4_160422_006_G01         L4
    ## F2S4_160422_006_H01         L4
    ## F2S4_160422_007_A01    L2/3 IT
    ## F2S4_160422_007_B01    L2/3 IT
    ## F2S4_160422_007_C01    L2/3 IT
    ## F2S4_160422_007_D01    L2/3 IT
    ## F2S4_160422_007_E01    L2/3 IT
    ## F2S4_160422_007_F01    L2/3 IT
    ## F2S4_160422_007_G01        Vip
    ## F2S4_160422_007_H01    L2/3 IT
    ## F2S4_160422_008_A01    L2/3 IT
    ## F2S4_160422_008_B01    L2/3 IT
    ## F2S4_160422_008_C01    L2/3 IT
    ## F2S4_160422_008_D01    L2/3 IT
    ## F2S4_160422_008_E01    L2/3 IT
    ## F2S4_160422_008_F01    L2/3 IT
    ## F2S4_160422_008_G01      Lamp5
    ## F2S4_160422_008_H01    L2/3 IT
    ## F2S4_160422_009_A01    L2/3 IT
    ## F2S4_160422_009_B01      Lamp5
    ## F2S4_160422_009_C01        Vip
    ## F2S4_160422_009_D01    L2/3 IT
    ## F2S4_160422_009_E01    L2/3 IT
    ## F2S4_160422_009_F01    L2/3 IT
    ## F2S4_160422_009_G01    L2/3 IT
    ## F2S4_160422_009_H01        Vip
    ## F2S4_160422_010_A01      Lamp5
    ## F2S4_160422_010_B01        Vip
    ## F2S4_160422_010_C01    L2/3 IT
    ## F2S4_160422_010_D01      Lamp5
    ## F2S4_160422_010_E01    L2/3 IT
    ## F2S4_160422_010_F01      Lamp5
    ## F2S4_160422_010_G01    L2/3 IT
    ## F2S4_160422_010_H01    L2/3 IT
    ## F2S4_160422_011_A01    L2/3 IT
    ## F2S4_160422_011_B01    L2/3 IT
    ## F2S4_160422_011_C01    L2/3 IT
    ## F2S4_160422_011_D01    L2/3 IT
    ## F2S4_160422_011_E01    L2/3 IT
    ## F2S4_160422_011_F01        Vip
    ## F2S4_160422_011_G01    L2/3 IT
    ## F2S4_160422_011_H01    L2/3 IT
    ## F2S4_160422_012_A01    L2/3 IT
    ## F2S4_160422_012_B01    L2/3 IT
    ## F2S4_160422_012_D01    L2/3 IT
    ## F2S4_160422_012_E01        Vip
    ## F2S4_160422_012_F01    L2/3 IT
    ## F2S4_160422_012_G01        Vip
    ## F2S4_160422_012_H01    L2/3 IT
    ## F2S4_160422_013_A01      L6 CT
    ## F2S4_160422_013_B01      L6 IT
    ## F2S4_160422_013_C01      L6 CT
    ## F2S4_160422_013_D01      L6 IT
    ## F2S4_160422_013_E01      L6 CT
    ## F2S4_160422_013_F01      L6 IT
    ## F2S4_160422_013_G01      L6 IT
    ## F2S4_160422_013_H01      L6 IT
    ## F2S4_160422_014_A01      L6 CT
    ## F2S4_160422_014_B01      L6 IT
    ## F2S4_160422_014_C01      L6 IT
    ## F2S4_160422_014_D01      L6 IT
    ## F2S4_160422_014_E01      L6 IT
    ## F2S4_160422_014_F01      L6 IT
    ## F2S4_160422_014_G01      L6 IT
    ## F2S4_160422_014_H01      L6 IT
    ## F2S4_160422_015_A01      L6 IT
    ## F2S4_160422_015_B01      L6 CT
    ## F2S4_160422_015_C01      L6 IT
    ## F2S4_160422_015_D01      L6 IT
    ## F2S4_160422_015_E01      L6 IT
    ## F2S4_160422_015_F01      L6 IT
    ## F2S4_160422_015_G01      L6 IT
    ## F2S4_160422_015_H01        Sst
    ## F2S4_160422_016_A01      L6 IT
    ## F2S4_160422_016_B01      L6 IT
    ## F2S4_160422_016_C01        Vip
    ## F2S4_160422_016_D01      L6 IT
    ## F2S4_160422_016_E01        Sst
    ## F2S4_160422_016_F01      L6 IT
    ## F2S4_160422_016_G01      L6 IT
    ## F2S4_160422_016_H01        Sst
    ## F2S4_160422_017_A01      L6 IT
    ## F2S4_160422_017_B01      L6 IT
    ## F2S4_160422_017_C01      L6 CT
    ## F2S4_160422_017_D01      L6 IT
    ## F2S4_160422_017_E01      L6 IT
    ## F2S4_160422_017_F01      L6 IT
    ## F2S4_160422_017_G01      L6 IT
    ## F2S4_160422_017_H01      L6 IT
    ## F2S4_160422_018_B01       Sncg
    ## F2S4_160422_018_C01      L6 CT
    ## F2S4_160422_018_D01      L6 IT
    ## F2S4_160422_018_E01      L6 IT
    ## F2S4_160422_018_F01      L6 IT
    ## F2S4_160422_018_G01      L6 CT
    ## F2S4_160422_018_H01      L6 CT
    ## F2S4_160422_019_A01      L5 IT
    ## F2S4_160422_019_B01      L5 IT
    ## F2S4_160422_019_C01        Vip
    ## F2S4_160422_019_D01         NP
    ## F2S4_160422_019_E01      L6 IT
    ## F2S4_160422_019_F01      L6 IT
    ## F2S4_160422_019_G01        Vip
    ## F2S4_160422_020_A01      L6 IT
    ## F2S4_160422_020_B01      L6 IT
    ## F2S4_160422_020_C01      L6 IT
    ## F2S4_160422_020_D01        Sst
    ## F2S4_160422_020_E01        Vip
    ## F2S4_160422_020_F01         NP
    ## F2S4_160422_020_G01        Vip
    ## F2S4_160422_020_H01      L5 PT
    ## F2S4_160422_021_A01      L5 IT
    ## F2S4_160422_021_B01         NP
    ## F2S4_160422_021_C01       Sncg
    ## F2S4_160422_021_D01      L6 IT
    ## F2S4_160422_021_E01         NP
    ## F2S4_160422_021_F01         NP
    ## F2S4_160422_021_G01         NP
    ## F2S4_160422_021_H01      L5 IT
    ## F2S4_160422_022_A01        Vip
    ## F2S4_160422_022_B01         NP
    ## F2S4_160422_022_C01        Sst
    ## F2S4_160422_022_D01      L5 IT
    ## F2S4_160422_022_E01         NP
    ## F2S4_160422_022_F01        Vip
    ## F2S4_160422_022_G01      L6 IT
    ## F2S4_160422_022_H01      L6 IT
    ## F2S4_160422_023_A01         NP
    ## F2S4_160422_023_B01         NP
    ## F2S4_160422_023_C01      L6 IT
    ## F2S4_160422_023_D01         NP
    ## F2S4_160422_023_E01      L6 IT
    ## F2S4_160422_023_F01        Sst
    ## F2S4_160422_023_G01      L5 PT
    ## F2S4_160422_023_H01      L6 CT
    ## F2S4_160422_024_A01         NP
    ## F2S4_160422_024_B01        Vip
    ## F2S4_160422_024_C01        Vip
    ## F2S4_160422_024_D01         NP
    ## F2S4_160422_024_E01      L6 IT
    ## F2S4_160422_024_F01      L6 IT
    ## F2S4_160422_024_G01         NP
    ## F2S4_160422_024_H01      L6 IT
    ## F2S4_160426_001_A01      L6 IT
    ## F2S4_160426_001_B01      L6 IT
    ## F2S4_160426_001_C01      L6 IT
    ## F2S4_160426_001_D01      L6 IT
    ## F2S4_160426_001_E01      L5 IT
    ## F2S4_160426_001_F01      L6 IT
    ## F2S4_160426_001_G01      L5 IT
    ## F2S4_160426_001_H01         NP
    ## F2S4_160426_002_A01      L5 PT
    ## F2S4_160426_002_B01      L5 PT
    ## F2S4_160426_002_D01         NP
    ## F2S4_160426_002_G01      L6 IT
    ## F2S4_160426_002_H01      L5 IT
    ## F2S4_160426_003_A01      L6 IT
    ## F2S4_160426_003_B01      L6 IT
    ## F2S4_160426_003_C01      L5 PT
    ## F2S4_160426_003_D01      L6 IT
    ## F2S4_160426_003_E01      L6 IT
    ## F2S4_160426_003_F01      L5 IT
    ## F2S4_160426_003_G01         NP
    ## F2S4_160426_003_H01      L5 PT
    ## F2S4_160426_004_A01      L6 IT
    ## F2S4_160426_004_B01      L6 IT
    ## F2S4_160426_004_C01      L5 PT
    ## F2S4_160426_004_D01      L5 PT
    ## F2S4_160426_004_E01      L6 IT
    ## F2S4_160426_004_F01      L6 IT
    ## F2S4_160426_004_G01         NP
    ## F2S4_160426_004_H01      L6 IT
    ## F2S4_160426_005_A01      L6 IT
    ## F2S4_160426_005_B01      L6 IT
    ## F2S4_160426_005_C01         NP
    ## F2S4_160426_005_D01         NP
    ## F2S4_160426_005_E01      L5 IT
    ## F2S4_160426_005_F01      L6 IT
    ## F2S4_160426_005_G01      L6 IT
    ## F2S4_160426_005_H01      L6 IT
    ## F2S4_160426_006_A01      L6 IT
    ## F2S4_160426_006_B01      L5 PT
    ## F2S4_160426_006_C01      L6 IT
    ## F2S4_160426_006_D01      L6 IT
    ## F2S4_160426_006_E01      L6 IT
    ## F2S4_160426_006_F01      L5 IT
    ## F2S4_160426_006_G01      L6 IT
    ## F2S4_160426_006_H01      L6 IT
    ## F2S4_160426_007_A01      L6 IT
    ## F2S4_160426_007_B01      L6 IT
    ## F2S4_160426_007_C01      L6 CT
    ## F2S4_160426_007_D01      L6 IT
    ## F2S4_160426_007_E01      L6 IT
    ## F2S4_160426_007_F01      L6 IT
    ## F2S4_160426_007_G01      L6 CT
    ## F2S4_160426_007_H01      L6 IT
    ## F2S4_160426_008_A01      L6 IT
    ## F2S4_160426_008_B01      L6 IT
    ## F2S4_160426_008_C01      L6 IT
    ## F2S4_160426_008_D01      L6 IT
    ## F2S4_160426_008_E01      L6 IT
    ## F2S4_160426_008_F01        L6b
    ## F2S4_160426_008_G01      L6 CT
    ## F2S4_160426_008_H01        L6b
    ## F2S4_160426_009_A01      L6 IT
    ## F2S4_160426_009_B01      L6 IT
    ## F2S4_160426_009_C01      L6 IT
    ## F2S4_160426_009_D01      L6 IT
    ## F2S4_160426_009_E01      L6 CT
    ## F2S4_160426_009_F01      L6 IT
    ## F2S4_160426_009_G01      L6 IT
    ## F2S4_160426_009_H01      L6 IT
    ## F2S4_160426_010_A01      L6 IT
    ## F2S4_160426_010_B01      L6 IT
    ## F2S4_160426_010_C01      L6 CT
    ## F2S4_160426_010_D01      L6 IT
    ## F2S4_160426_010_E01      L6 IT
    ## F2S4_160426_010_F01      L6 CT
    ## F2S4_160426_010_G01      L6 IT
    ## F2S4_160426_010_H01      L6 IT
    ## F2S4_160426_011_A01      L6 IT
    ## F2S4_160426_011_B01      L6 IT
    ## F2S4_160426_011_C01      L6 IT
    ## F2S4_160426_011_D01      L6 IT
    ## F2S4_160426_011_E01      L6 IT
    ## F2S4_160426_011_F01      L6 CT
    ## F2S4_160426_011_G01      L6 CT
    ## F2S4_160426_011_H01      L6 IT
    ## F2S4_160426_012_A01      L6 CT
    ## F2S4_160426_012_B01      L6 IT
    ## F2S4_160426_012_D01      L6 IT
    ## F2S4_160426_012_E01      L6 IT
    ## F2S4_160426_012_F01      L6 IT
    ## F2S4_160426_012_G01      L6 IT
    ## F2S4_160426_013_A01    L2/3 IT
    ## F2S4_160426_013_B01    L2/3 IT
    ## F2S4_160426_013_C01    L2/3 IT
    ## F2S4_160426_013_D01    L2/3 IT
    ## F2S4_160426_013_E01    L2/3 IT
    ## F2S4_160426_013_F01    L2/3 IT
    ## F2S4_160426_013_G01    L2/3 IT
    ## F2S4_160426_013_H01         L4
    ## F2S4_160426_014_A01    L2/3 IT
    ## F2S4_160426_014_B01    L2/3 IT
    ## F2S4_160426_014_C01    L2/3 IT
    ## F2S4_160426_014_D01    L2/3 IT
    ## F2S4_160426_014_E01    L2/3 IT
    ## F2S4_160426_014_F01    L2/3 IT
    ## F2S4_160426_014_G01    L2/3 IT
    ## F2S4_160426_014_H01    L2/3 IT
    ## F2S4_160426_015_A01    L2/3 IT
    ## F2S4_160426_015_B01    L2/3 IT
    ## F2S4_160426_015_C01    L2/3 IT
    ## F2S4_160426_015_D01    L2/3 IT
    ## F2S4_160426_015_E01    L2/3 IT
    ## F2S4_160426_015_F01    L2/3 IT
    ## F2S4_160426_015_G01    L2/3 IT
    ## F2S4_160426_015_H01    L2/3 IT
    ## F2S4_160426_016_A01    L2/3 IT
    ## F2S4_160426_016_B01    L2/3 IT
    ## F2S4_160426_016_C01    L2/3 IT
    ## F2S4_160426_016_D01    L2/3 IT
    ## F2S4_160426_016_E01    L2/3 IT
    ## F2S4_160426_016_F01    L2/3 IT
    ## F2S4_160426_016_G01    L2/3 IT
    ## F2S4_160426_016_H01    L2/3 IT
    ## F2S4_160426_017_A01    L2/3 IT
    ## F2S4_160426_017_B01    L2/3 IT
    ## F2S4_160426_017_C01    L2/3 IT
    ## F2S4_160426_017_D01    L2/3 IT
    ## F2S4_160426_017_E01    L2/3 IT
    ## F2S4_160426_017_F01    L2/3 IT
    ## F2S4_160426_017_G01    L2/3 IT
    ## F2S4_160426_017_H01    L2/3 IT
    ## F2S4_160426_018_A01         L4
    ## F2S4_160426_018_B01    L2/3 IT
    ## F2S4_160426_018_C01    L2/3 IT
    ## F2S4_160426_018_D01    L2/3 IT
    ## F2S4_160426_018_E01    L2/3 IT
    ## F2S4_160426_018_F01    L2/3 IT
    ## F2S4_160426_018_G01    L2/3 IT
    ## F2S4_160426_018_H01    L2/3 IT
    ## F2S4_160426_019_A01         L4
    ## F2S4_160426_019_B01         L4
    ## F2S4_160426_019_C01         L4
    ## F2S4_160426_019_D01         L4
    ## F2S4_160426_019_E01         L4
    ## F2S4_160426_019_F01         L4
    ## F2S4_160426_019_G01         L4
    ## F2S4_160426_019_H01      L5 IT
    ## F2S4_160426_020_A01         L4
    ## F2S4_160426_020_B01         L4
    ## F2S4_160426_020_C01         L4
    ## F2S4_160426_020_D01         L4
    ## F2S4_160426_020_E01         L4
    ## F2S4_160426_020_F01         L4
    ## F2S4_160426_020_G01         L4
    ## F2S4_160426_021_A01         L4
    ## F2S4_160426_021_B01         L4
    ## F2S4_160426_021_C01         L4
    ## F2S4_160426_021_D01         L4
    ## F2S4_160426_021_E01         L4
    ## F2S4_160426_021_F01         L4
    ## F2S4_160426_021_G01         L4
    ## F2S4_160426_021_H01         L4
    ## F2S4_160426_022_A01         L4
    ## F2S4_160426_022_B01         L4
    ## F2S4_160426_022_C01      L5 IT
    ## F2S4_160426_022_E01         L4
    ## F2S4_160426_022_F01         L4
    ## F2S4_160426_022_G01         L4
    ## F2S4_160426_023_B01         L4
    ## F2S4_160426_023_C01         L4
    ## F2S4_160426_023_E01      L5 IT
    ## F2S4_160426_023_F01         L4
    ## F2S4_160426_023_G01         L4
    ## F2S4_160426_023_H01         L4
    ## F2S4_160426_024_A01         L4
    ## F2S4_160426_024_C01         L4
    ## F2S4_160426_024_D01         L4
    ## F2S4_160426_024_E01         L4
    ## F2S4_160426_024_F01         L4
    ## F2S4_160426_024_G01         L4
    ## F2S4_160426_024_H01         L4
    ## F2S4_160428_001_A01         NP
    ## F2S4_160428_001_B01         NP
    ## F2S4_160428_001_C01         NP
    ## F2S4_160428_001_D01      L5 IT
    ## F2S4_160428_001_E01         NP
    ## F2S4_160428_001_F01      L6 IT
    ## F2S4_160428_001_H01      L5 IT
    ## F2S4_160428_002_A01         NP
    ## F2S4_160428_002_B01         NP
    ## F2S4_160428_002_C01      L5 PT
    ## F2S4_160428_002_D01         NP
    ## F2S4_160428_002_E01      L6 IT
    ## F2S4_160428_002_F01         NP
    ## F2S4_160428_002_G01         NP
    ## F2S4_160428_002_H01      L5 IT
    ## F2S4_160428_003_B01         NP
    ## F2S4_160428_003_C01         NP
    ## F2S4_160428_003_D01         NP
    ## F2S4_160428_003_E01      L6 IT
    ## F2S4_160428_003_F01         NP
    ## F2S4_160428_003_G01         NP
    ## F2S4_160428_003_H01         NP
    ## F2S4_160428_004_A01         L4
    ## F2S4_160428_004_B01    L2/3 IT
    ## F2S4_160428_004_C01         NP
    ## F2S4_160428_004_D01         NP
    ## F2S4_160428_004_E01         NP
    ## F2S4_160428_004_F01         NP
    ## F2S4_160428_004_G01         NP
    ## F2S4_160428_004_H01         NP
    ## F2S4_160428_005_A01    L2/3 IT
    ## F2S4_160428_005_B01    L2/3 IT
    ## F2S4_160428_005_C01    L2/3 IT
    ## F2S4_160428_005_D01    L2/3 IT
    ## F2S4_160428_005_E01    L2/3 IT
    ## F2S4_160428_005_F01    L2/3 IT
    ## F2S4_160428_005_G01    L2/3 IT
    ## F2S4_160428_005_H01    L2/3 IT
    ## F2S4_160428_006_A01    L2/3 IT
    ## F2S4_160428_006_B01    L2/3 IT
    ## F2S4_160428_006_C01    L2/3 IT
    ## F2S4_160428_006_D01    L2/3 IT
    ## F2S4_160428_006_E01    L2/3 IT
    ## F2S4_160428_006_F01    L2/3 IT
    ## F2S4_160428_006_G01    L2/3 IT
    ## F2S4_160428_006_H01         L4
    ## F2S4_160428_007_A01    L2/3 IT
    ## F2S4_160428_007_B01    L2/3 IT
    ## F2S4_160428_007_C01    L2/3 IT
    ## F2S4_160428_007_D01    L2/3 IT
    ## F2S4_160428_007_E01    L2/3 IT
    ## F2S4_160428_007_F01    L2/3 IT
    ## F2S4_160428_007_G01    L2/3 IT
    ## F2S4_160428_007_H01    L2/3 IT
    ## F2S4_160428_008_A01      L6 CT
    ## F2S4_160428_008_B01      L6 CT
    ## F2S4_160428_008_C01      L6 IT
    ## F2S4_160428_008_D01      L6 IT
    ## F2S4_160428_008_E01      L6 CT
    ## F2S4_160428_008_F01      L6 IT
    ## F2S4_160428_008_G01    L2/3 IT
    ## F2S4_160428_008_H01    L2/3 IT
    ## F2S4_160428_009_A01      L6 CT
    ## F2S4_160428_009_B01      L6 CT
    ## F2S4_160428_009_C01        L6b
    ## F2S4_160428_009_D01        L6b
    ## F2S4_160428_009_E01      L6 IT
    ## F2S4_160428_009_F01      L6 IT
    ## F2S4_160428_009_G01      L6 IT
    ## F2S4_160428_009_H01      L6 CT
    ## F2S4_160428_010_A01      L6 IT
    ## F2S4_160428_010_B01      L6 IT
    ## F2S4_160428_010_C01      L6 CT
    ## F2S4_160428_010_D01      L6 CT
    ## F2S4_160428_010_E01      L6 CT
    ## F2S4_160428_010_F01      L6 CT
    ## F2S4_160428_010_G01      L6 IT
    ## F2S4_160428_011_A01      L6 CT
    ## F2S4_160428_011_B01      L6 IT
    ## F2S4_160428_011_C01      L6 CT
    ## F2S4_160428_011_D01      L6 CT
    ## F2S4_160428_011_H01      L6 IT
    ## F2S4_160428_012_A01      L6 IT
    ## F2S4_160428_012_B01      L6 IT
    ## F2S4_160428_012_C01      L6 CT
    ## F2S4_160428_012_D01      L6 CT
    ## F2S4_160428_012_E01      L6 IT
    ## F2S4_160428_012_F01      L6 IT
    ## F2S4_160428_012_G01      L6 IT
    ## F2S4_160428_012_H01      L6 IT
    ## F2S4_160428_013_A01      L6 IT
    ## F2S4_160428_013_B01      L6 CT
    ## F2S4_160428_013_C01      L6 IT
    ## F2S4_160428_013_D01      L6 IT
    ## F2S4_160428_013_E01      L6 IT
    ## F2S4_160428_013_F01      L6 IT
    ## F2S4_160428_013_G01      L6 IT
    ## F2S4_160428_013_H01      L6 IT
    ## F2S4_160428_014_A01      L6 IT
    ## F2S4_160428_014_B01      L6 CT
    ## F2S4_160428_014_C01      L6 IT
    ## F2S4_160428_014_D01      L6 CT
    ## F2S4_160428_014_E01      L6 IT
    ## F2S4_160428_014_F01      L6 IT
    ## F2S4_160428_014_H01      L6 IT
    ## F2S4_160428_015_A01         L4
    ## F2S4_160428_015_B01         L4
    ## F2S4_160428_015_C01         L4
    ## F2S4_160428_015_D01         L4
    ## F2S4_160428_015_E01         L4
    ## F2S4_160428_015_F01         L4
    ## F2S4_160428_015_G01         L4
    ## F2S4_160428_015_H01         L4
    ## F2S4_160428_016_A01         L4
    ## F2S4_160428_016_B01      L5 IT
    ## F2S4_160428_016_C01         L4
    ## F2S4_160428_016_D01         L4
    ## F2S4_160428_016_E01         L4
    ## F2S4_160428_016_F01         L4
    ## F2S4_160428_016_G01         L4
    ## F2S4_160428_016_H01         L4
    ## F2S4_160428_017_A01         L4
    ## F2S4_160428_017_B01         L4
    ## F2S4_160428_017_C01         L4
    ## F2S4_160428_017_D01         L4
    ## F2S4_160428_017_E01         L4
    ## F2S4_160428_017_F01         L4
    ## F2S4_160428_017_G01         L4
    ## F2S4_160428_017_H01         L4
    ## F2S4_160428_018_A01         L4
    ## F2S4_160428_018_B01         L4
    ## F2S4_160428_018_C01         L4
    ## F2S4_160428_018_D01         L4
    ## F2S4_160428_018_E01         L4
    ## F2S4_160428_018_F01         L4
    ## F2S4_160428_018_G01         L4
    ## F2S4_160428_018_H01         L4
    ## F2S4_160428_019_A01         L4
    ## F2S4_160428_019_B01         L4
    ## F2S4_160428_019_C01         L4
    ## F2S4_160428_019_D01         L4
    ## F2S4_160428_019_E01         L4
    ## F2S4_160428_019_F01         L4
    ## F2S4_160428_019_G01         L4
    ## F2S4_160428_019_H01         L4
    ## F2S4_160428_020_A01      L6 IT
    ## F2S4_160428_020_B01      L6 CT
    ## F2S4_160428_020_C01      L6 IT
    ## F2S4_160428_020_D01      L6 IT
    ## F2S4_160428_020_E01      L6 IT
    ## F2S4_160428_020_F01      L6 IT
    ## F2S4_160428_020_G01      L6 CT
    ## F2S4_160428_020_H01         L4
    ## F2S4_160428_021_A01      L6 IT
    ## F2S4_160428_021_B01      L6 IT
    ## F2S4_160428_021_C01      L6 IT
    ## F2S4_160428_021_D01      L6 IT
    ## F2S4_160428_021_E01      L6 CT
    ## F2S4_160428_021_F01      L6 IT
    ## F2S4_160428_021_G01      L6 IT
    ## F2S4_160428_021_H01      L6 CT
    ## F2S4_160428_022_A01      L6 IT
    ## F2S4_160428_022_C01      L6 IT
    ## F2S4_160428_022_D01      L6 IT
    ## F2S4_160428_022_E01      L6 CT
    ## F2S4_160428_022_F01      L6 IT
    ## F2S4_160428_022_G01      L6 IT
    ## F2S4_160428_022_H01      L6 IT
    ## F2S4_160428_023_A01      L6 IT
    ## F2S4_160428_023_B01      L6 IT
    ## F2S4_160428_023_C01      L6 IT
    ## F2S4_160428_023_D01      L6 IT
    ## F2S4_160428_023_E01      L6 IT
    ## F2S4_160428_023_F01      L6 IT
    ## F2S4_160428_023_G01      L6 IT
    ## F2S4_160428_023_H01      L6 IT
    ## F2S4_160428_024_A01      L6 IT
    ## F2S4_160428_024_B01      L6 IT
    ## F2S4_160428_024_C01      L6 IT
    ## F2S4_160428_024_D01      L6 IT
    ## F2S4_160428_024_E01      L6 IT
    ## F2S4_160428_024_F01      L6 CT
    ## F2S4_160428_024_G01      L6 IT
    ## F2S4_160428_024_H01      L6 IT
    ## F2S4_160429_001_A01      Lamp5
    ## F2S4_160429_001_B01      Lamp5
    ## F2S4_160429_001_D01      Lamp5
    ## F2S4_160429_001_E01      Lamp5
    ## F2S4_160429_001_F01      Lamp5
    ## F2S4_160429_001_G01      Lamp5
    ## F2S4_160429_001_H01      Lamp5
    ## F2S4_160429_002_A01      Lamp5
    ## F2S4_160429_002_B01      Lamp5
    ## F2S4_160429_002_C01      Lamp5
    ## F2S4_160429_002_D01      Lamp5
    ## F2S4_160429_002_E01      Lamp5
    ## F2S4_160429_002_F01      Lamp5
    ## F2S4_160429_002_G01      Lamp5
    ## F2S4_160429_002_H01      Lamp5
    ## F2S4_160429_003_A01      Lamp5
    ## F2S4_160429_003_B01      Lamp5
    ## F2S4_160429_003_C01      Lamp5
    ## F2S4_160429_003_D01      Lamp5
    ## F2S4_160429_003_E01        Vip
    ## F2S4_160429_003_F01        Vip
    ## F2S4_160429_003_G01        Vip
    ## F2S4_160429_003_H01        Sst
    ## F2S4_160429_004_A01        Vip
    ## F2S4_160429_004_B01        Sst
    ## F2S4_160429_004_C01      Lamp5
    ## F2S4_160429_004_D01        Sst
    ## F2S4_160429_004_E01        Sst
    ## F2S4_160429_004_F01        Sst
    ## F2S4_160429_004_G01        Sst
    ## F2S4_160429_004_H01      Pvalb
    ## F2S4_160429_005_A01        Sst
    ## F2S4_160429_005_B01        Sst
    ## F2S4_160429_005_C01      Lamp5
    ## F2S4_160429_005_D01      Pvalb
    ## F2S4_160429_005_E01      Lamp5
    ## F2S4_160429_005_F01        Sst
    ## F2S4_160429_005_G01        Sst
    ## F2S4_160429_005_H01        Sst
    ## F2S4_160429_006_A01        Sst
    ## F2S4_160429_006_B01        Sst
    ## F2S4_160429_006_C01        Vip
    ## F2S4_160429_006_D01      Pvalb
    ## F2S4_160429_006_E01      Pvalb
    ## F2S4_160429_006_F01      Lamp5
    ## F2S4_160429_006_G01      Lamp5
    ## F2S4_160429_006_H01        Sst
    ## F2S4_160429_007_A01      Pvalb
    ## F2S4_160429_007_B01        Vip
    ## F2S4_160429_007_C01      Lamp5
    ## F2S4_160429_007_D01        Vip
    ## F2S4_160429_007_E01        Vip
    ## F2S4_160429_007_F01      Lamp5
    ## F2S4_160429_007_G01      Lamp5
    ## F2S4_160429_007_H01        Vip
    ## F2S4_160429_008_A01      Lamp5
    ## F2S4_160429_008_B01       Sncg
    ## F2S4_160429_008_C01        Vip
    ## F2S4_160429_008_D01        Sst
    ## F2S4_160429_008_E01      Lamp5
    ## F2S4_160429_008_F01        Vip
    ## F2S4_160429_008_G01      Lamp5
    ## F2S4_160429_008_H01      Lamp5
    ## F2S4_160429_009_A01      Lamp5
    ## F2S4_160429_009_B01        Vip
    ## F2S4_160429_009_C01      Lamp5
    ## F2S4_160429_009_D01        Vip
    ## F2S4_160429_009_E01        Vip
    ## F2S4_160429_009_F01      Lamp5
    ## F2S4_160429_009_G01      Lamp5
    ## F2S4_160429_009_H01        Vip
    ## F2S4_160429_010_A01        Vip
    ## F2S4_160429_010_B01        Vip
    ## F2S4_160429_010_C01        Vip
    ## F2S4_160429_010_D01      Pvalb
    ## F2S4_160429_010_E01        Sst
    ## F2S4_160429_010_F01        Sst
    ## F2S4_160429_010_G01      Lamp5
    ## F2S4_160429_010_H01        Sst
    ## F2S4_160429_011_A01       Endo
    ## F2S4_160429_011_B01        Vip
    ## F2S4_160429_011_C01        Sst
    ## F2S4_160429_011_D01      Lamp5
    ## F2S4_160429_011_E01      Pvalb
    ## F2S4_160429_011_F01        Sst
    ## F2S4_160429_011_G01        Vip
    ## F2S4_160429_011_H01        Sst
    ## F2S4_160429_012_A01        Sst
    ## F2S4_160429_012_B01        Vip
    ## F2S4_160429_012_C01        Sst
    ## F2S4_160429_012_D01      Pvalb
    ## F2S4_160429_012_E01        Sst
    ## F2S4_160429_012_F01      Pvalb
    ## F2S4_160429_012_G01        Sst
    ## F2S4_160429_012_H01        Sst
    ## F2S4_160429_013_A01      Lamp5
    ## F2S4_160429_013_B01        Vip
    ## F2S4_160429_013_C01        Sst
    ## F2S4_160429_013_D01      Pvalb
    ## F2S4_160429_013_E01        Sst
    ## F2S4_160429_013_F01       Endo
    ## F2S4_160429_013_G01        Sst
    ## F2S4_160429_013_H01      Pvalb
    ## F2S4_160429_014_A01      Pvalb
    ## F2S4_160429_014_B01        Sst
    ## F2S4_160429_014_C01        Sst
    ## F2S4_160429_014_D01        Sst
    ## F2S4_160429_014_E01        Sst
    ## F2S4_160429_014_F01      Pvalb
    ## F2S4_160429_014_G01        Vip
    ## F2S4_160429_014_H01        Sst
    ## F2S4_160429_015_A01        Sst
    ## F2S4_160429_015_B01      Lamp5
    ## F2S4_160429_015_C01        Sst
    ## F2S4_160429_015_D01      Pvalb
    ## F2S4_160429_015_E01        Vip
    ## F2S4_160429_015_F01        Vip
    ## F2S4_160429_015_G01      Pvalb
    ## F2S4_160429_015_H01        Vip
    ## F2S4_160429_016_A01        Vip
    ## F2S4_160429_016_B01        Sst
    ## F2S4_160429_016_C01        Vip
    ## F2S4_160429_016_D01      Lamp5
    ## F2S4_160429_016_E01        Vip
    ## F2S4_160429_016_F01      Lamp5
    ## F2S4_160429_016_G01        Vip
    ## F2S4_160429_016_H01        Vip
    ## F2S4_160429_017_A01        Vip
    ## F2S4_160429_017_B01        Vip
    ## F2S4_160429_017_C01      Lamp5
    ## F2S4_160429_017_D01        Vip
    ## F2S4_160429_017_E01        Vip
    ## F2S4_160429_017_F01        Vip
    ## F2S4_160429_017_G01       Sncg
    ## F2S4_160429_017_H01      Lamp5
    ## F2S4_160429_018_A01        Sst
    ## F2S4_160429_018_B01        Vip
    ## F2S4_160429_018_C01        Vip
    ## F2S4_160429_018_D01      Meis2
    ## F2S4_160429_018_E01        Vip
    ## F2S4_160429_018_F01        Sst
    ## F2S4_160429_018_G01        Vip
    ## F2S4_160429_018_H01        Vip
    ## F2S4_160429_019_A01        Vip
    ## F2S4_160429_019_B01        Vip
    ## F2S4_160429_019_C01        Vip
    ## F2S4_160429_019_D01        Vip
    ## F2S4_160429_019_E01        Vip
    ## F2S4_160429_019_F01        Vip
    ## F2S4_160429_019_G01        Vip
    ## F2S4_160429_019_H01        Sst
    ## F2S4_160429_020_A01      Lamp5
    ## F2S4_160429_020_B01        Vip
    ## F2S4_160429_020_C01        Vip
    ## F2S4_160429_020_D01      Lamp5
    ## F2S4_160429_020_E01        Vip
    ## F2S4_160429_020_F01       Sncg
    ## F2S4_160429_020_G01      Lamp5
    ## F2S4_160429_020_H01      Pvalb
    ## F2S4_160429_021_A01        Vip
    ## F2S4_160429_021_B01      Pvalb
    ## F2S4_160429_021_C01      Lamp5
    ## F2S4_160429_021_D01       Sncg
    ## F2S4_160429_021_E01        Vip
    ## F2S4_160429_021_F01      Lamp5
    ## F2S4_160429_021_G01      Lamp5
    ## F2S4_160429_021_H01        Vip
    ## F2S4_160429_022_A01        Sst
    ## F2S4_160429_022_B01        Sst
    ## F2S4_160429_022_C01        Sst
    ## F2S4_160429_022_D01        Sst
    ## F2S4_160429_022_F01        Sst
    ## F2S4_160429_022_G01        Vip
    ## F2S4_160429_022_H01        Sst
    ## F2S4_160502_001_A01        Vip
    ## F2S4_160502_001_B01        Vip
    ## F2S4_160502_001_C01      Lamp5
    ## F2S4_160502_001_D01       Sncg
    ## F2S4_160502_001_E01      Lamp5
    ## F2S4_160502_001_F01        Vip
    ## F2S4_160502_001_G01      Pvalb
    ## F2S4_160502_001_H01        Vip
    ## F2S4_160502_002_A01      Lamp5
    ## F2S4_160502_002_B01      Pvalb
    ## F2S4_160502_002_C01        Vip
    ## F2S4_160502_002_D01        Vip
    ## F2S4_160502_002_E01        Vip
    ## F2S4_160502_002_F01        Vip
    ## F2S4_160502_002_G01        Sst
    ## F2S4_160502_002_H01        Vip
    ## F2S4_160502_003_A01        Vip
    ## F2S4_160502_003_B01        Vip
    ## F2S4_160502_003_C01        Vip
    ## F2S4_160502_003_D01        Vip
    ## F2S4_160502_003_E01        Vip
    ## F2S4_160502_003_F01        Vip
    ## F2S4_160502_003_G01       Sncg
    ## F2S4_160502_003_H01        Vip
    ## F2S4_160502_004_A01        Vip
    ## F2S4_160502_004_B01        Vip
    ## F2S4_160502_004_C01      Lamp5
    ## F2S4_160502_004_D01      Lamp5
    ## F2S4_160502_004_E01       Sncg
    ## F2S4_160502_004_F01        Vip
    ## F2S4_160502_004_G01        Vip
    ## F2S4_160502_004_H01        Vip
    ## F2S4_160502_005_A01        Vip
    ## F2S4_160502_005_B01        Vip
    ## F2S4_160502_005_C01      Pvalb
    ## F2S4_160502_005_D01        Vip
    ## F2S4_160502_005_E01        Vip
    ## F2S4_160502_005_F01        Vip
    ## F2S4_160502_005_G01        Sst
    ## F2S4_160502_005_H01        Vip
    ## F2S4_160502_006_A01      Lamp5
    ## F2S4_160502_006_B01      Lamp5
    ## F2S4_160502_006_C01        Vip
    ## F2S4_160502_006_E01        Vip
    ## F2S4_160502_006_F01        Sst
    ## F2S4_160502_006_G01        Vip
    ## F2S4_160502_006_H01      Lamp5
    ## F2S4_160502_007_A01      Lamp5
    ## F2S4_160502_007_B01      Lamp5
    ## F2S4_160502_007_C01      Lamp5
    ## F2S4_160502_007_D01      Lamp5
    ## F2S4_160502_007_E01      Lamp5
    ## F2S4_160502_007_F01      Lamp5
    ## F2S4_160502_007_G01      Lamp5
    ## F2S4_160502_007_H01      Lamp5
    ## F2S4_160502_008_A01        Vip
    ## F2S4_160502_008_B01       Sncg
    ## F2S4_160502_008_C01        Vip
    ## F2S4_160502_008_D01      Lamp5
    ## F2S4_160502_008_E01        Vip
    ## F2S4_160502_008_F01      Lamp5
    ## F2S4_160502_008_G01      Lamp5
    ## F2S4_160502_008_H01      Lamp5
    ## F2S4_160502_009_A01        Vip
    ## F2S4_160502_009_B01        Vip
    ## F2S4_160502_009_C01        Vip
    ## F2S4_160502_009_D01      Pvalb
    ## F2S4_160502_009_E01        Vip
    ## F2S4_160502_009_F01        Vip
    ## F2S4_160502_009_G01      Lamp5
    ## F2S4_160502_009_H01        Vip
    ## F2S4_160502_010_A01      Pvalb
    ## F2S4_160502_010_B01      Lamp5
    ## F2S4_160502_010_C01        Vip
    ## F2S4_160502_010_D01        Vip
    ## F2S4_160502_010_E01        Vip
    ## F2S4_160502_010_F01      Lamp5
    ## F2S4_160502_010_G01        Vip
    ## F2S4_160502_010_H01        Vip
    ## F2S4_160502_011_A01        Sst
    ## F2S4_160502_011_B01        Sst
    ## F2S4_160502_011_C01        Vip
    ## F2S4_160502_011_D01        Vip
    ## F2S4_160502_011_E01        Sst
    ## F2S4_160502_011_F01        Vip
    ## F2S4_160502_011_G01      Pvalb
    ## F2S4_160502_011_H01      Lamp5
    ## F2S4_160502_012_A01        Vip
    ## F2S4_160502_012_B01        Vip
    ## F2S4_160502_012_C01        Vip
    ## F2S4_160502_012_D01        Vip
    ## F2S4_160502_012_E01       Sncg
    ## F2S4_160502_012_F01        Vip
    ## F2S4_160502_012_G01      Lamp5
    ## F2S4_160502_012_H01        Vip
    ## F2S4_160502_013_A01        Vip
    ## F2S4_160502_013_B01      Lamp5
    ## F2S4_160502_013_C01        Vip
    ## F2S4_160502_013_D01        Vip
    ## F2S4_160502_013_E01      Lamp5
    ## F2S4_160502_013_F01      Lamp5
    ## F2S4_160502_013_G01        Vip
    ## F2S4_160502_013_H01        Vip
    ## F2S4_160502_014_A01        Vip
    ## F2S4_160502_014_B01      Lamp5
    ## F2S4_160502_014_C01        Vip
    ## F2S4_160502_014_D01      Lamp5
    ## F2S4_160502_014_E01        Vip
    ## F2S4_160502_014_F01        Vip
    ## F2S4_160502_014_G01        Vip
    ## F2S4_160502_014_H01        Vip
    ## F2S4_160502_015_A01        Vip
    ## F2S4_160502_015_B01      Pvalb
    ## F2S4_160502_015_C01        Sst
    ## F2S4_160502_015_D01        Sst
    ## F2S4_160502_015_E01      Lamp5
    ## F2S4_160502_015_F01        Sst
    ## F2S4_160502_015_G01        Sst
    ## F2S4_160502_015_H01        Sst
    ## F2S4_160502_016_A01        Sst
    ## F2S4_160502_016_B01        Vip
    ## F2S4_160502_016_C01        Sst
    ## F2S4_160502_016_D01      Pvalb
    ## F2S4_160502_016_E01      Meis2
    ## F2S4_160502_016_F01      Pvalb
    ## F2S4_160502_016_G01        Vip
    ## F2S4_160502_016_H01        Sst
    ## F2S4_160502_017_A01        Vip
    ## F2S4_160502_017_B01      Pvalb
    ## F2S4_160502_017_C01        Sst
    ## F2S4_160502_017_D01        Vip
    ## F2S4_160502_017_E01      Lamp5
    ## F2S4_160502_017_F01      Lamp5
    ## F2S4_160502_017_G01      Lamp5
    ## F2S4_160502_017_H01      Lamp5
    ## F2S4_160502_018_A01      Lamp5
    ## F2S4_160502_018_B01      Lamp5
    ## F2S4_160502_018_C01      Lamp5
    ## F2S4_160502_018_D01      Meis2
    ## F2S4_160502_018_E01      Lamp5
    ## F2S4_160502_018_F01      Lamp5
    ## F2S4_160502_018_G01      Lamp5
    ## F2S4_160502_018_H01      Lamp5
    ## F2S4_160502_019_A01      Lamp5
    ## F2S4_160502_019_B01      Lamp5
    ## F2S4_160502_019_C01        Vip
    ## F2S4_160502_019_D01      Lamp5
    ## F2S4_160502_019_E01      Lamp5
    ## F2S4_160502_019_F01      Lamp5
    ## F2S4_160502_019_G01      Lamp5
    ## F2S4_160502_019_H01      Lamp5
    ## F2S4_160502_020_A01      Lamp5
    ## F2S4_160502_020_B01      Lamp5
    ## F2S4_160502_020_C01      Lamp5
    ## F2S4_160502_020_D01        Vip
    ## F2S4_160502_020_F01        Vip
    ## F2S4_160502_020_G01      Pvalb
    ## F2S4_160502_020_H01        Sst
    ## F2S4_160504_001_A01        Sst
    ## F2S4_160504_001_B01        Sst
    ## F2S4_160504_001_C01        Sst
    ## F2S4_160504_001_D01        Sst
    ## F2S4_160504_001_E01        Sst
    ## F2S4_160504_001_F01      Lamp5
    ## F2S4_160504_001_G01      Lamp5
    ## F2S4_160504_001_H01      Lamp5
    ## F2S4_160504_002_A01      Lamp5
    ## F2S4_160504_002_B01      Lamp5
    ## F2S4_160504_002_C01        Vip
    ## F2S4_160504_002_D01      Lamp5
    ## F2S4_160504_002_E01      Lamp5
    ## F2S4_160504_002_F01       Endo
    ## F2S4_160504_002_G01       Endo
    ## F2S4_160504_002_H01      Lamp5
    ## F2S4_160504_003_A01      Lamp5
    ## F2S4_160504_003_B01      Lamp5
    ## F2S4_160504_003_C01        Vip
    ## F2S4_160504_003_D01      Lamp5
    ## F2S4_160504_003_E01      Lamp5
    ## F2S4_160504_003_F01      Lamp5
    ## F2S4_160504_003_G01      Lamp5
    ## F2S4_160504_003_H01      Lamp5
    ## F2S4_160504_004_A01        Vip
    ## F2S4_160504_004_B01      Lamp5
    ## F2S4_160504_004_C01      Lamp5
    ## F2S4_160504_004_D01        Vip
    ## F2S4_160504_004_E01      Lamp5
    ## F2S4_160504_004_F01      Lamp5
    ## F2S4_160504_004_G01      Lamp5
    ## F2S4_160504_004_H01        Vip
    ## F2S4_160504_005_A01       Endo
    ## F2S4_160504_005_B01      Lamp5
    ## F2S4_160504_005_C01      Lamp5
    ## F2S4_160504_005_D01      Lamp5
    ## F2S4_160504_005_E01      Lamp5
    ## F2S4_160504_005_F01       Endo
    ## F2S4_160504_005_G01      Lamp5
    ## F2S4_160504_005_H01      Lamp5
    ## F2S4_160504_006_A01      Lamp5
    ## F2S4_160504_006_B01       Endo
    ## F2S4_160504_006_C01      Lamp5
    ## F2S4_160504_006_D01       Endo
    ## F2S4_160504_006_E01        Vip
    ## F2S4_160504_006_F01       Endo
    ## F2S4_160504_006_G01      Lamp5
    ## F2S4_160504_007_A01        Sst
    ## F2S4_160504_007_B01      Pvalb
    ## F2S4_160504_007_C01       Endo
    ## F2S4_160504_007_D01       Endo
    ## F2S4_160504_007_E01       Endo
    ## F2S4_160504_007_F01      Lamp5
    ## F2S4_160504_007_G01      Pvalb
    ## F2S4_160504_007_H01        Sst
    ## F2S4_160504_008_A01       Endo
    ## F2S4_160504_008_B01        Sst
    ## F2S4_160504_008_C01        Vip
    ## F2S4_160504_008_D01        Sst
    ## F2S4_160504_008_E01        Sst
    ## F2S4_160504_008_F01      Lamp5
    ## F2S4_160504_008_H01      Lamp5
    ## F2S4_160504_009_A01      Lamp5
    ## F2S4_160504_009_B01      Lamp5
    ## F2S4_160504_009_C01      Lamp5
    ## F2S4_160504_009_D01      Lamp5
    ## F2S4_160504_009_E01      Lamp5
    ## F2S4_160504_009_F01      Lamp5
    ## F2S4_160504_009_G01      Lamp5
    ## F2S4_160504_009_H01      Lamp5
    ## F2S4_160504_010_A01       Endo
    ## F2S4_160504_010_B01      Lamp5
    ## F2S4_160504_010_C01      Lamp5
    ## F2S4_160504_010_D01      Lamp5
    ## F2S4_160504_010_E01       Sncg
    ## F2S4_160504_010_F01      Lamp5
    ## F2S4_160504_010_G01      Lamp5
    ## F2S4_160504_010_H01        Vip
    ## F2S4_160504_011_A01      Pvalb
    ## F2S4_160504_011_B01        Sst
    ## F2S4_160504_011_C01      Pvalb
    ## F2S4_160504_011_D01        Sst
    ## F2S4_160504_011_E01        Sst
    ## F2S4_160504_011_F01        Vip
    ## F2S4_160504_011_G01      Meis2
    ## F2S4_160512_001_A01      Pvalb
    ## F2S4_160512_001_B01      Pvalb
    ## F2S4_160512_001_C01      Pvalb
    ## F2S4_160512_001_D01      Pvalb
    ## F2S4_160512_001_E01      Pvalb
    ## F2S4_160512_001_F01      Pvalb
    ## F2S4_160512_001_G01      Pvalb
    ## F2S4_160512_001_H01      Pvalb
    ## F2S4_160512_002_A01      Pvalb
    ## F2S4_160512_002_C01      Pvalb
    ## F2S4_160512_002_D01      Pvalb
    ## F2S4_160512_002_E01      Pvalb
    ## F2S4_160512_002_F01      Pvalb
    ## F2S4_160512_002_G01      Pvalb
    ## F2S4_160512_002_H01      Pvalb
    ## F2S4_160512_003_A01      Pvalb
    ## F2S4_160512_003_B01      Pvalb
    ## F2S4_160512_003_C01      Pvalb
    ## F2S4_160512_003_D01      Pvalb
    ## F2S4_160512_003_E01      Pvalb
    ## F2S4_160512_003_F01      Pvalb
    ## F2S4_160516_001_A01         NP
    ## F2S4_160516_001_B01         NP
    ## F2S4_160516_001_C01         NP
    ## F2S4_160516_001_D01         NP
    ## F2S4_160516_001_E01      L6 IT
    ## F2S4_160516_001_F01        Sst
    ## F2S4_160516_001_G01         NP
    ## F2S4_160516_001_H01        Sst
    ## F2S4_160516_002_A01         NP
    ## F2S4_160516_002_B01         NP
    ## F2S4_160516_002_C01      L6 CT
    ## F2S4_160516_002_D01        Sst
    ## F2S4_160516_002_E01         NP
    ## F2S4_160516_002_F01      L6 IT
    ## F2S4_160516_002_G01        Sst
    ## F2S4_160516_002_H01      L5 IT
    ## F2S4_160516_003_A01        Sst
    ## F2S4_160516_003_B01        Vip
    ## F2S4_160516_003_C01      L5 PT
    ## F2S4_160516_003_D01      L5 PT
    ## F2S4_160516_003_E01         NP
    ## F2S4_160516_003_F01        Sst
    ## F2S4_160516_003_G01      L5 IT
    ## F2S4_160516_003_H01         NP
    ## F2S4_160516_004_A01      L6 IT
    ## F2S4_160516_004_B01         NP
    ## F2S4_160516_004_C01      L5 IT
    ## F2S4_160516_004_D01         NP
    ## F2S4_160516_004_E01         NP
    ## F2S4_160516_004_F01         NP
    ## F2S4_160516_004_G01      L6 IT
    ## F2S4_160516_004_H01      L5 IT
    ## F2S4_160516_005_A01         NP
    ## F2S4_160516_005_B01        Sst
    ## F2S4_160516_005_C01        Sst
    ## F2S4_160516_005_D01      L5 IT
    ## F2S4_160516_005_E01      L5 PT
    ## F2S4_160516_005_F01         NP
    ## F2S4_160516_005_G01        Vip
    ## F2S4_160516_005_H01        Sst
    ## F2S4_160516_006_A01      L5 IT
    ## F2S4_160516_006_B01         NP
    ## F2S4_160516_006_C01         NP
    ## F2S4_160516_006_D01        Sst
    ## F2S4_160516_006_E01      L6 IT
    ## F2S4_160516_006_F01      L5 IT
    ## F2S4_160516_006_G01        Sst
    ## F2S4_160516_006_H01      L6 IT
    ## F2S4_160516_007_A01      L6 IT
    ## F2S4_160516_007_B01         NP
    ## F2S4_160516_007_C01      L5 PT
    ## F2S4_160516_007_D01        Sst
    ## F2S4_160516_007_E01        Sst
    ## F2S4_160516_007_F01         NP
    ## F2S4_160516_007_G01         NP
    ## F2S4_160516_007_H01      L6 IT
    ## F2S4_160516_008_A01    L2/3 IT
    ## F2S4_160516_008_B01    L2/3 IT
    ## F2S4_160516_008_C01    L2/3 IT
    ## F2S4_160516_008_D01      Lamp5
    ## F2S4_160516_008_E01    L2/3 IT
    ## F2S4_160516_008_F01        Vip
    ## F2S4_160516_008_G01        Vip
    ## F2S4_160516_008_H01      Lamp5
    ## F2S4_160516_009_A01    L2/3 IT
    ## F2S4_160516_009_B01        Vip
    ## F2S4_160516_009_C01        Vip
    ## F2S4_160516_009_D01    L2/3 IT
    ## F2S4_160516_009_E01    L2/3 IT
    ## F2S4_160516_009_F01    L2/3 IT
    ## F2S4_160516_009_G01    L2/3 IT
    ## F2S4_160516_009_H01    L2/3 IT
    ## F2S4_160516_010_A01    L2/3 IT
    ## F2S4_160516_010_B01        Vip
    ## F2S4_160516_010_C01    L2/3 IT
    ## F2S4_160516_010_D01      Lamp5
    ## F2S4_160516_010_E01    L2/3 IT
    ## F2S4_160516_010_F01        Sst
    ## F2S4_160516_010_G01        Vip
    ## F2S4_160516_010_H01    L2/3 IT
    ## F2S4_160516_011_A01    L2/3 IT
    ## F2S4_160516_011_B01    L2/3 IT
    ## F2S4_160516_011_C01    L2/3 IT
    ## F2S4_160516_011_D01         L4
    ## F2S4_160516_011_E01         L4
    ## F2S4_160516_011_F01      Lamp5
    ## F2S4_160516_011_G01         L4
    ## F2S4_160516_011_H01    L2/3 IT
    ## F2S4_160516_012_A01    L2/3 IT
    ## F2S4_160516_012_B01    L2/3 IT
    ## F2S4_160516_012_C01    L2/3 IT
    ## F2S4_160516_012_D01    L2/3 IT
    ## F2S4_160516_012_E01    L2/3 IT
    ## F2S4_160516_012_F01    L2/3 IT
    ## F2S4_160516_012_G01        Vip
    ## F2S4_160516_012_H01      Lamp5
    ## F2S4_160516_013_A01    L2/3 IT
    ## F2S4_160516_013_B01    L2/3 IT
    ## F2S4_160516_013_C01        Vip
    ## F2S4_160516_013_D01      Lamp5
    ## F2S4_160516_013_E01    L2/3 IT
    ## F2S4_160516_013_F01        Vip
    ## F2S4_160516_013_G01    L2/3 IT
    ## F2S4_160516_013_H01    L2/3 IT
    ## F2S4_160516_014_A01    L2/3 IT
    ## F2S4_160516_014_B01    L2/3 IT
    ## F2S4_160516_014_C01    L2/3 IT
    ## F2S4_160516_014_D01        Vip
    ## F2S4_160516_014_E01      Lamp5
    ## F2S4_160516_014_F01    L2/3 IT
    ## F2S4_160516_014_G01      Lamp5
    ## F2S4_160516_014_H01      Lamp5
    ## F2S4_160516_015_A01        Sst
    ## F2S4_160516_015_B01    L2/3 IT
    ## F2S4_160516_015_C01      Lamp5
    ## F2S4_160516_015_D01    L2/3 IT
    ## F2S4_160516_015_E01    L2/3 IT
    ## F2S4_160516_015_F01        Vip
    ## F2S4_160516_015_G01        Vip
    ## F2S4_160516_015_H01        Vip
    ## F2S4_160516_016_A01    L2/3 IT
    ## F2S4_160516_016_B01    L2/3 IT
    ## F2S4_160516_016_C01    L2/3 IT
    ## F2S4_160516_016_D01    L2/3 IT
    ## F2S4_160516_016_E01    L2/3 IT
    ## F2S4_160516_016_F01    L2/3 IT
    ## F2S4_160516_016_G01        Vip
    ## F2S4_160516_016_H01    L2/3 IT
    ## F2S4_160516_017_A01        Vip
    ## F2S4_160516_017_B01    L2/3 IT
    ## F2S4_160516_017_C01    L2/3 IT
    ## F2S4_160516_017_D01        Vip
    ## F2S4_160516_017_E01    L2/3 IT
    ## F2S4_160516_017_F01    L2/3 IT
    ## F2S4_160516_017_G01        Vip
    ## F2S4_160516_017_H01    L2/3 IT
    ## F2S4_160516_018_A01        Vip
    ## F2S4_160516_018_B01    L2/3 IT
    ## F2S4_160516_018_C01        Sst
    ## F2S4_160516_018_D01        Vip
    ## F2S4_160516_018_E01    L2/3 IT
    ## F2S4_160516_018_F01         L4
    ## F2S4_160516_018_G01    L2/3 IT
    ## F2S4_160516_018_H01    L2/3 IT
    ## F2S4_160516_019_A01        Vip
    ## F2S4_160516_019_B01         L4
    ## F2S4_160516_019_C01    L2/3 IT
    ## F2S4_160516_019_D01        Vip
    ## F2S4_160516_019_E01    L2/3 IT
    ## F2S4_160516_019_F01        Vip
    ## F2S4_160516_019_G01         L4
    ## F2S4_160516_019_H01      Lamp5
    ## F2S4_160516_020_A01    L2/3 IT
    ## F2S4_160516_020_B01    L2/3 IT
    ## F2S4_160516_020_C01        Vip
    ## F2S4_160516_020_D01    L2/3 IT
    ## F2S4_160516_020_E01    L2/3 IT
    ## F2S4_160516_020_F01    L2/3 IT
    ## F2S4_160516_020_G01         L4
    ## F2S4_160516_020_H01        Sst
    ## F2S4_160516_021_A01         L4
    ## F2S4_160516_021_B01    L2/3 IT
    ## F2S4_160516_021_C01    L2/3 IT
    ## F2S4_160516_021_D01    L2/3 IT
    ## F2S4_160516_021_E01         L4
    ## F2S4_160516_021_F01    L2/3 IT
    ## F2S4_160516_021_G01      Lamp5
    ## F2S4_160516_021_H01    L2/3 IT
    ## F2S4_160516_022_A01         L4
    ## F2S4_160516_022_B01         L4
    ## F2S4_160516_022_C01         L4
    ## F2S4_160516_022_D01         L4
    ## F2S4_160516_022_E01      L5 IT
    ## F2S4_160516_022_G01         L4
    ## F2S4_160516_023_A01         L4
    ## F2S4_160516_023_B01         L4
    ## F2S4_160516_023_C01      L5 IT
    ## F2S4_160516_023_D01         L4
    ## F2S4_160516_023_E01         L4
    ## F2S4_160516_023_G01      L5 IT
    ## F2S4_160516_023_H01         L4
    ## F2S4_160516_024_A01      L5 IT
    ## F2S4_160516_024_B01         L4
    ## F2S4_160516_024_C01         L4
    ## F2S4_160516_024_D01         L4
    ## F2S4_160516_024_F01         L4
    ## F2S4_160516_024_G01         L4
    ## F2S4_160516_024_H01         L4
    ## F2S4_160516_025_B01         L4
    ## F2S4_160516_025_C01         L4
    ## F2S4_160516_025_D01         L4
    ## F2S4_160516_025_E01         L4
    ## F2S4_160516_025_F01         L4
    ## F2S4_160516_025_G01         L4
    ## F2S4_160516_025_H01         L4
    ## F2S4_160516_026_A01         L4
    ## F2S4_160516_026_B01         L4
    ## F2S4_160516_026_D01         L4
    ## F2S4_160516_026_E01         L4
    ## F2S4_160516_026_F01         L4
    ## F2S4_160516_026_G01         L4
    ## F2S4_160516_026_H01         L4
    ## F2S4_160516_027_A01         L4
    ## F2S4_160516_027_B01         L4
    ## F2S4_160516_027_C01         L4
    ## F2S4_160516_027_D01         L4
    ## F2S4_160516_027_G01         L4
    ## F2S4_160516_027_H01         L4
    ## F2S4_160516_028_A01         L4
    ## F2S4_160516_028_B01         L4
    ## F2S4_160516_028_C01         L4
    ## F2S4_160516_028_D01         L4
    ## F2S4_160516_028_E01         L4
    ## F2S4_160516_028_F01         L4
    ## F2S4_160516_028_H01         L4
    ## F2S4_160516_029_A01      L6 IT
    ## F2S4_160516_029_B01      L6 CT
    ## F2S4_160516_029_C01      L6 IT
    ## F2S4_160516_029_D01      L6 IT
    ## F2S4_160516_029_E01        Vip
    ## F2S4_160516_029_F01      L6 IT
    ## F2S4_160516_029_G01      L6 IT
    ## F2S4_160516_029_H01      L6 IT
    ## F2S4_160516_030_A01      L6 IT
    ## F2S4_160516_030_B01      L6 IT
    ## F2S4_160516_030_C01      L6 IT
    ## F2S4_160516_030_D01      L6 CT
    ## F2S4_160516_030_E01      L6 IT
    ## F2S4_160516_030_F01      L6 CT
    ## F2S4_160516_030_G01      L6 IT
    ## F2S4_160516_030_H01      L6 IT
    ## F2S4_160516_031_A01      L6 IT
    ## F2S4_160516_031_B01      L6 IT
    ## F2S4_160516_031_C01      L6 IT
    ## F2S4_160516_031_D01      L6 CT
    ## F2S4_160516_031_E01      L6 IT
    ## F2S4_160516_031_F01      L6 IT
    ## F2S4_160516_031_G01      L6 IT
    ## F2S4_160516_032_A01      L6 IT
    ## F2S4_160516_032_B01      L6 IT
    ## F2S4_160516_032_C01      L6 IT
    ## F2S4_160516_032_D01      L6 IT
    ## F2S4_160516_032_E01      L6 CT
    ## F2S4_160516_032_F01      L6 IT
    ## F2S4_160516_032_G01      L6 CT
    ## F2S4_160516_032_H01      L5 IT
    ## F2S4_160516_033_A01      L6 CT
    ## F2S4_160516_033_B01      L6 IT
    ## F2S4_160516_033_C01      L6 IT
    ## F2S4_160516_033_D01      L6 CT
    ## F2S4_160516_033_E01      L6 CT
    ## F2S4_160516_033_F01      L6 IT
    ## F2S4_160516_033_G01      L6 CT
    ## F2S4_160516_033_H01      L6 IT
    ## F2S4_160516_034_A01      L6 CT
    ## F2S4_160516_034_B01      L6 CT
    ## F2S4_160516_034_C01      L6 CT
    ## F2S4_160516_034_D01      L6 IT
    ## F2S4_160516_034_E01      L6 IT
    ## F2S4_160516_034_F01      L6 IT
    ## F2S4_160516_034_G01      L6 IT
    ## F2S4_160516_034_H01      L6 CT
    ## F2S4_160516_035_A01      L6 IT
    ## F2S4_160516_035_B01      L6 CT
    ## F2S4_160516_035_C01      L6 IT
    ## F2S4_160516_035_D01      L6 IT
    ## F2S4_160516_035_E01      L6 CT
    ## F2S4_160516_035_F01      L6 IT
    ## F2S4_160516_035_G01      L6 IT
    ## F2S4_160516_035_H01        Vip
    ## F2S4_160516_036_A01      L6 CT
    ## F2S4_160516_036_B01      L6 CT
    ## F2S4_160516_036_C01      L6 IT
    ## F2S4_160516_036_D01      L6 IT
    ## F2S4_160516_036_E01      L6 IT
    ## F2S4_160516_036_F01      L6 IT
    ## F2S4_160516_036_G01      L6 CT
    ## F2S4_160516_036_H01      L6 IT
    ## F2S4_160516_037_A01      L6 IT
    ## F2S4_160516_037_B01      L6 IT
    ## F2S4_160516_037_C01      L6 CT
    ## F2S4_160516_037_D01      L6 CT
    ## F2S4_160516_037_E01        Vip
    ## F2S4_160516_037_F01      L6 IT
    ## F2S4_160516_037_G01      L6 IT
    ## F2S4_160516_037_H01      L6 CT
    ## F2S4_160516_038_A01      L6 IT
    ## F2S4_160516_038_B01      L6 IT
    ## F2S4_160516_038_C01      L6 IT
    ## F2S4_160516_038_D01      L6 IT
    ## F2S4_160516_038_E01      L6 IT
    ## F2S4_160516_038_F01      L6 IT
    ## F2S4_160516_038_G01      L6 CT
    ## F2S4_160516_038_H01      L6 IT
    ## F2S4_160516_039_A01      L6 IT
    ## F2S4_160516_039_B01      L6 IT
    ## F2S4_160516_039_C01      L6 IT
    ## F2S4_160516_039_E01      L6 IT
    ## F2S4_160516_039_F01      L6 IT
    ## F2S4_160516_039_G01      L6 IT
    ## F2S4_160516_039_H01      L6 IT
    ## F2S4_160516_040_A01      L6 CT
    ## F2S4_160516_040_B01      L6 IT
    ## F2S4_160516_040_C01      L6 CT
    ## F2S4_160516_040_D01      L6 IT
    ## F2S4_160516_040_E01      L6 IT
    ## F2S4_160516_040_F01        Vip
    ## F2S4_160516_040_G01      L6 CT
    ## F2S4_160516_040_H01      L6 CT
    ## F2S4_160516_041_A01      L6 CT
    ## F2S4_160516_041_B01      L6 IT
    ## F2S4_160516_041_C01      L6 CT
    ## F2S4_160516_041_D01      L6 CT
    ## F2S4_160516_041_E01        Vip
    ## F2S4_160516_041_F01      L6 IT
    ## F2S4_160516_041_G01      L6 IT
    ## F2S4_160516_041_H01        Sst
    ## F2S4_160516_042_A01      L6 IT
    ## F2S4_160516_042_B01      L6 CT
    ## F2S4_160516_042_C01      L6 IT
    ## F2S4_160516_042_D01      L6 CT
    ## F2S4_160516_042_E01      L6 IT
    ## F2S4_160516_042_F01      L6 IT
    ## F2S4_160516_042_G01      L6 IT
    ## F2S4_160516_042_H01      L6 IT
    ## F2S4_160608_001_A01      Pvalb
    ## F2S4_160608_001_B01      Pvalb
    ## F2S4_160608_001_C01      Pvalb
    ## F2S4_160608_001_D01      Pvalb
    ## F2S4_160608_001_E01      Pvalb
    ## F2S4_160608_001_G01      Pvalb
    ## F2S4_160608_001_H01      Pvalb
    ## F2S4_160608_002_A01      Pvalb
    ## F2S4_160608_002_B01      Pvalb
    ## F2S4_160608_002_D01      Pvalb
    ## F2S4_160608_002_F01      Pvalb
    ## F2S4_160608_002_G01      Pvalb
    ## F2S4_160608_002_H01      Pvalb
    ## F2S4_160608_003_A01      Pvalb
    ## F2S4_160608_003_B01      Pvalb
    ## F2S4_160608_003_C01      Pvalb
    ## F2S4_160608_003_D01      Pvalb
    ## F2S4_160608_003_E01      Pvalb
    ## F2S4_160608_003_G01      Pvalb
    ## F2S4_160608_003_H01      Pvalb
    ## F2S4_160608_004_A01      Pvalb
    ## F2S4_160608_004_B01      Pvalb
    ## F2S4_160608_004_D01      Pvalb
    ## F2S4_160608_004_E01      Pvalb
    ## F2S4_160608_004_F01      Pvalb
    ## F2S4_160608_004_G01      Pvalb
    ## F2S4_160608_004_H01      Pvalb
    ## F2S4_160609_001_A01      Pvalb
    ## F2S4_160609_001_B01      Pvalb
    ## F2S4_160609_001_C01      Pvalb
    ## F2S4_160609_001_D01      Pvalb
    ## F2S4_160609_001_E01      Pvalb
    ## F2S4_160609_001_F01      Pvalb
    ## F2S4_160609_001_G01      Pvalb
    ## F2S4_160609_001_H01      Pvalb
    ## F2S4_160610_001_C01         L4
    ## F2S4_160610_002_A01         L4
    ## F2S4_160610_002_B01      L5 IT
    ## F2S4_160610_002_D01         L4
    ## F2S4_160610_002_E01         L4
    ## F2S4_160610_002_G01      L5 IT
    ## F2S4_160610_003_D01         L4
    ## F2S4_160610_003_E01         L4
    ## F2S4_160610_003_F01      L5 IT
    ## F2S4_160610_003_G01         L4
    ## F2S4_160610_004_C01         L4
    ## F2S4_160613_001_A01      Pvalb
    ## F2S4_160613_001_B01      L5 IT
    ## F2S4_160613_001_D01      L5 IT
    ## F2S4_160613_001_E01      L5 IT
    ## F2S4_160613_001_G01         L4
    ## F2S4_160613_002_A01         L4
    ## F2S4_160613_002_B01         L4
    ## F2S4_160613_002_D01      L5 IT
    ## F2S4_160613_002_E01         L4
    ## F2S4_160613_002_F01      L5 IT
    ## F2S4_160613_002_G01         L4
    ## F2S4_160613_003_H01         L4
    ## F2S4_160614_001_A01        Sst
    ## F2S4_160614_001_B01        Vip
    ## F2S4_160614_001_C01        Sst
    ## F2S4_160614_001_D01        Vip
    ## F2S4_160614_001_E01      Pvalb
    ## F2S4_160614_001_F01        Vip
    ## F2S4_160614_001_G01        Sst
    ## F2S4_160614_001_H01        Sst
    ## F2S4_160614_002_A01      Pvalb
    ## F2S4_160614_002_B01        Sst
    ## F2S4_160614_002_C01      Pvalb
    ## F2S4_160614_002_D01        Sst
    ## F2S4_160614_002_E01      Pvalb
    ## F2S4_160614_002_F01        Sst
    ## F2S4_160614_002_G01        Sst
    ## F2S4_160614_002_H01        Vip
    ## F2S4_160614_003_A01      Pvalb
    ## F2S4_160614_003_B01        Sst
    ## F2S4_160614_003_C01      Pvalb
    ## F2S4_160614_003_D01        Sst
    ## F2S4_160614_003_E01        Vip
    ## F2S4_160614_003_F01        Sst
    ## F2S4_160614_003_G01      Pvalb
    ## F2S4_160614_003_H01      Pvalb
    ## F2S4_160614_004_A01      Pvalb
    ## F2S4_160614_004_B01      Pvalb
    ## F2S4_160614_004_C01        Sst
    ## F2S4_160614_004_D01      Pvalb
    ## F2S4_160614_004_E01        Vip
    ## F2S4_160614_004_F01        Sst
    ## F2S4_160614_004_G01        Sst
    ## F2S4_160614_004_H01        Sst
    ## F2S4_160614_005_A01        Sst
    ## F2S4_160614_005_B01        Sst
    ## F2S4_160614_005_C01        Sst
    ## F2S4_160614_005_D01        Sst
    ## F2S4_160614_005_E01        Sst
    ## F2S4_160614_005_F01        Sst
    ## F2S4_160614_005_G01        Sst
    ## F2S4_160614_005_H01        Sst
    ## F2S4_160614_006_A01        Sst
    ## F2S4_160614_006_B01        Sst
    ## F2S4_160614_006_C01        Sst
    ## F2S4_160614_006_D01        Sst
    ## F2S4_160614_006_E01        Sst
    ## F2S4_160614_006_F01      Pvalb
    ## F2S4_160614_006_G01        Sst
    ## F2S4_160614_006_H01        Sst
    ## F2S4_160614_007_A01        Sst
    ## F2S4_160614_007_B01        Sst
    ## F2S4_160614_007_C01        Sst
    ## F2S4_160614_007_D01        Sst
    ## F2S4_160614_007_E01        Sst
    ## F2S4_160614_007_F01        Sst
    ## F2S4_160614_007_G01        Sst
    ## F2S4_160614_007_H01        Sst
    ## F2S4_160614_008_A01        Sst
    ## F2S4_160614_008_B01        Sst
    ## F2S4_160614_008_C01        Sst
    ## F2S4_160614_008_D01        Sst
    ## F2S4_160614_008_E01        Sst
    ## F2S4_160614_008_F01        Sst
    ## F2S4_160614_008_G01        Sst
    ## F2S4_160614_008_H01        Sst
    ## F2S4_160614_009_F01        Sst
    ## F2S4_160614_009_G01        Sst
    ## F2S4_160614_009_H01        Sst
    ## F2S4_160615_001_B01      L5 IT
    ## F2S4_160615_001_C01         L4
    ## F2S4_160615_001_D01      L5 IT
    ## F2S4_160615_001_E01      L5 IT
    ## F2S4_160615_001_G01      L5 IT
    ## F2S4_160615_001_H01         L4
    ## F2S4_160615_002_A01      L5 IT
    ## F2S4_160615_002_B01        Sst
    ## F2S4_160615_002_C01      L6 CT
    ## F2S4_160615_002_D01      L6 IT
    ## F2S4_160615_002_E01      L5 PT
    ## F2S4_160615_002_F01      L5 IT
    ## F2S4_160615_002_G01      L5 IT
    ## F2S4_160615_002_H01         L4
    ## F2S4_160615_003_C01      L6 IT
    ## F2S4_160615_003_D01      L5 PT
    ## F2S4_160615_003_E01      L5 IT
    ## F2S4_160615_003_F01      L5 IT
    ## F2S4_160615_003_G01      L5 IT
    ## F2S4_160615_003_H01      L5 IT
    ## F2S4_160615_004_A01      L5 IT
    ## F2S4_160615_004_B01      L5 IT
    ## F2S4_160615_004_C01      L5 IT
    ## F2S4_160615_004_D01      L5 PT
    ## F2S4_160615_004_E01      L5 IT
    ## F2S4_160615_004_F01      L5 IT
    ## F2S4_160615_004_G01      L6 CT
    ## F2S4_160615_004_H01      L5 PT
    ## F2S4_160615_005_A01      L5 IT
    ## F2S4_160615_005_B01      L5 IT
    ## F2S4_160615_005_C01      L5 PT
    ## F2S4_160615_005_D01      L6 IT
    ## F2S4_160615_005_F01      L6 IT
    ## F2S4_160615_005_G01      L5 IT
    ## F2S4_160615_005_H01      L5 PT
    ## F2S4_160615_006_A01      L5 IT
    ## F2S4_160615_006_B01      L5 IT
    ## F2S4_160615_006_C01      L5 IT
    ## F2S4_160615_006_D01      L5 IT
    ## F2S4_160615_006_E01      L6 IT
    ## F2S4_160615_006_G01      L5 IT
    ## F2S4_160615_006_H01      L5 PT
    ## F2S4_160615_007_A01      L5 IT
    ## F2S4_160615_007_B01      L5 PT
    ## F2S4_160615_007_C01      L6 IT
    ## F2S4_160615_007_D01        Sst
    ## F2S4_160615_007_E01      L5 IT
    ## F2S4_160615_007_F01      L6 IT
    ## F2S4_160615_007_G01        Sst
    ## F2S4_160615_007_H01      L5 IT
    ## F2S4_160615_008_A01      L5 IT
    ## F2S4_160615_008_B01      L5 IT
    ## F2S4_160615_008_C01      L5 IT
    ## F2S4_160615_008_D01      L5 PT
    ## F2S4_160615_008_E01      L5 PT
    ## F2S4_160615_008_F01      L5 PT
    ## F2S4_160615_008_H01      L5 IT
    ## F2S4_160615_009_A01      L5 IT
    ## F2S4_160615_009_B01      L5 IT
    ## F2S4_160615_009_C01      L5 IT
    ## F2S4_160615_009_D01      L5 IT
    ## F2S4_160615_009_E01      L5 IT
    ## F2S4_160615_009_F01      L5 IT
    ## F2S4_160615_009_G01      L5 IT
    ## F2S4_160615_009_H01      L5 IT
    ## F2S4_160615_010_A01      L5 IT
    ## F2S4_160615_010_B01      L5 IT
    ## F2S4_160615_010_C01      L5 IT
    ## F2S4_160615_010_D01        Sst
    ## F2S4_160615_010_E01        Sst
    ## F2S4_160615_010_F01      L5 IT
    ## F2S4_160615_010_G01      L5 IT
    ## F2S4_160615_010_H01      L5 IT
    ## F2S4_160615_011_A01      L5 IT
    ## F2S4_160615_011_B01      L5 IT
    ## F2S4_160615_011_C01      L5 IT
    ## F2S4_160615_011_D01      L5 IT
    ## F2S4_160615_011_E01      L6 IT
    ## F2S4_160615_011_F01      L5 IT
    ## F2S4_160615_011_G01      L5 IT
    ## F2S4_160615_011_H01      L5 IT
    ## F2S4_160615_012_A01      Pvalb
    ## F2S4_160615_012_B01      Pvalb
    ## F2S4_160615_012_C01      Pvalb
    ## F2S4_160615_012_D01      Pvalb
    ## F2S4_160615_012_E01      Pvalb
    ## F2S4_160615_012_F01      Pvalb
    ## F2S4_160615_012_G01      Pvalb
    ## F2S4_160615_012_H01      Pvalb
    ## F2S4_160615_013_A01      Pvalb
    ## F2S4_160615_013_B01      Pvalb
    ## F2S4_160615_013_C01      Pvalb
    ## F2S4_160615_013_D01      Pvalb
    ## F2S4_160615_013_E01      Pvalb
    ## F2S4_160615_013_F01      Pvalb
    ## F2S4_160615_013_G01      Pvalb
    ## F2S4_160615_013_H01      Pvalb
    ## F2S4_160615_014_A01      Pvalb
    ## F2S4_160615_014_B01      Pvalb
    ## F2S4_160615_014_C01      Pvalb
    ## F2S4_160615_014_D01      Pvalb
    ## F2S4_160615_014_E01      Pvalb
    ## F2S4_160615_014_F01      Pvalb
    ## F2S4_160615_014_G01      Pvalb
    ## F2S4_160615_014_H01      Pvalb
    ## F2S4_160615_015_A01      Pvalb
    ## F2S4_160615_015_B01      Pvalb
    ## F2S4_160616_001_A01      L5 IT
    ## F2S4_160616_001_D01      L5 IT
    ## F2S4_160616_001_E01      L5 PT
    ## F2S4_160616_001_F01      L5 IT
    ## F2S4_160616_001_G01        Sst
    ## F2S4_160616_001_H01      L5 PT
    ## F2S4_160616_002_A01      L5 IT
    ## F2S4_160616_002_B01      L5 PT
    ## F2S4_160616_002_C01        Sst
    ## F2S4_160616_002_D01      L5 PT
    ## F2S4_160616_002_E01      L5 IT
    ## F2S4_160616_002_F01      L5 IT
    ## F2S4_160616_002_G01      L5 IT
    ## F2S4_160616_002_H01      L6 IT
    ## F2S4_160616_003_A01      L5 PT
    ## F2S4_160616_003_B01      L5 IT
    ## F2S4_160616_003_C01      L5 PT
    ## F2S4_160616_003_D01      L5 IT
    ## F2S4_160616_003_E01      L5 PT
    ## F2S4_160616_003_F01      L5 PT
    ## F2S4_160616_003_G01      L5 PT
    ## F2S4_160616_003_H01      L5 IT
    ## F2S4_160616_004_A01      L5 IT
    ## F2S4_160616_004_B01      L5 PT
    ## F2S4_160616_004_C01         NP
    ## F2S4_160616_004_D01      L5 PT
    ## F2S4_160616_004_E01      L5 PT
    ## F2S4_160616_004_F01      L5 PT
    ## F2S4_160616_004_G01      L5 IT
    ## F2S4_160616_004_H01      L5 PT
    ## F2S4_160616_005_A01      L5 IT
    ## F2S4_160616_005_B01      L5 IT
    ## F2S4_160616_005_C01      L5 PT
    ## F2S4_160616_005_D01      L5 IT
    ## F2S4_160616_005_E01      L6 CT
    ## F2S4_160616_005_F01      L6 CT
    ## F2S4_160616_005_G01      L6 CT
    ## F2S4_160616_005_H01      L6 CT
    ## F2S4_160616_006_A01        Sst
    ## F2S4_160616_006_B01      L6 IT
    ## F2S4_160616_006_C01      L5 PT
    ## F2S4_160616_006_D01      L5 PT
    ## F2S4_160616_006_E01      L5 IT
    ## F2S4_160616_006_F01      L5 PT
    ## F2S4_160616_006_G01      L5 IT
    ## F2S4_160616_006_H01      L5 IT
    ## F2S4_160616_007_A01      L5 IT
    ## F2S4_160616_007_B01      L5 IT
    ## F2S4_160616_007_C01      L5 PT
    ## F2S4_160616_007_D01      L5 IT
    ## F2S4_160616_007_E01      L5 IT
    ## F2S4_160616_007_F01        Vip
    ## F2S4_160616_007_G01      L5 IT
    ## F2S4_160616_007_H01      L5 IT
    ## F2S4_160616_008_A01      L5 IT
    ## F2S4_160616_008_B01      L5 IT
    ## F2S4_160616_008_C01      L5 IT
    ## F2S4_160616_008_D01      L5 IT
    ## F2S4_160616_008_E01      L5 IT
    ## F2S4_160616_008_F01      L5 IT
    ## F2S4_160616_008_G01      L5 IT
    ## F2S4_160616_008_H01      L5 IT
    ## F2S4_160616_009_A01        Sst
    ## F2S4_160616_009_B01      L5 IT
    ## F2S4_160616_009_C01      L5 PT
    ## F2S4_160616_009_D01      L5 IT
    ## F2S4_160616_009_E01      L5 IT
    ## F2S4_160616_009_F01      L5 IT
    ## F2S4_160616_009_G01      L5 PT
    ## F2S4_160616_009_H01      L5 PT
    ## F2S4_160616_010_A01      L5 IT
    ## F2S4_160616_010_B01      L5 PT
    ## F2S4_160616_010_C01      L5 PT
    ## F2S4_160616_010_D01      L5 IT
    ## F2S4_160616_010_E01      L5 IT
    ## F2S4_160616_010_F01      L5 IT
    ## F2S4_160616_010_G01      L5 IT
    ## F2S4_160616_010_H01      L5 IT
    ## F2S4_160616_011_A01      L5 PT
    ## F2S4_160616_011_B01      L5 PT
    ## F2S4_160616_011_C01      L5 IT
    ## F2S4_160616_011_D01      Lamp5
    ## F2S4_160616_011_E01        Sst
    ## F2S4_160616_011_F01      Lamp5
    ## F2S4_160616_011_G01      Lamp5
    ## F2S4_160616_011_H01        Vip
    ## F2S4_160616_012_A01        Vip
    ## F2S4_160616_012_C01        Vip
    ## F2S4_160616_012_D01        Vip
    ## F2S4_160616_012_E01        Vip
    ## F2S4_160616_012_F01      Lamp5
    ## F2S4_160616_012_H01        Vip
    ## F2S4_160616_013_A01        Sst
    ## F2S4_160616_013_B01        Sst
    ## F2S4_160616_013_E01        Sst
    ## F2S4_160616_013_F01        Sst
    ## F2S4_160616_013_H01      Lamp5
    ## F2S4_160616_014_A01        Sst
    ## F2S4_160616_014_B01        Sst
    ## F2S4_160616_014_C01        Sst
    ## F2S4_160616_014_D01        Vip
    ## F2S4_160616_014_E01      Lamp5
    ## F2S4_160616_014_F01        Vip
    ## F2S4_160616_014_G01       Endo
    ## F2S4_160616_014_H01        Sst
    ## F2S4_160616_015_A01        Vip
    ## F2S4_160616_015_B01       Endo
    ## F2S4_160616_015_D01        Vip
    ## F2S4_160616_015_E01        Vip
    ## F2S4_160616_015_F01      Lamp5
    ## F2S4_160616_015_G01        Vip
    ## F2S4_160616_015_H01      Lamp5
    ## F2S4_160617_001_A01        L6b
    ## F2S4_160617_001_C01        L6b
    ## F2S4_160617_001_D01       Endo
    ## F2S4_160617_001_E01        L6b
    ## F2S4_160617_001_F01        L6b
    ## F2S4_160617_001_G01        L6b
    ## F2S4_160617_001_H01        L6b
    ## F2S4_160617_002_A01        L6b
    ## F2S4_160617_002_B01       Endo
    ## F2S4_160617_002_D01        L6b
    ## F2S4_160617_002_E01        L6b
    ## F2S4_160617_002_G01        L6b
    ## F2S4_160617_002_H01       Endo
    ## F2S4_160617_003_A01        L6b
    ## F2S4_160617_003_B01        L6b
    ## F2S4_160617_003_C01        L6b
    ## F2S4_160617_003_D01        L6b
    ## F2S4_160617_003_E01        L6b
    ## F2S4_160617_003_F01        L6b
    ## F2S4_160617_003_G01        L6b
    ## F2S4_160617_003_H01       Endo
    ## F2S4_160617_004_A01        L6b
    ## F2S4_160617_004_B01        L6b
    ## F2S4_160617_004_C01        L6b
    ## F2S4_160617_004_D01        L6b
    ## F2S4_160617_004_E01        L6b
    ## F2S4_160617_004_F01        L6b
    ## F2S4_160617_004_G01        L6b
    ## F2S4_160617_004_H01        L6b
    ## F2S4_160617_005_F01       Endo
    ## F2S4_160617_005_G01       Endo
    ## F2S4_160620_038_A01        L6b
    ## F2S4_160620_038_B01        L6b
    ## F2S4_160620_038_D01       VLMC
    ## F2S4_160620_038_E01        L6b
    ## F2S4_160620_038_F01        L6b
    ## F2S4_160620_038_G01        L6b
    ## F2S4_160620_038_H01        L6b
    ## F2S4_160621_001_A01      Pvalb
    ## F2S4_160621_001_B01      Pvalb
    ## F2S4_160621_001_C01      Pvalb
    ## F2S4_160621_001_D01      Pvalb
    ## F2S4_160621_001_E01      Pvalb
    ## F2S4_160621_001_F01      Pvalb
    ## F2S4_160621_001_G01      Pvalb
    ## F2S4_160621_001_H01      Pvalb
    ## F2S4_160621_002_A01      Pvalb
    ## F2S4_160621_002_B01      Pvalb
    ## F2S4_160621_002_C01      Pvalb
    ## F2S4_160621_002_D01      Pvalb
    ## F2S4_160621_002_E01      Pvalb
    ## F2S4_160621_002_F01      Pvalb
    ## F2S4_160621_002_G01      Pvalb
    ## F2S4_160621_002_H01        Sst
    ## F2S4_160621_003_A01      Pvalb
    ## F2S4_160621_003_B01      Pvalb
    ## F2S4_160621_003_C01      Pvalb
    ## F2S4_160621_003_D01      Pvalb
    ## F2S4_160621_003_E01      Pvalb
    ## F2S4_160621_003_F01      Pvalb
    ## F2S4_160621_003_G01      Pvalb
    ## F2S4_160621_003_H01      Pvalb
    ## F2S4_160621_004_A01      Pvalb
    ## F2S4_160621_004_B01      Pvalb
    ## F2S4_160621_004_C01      Pvalb
    ## F2S4_160621_004_E01      Pvalb
    ## F2S4_160621_004_F01      Pvalb
    ## F2S4_160621_004_G01        Sst
    ## F2S4_160621_005_A01      Pvalb
    ## F2S4_160621_005_B01      Pvalb
    ## F2S4_160621_005_C01      Pvalb
    ## F2S4_160621_005_D01      Pvalb
    ## F2S4_160621_005_E01      Pvalb
    ## F2S4_160621_005_F01      Pvalb
    ## F2S4_160621_005_G01        Sst
    ## F2S4_160621_005_H01      Pvalb
    ## F2S4_160621_006_A01      Pvalb
    ## F2S4_160621_006_B01      Pvalb
    ## F2S4_160621_006_C01      Pvalb
    ## F2S4_160621_006_D01      Pvalb
    ## F2S4_160621_006_E01      Pvalb
    ## F2S4_160621_006_F01      Pvalb
    ## F2S4_160621_006_G01      Pvalb
    ## F2S4_160621_006_H01      Pvalb
    ## F2S4_160621_007_A01      Pvalb
    ## F2S4_160621_007_E01        Sst
    ## F2S4_160621_007_F01      Pvalb
    ## F2S4_160621_007_G01      Pvalb
    ## F2S4_160621_007_H01      Pvalb
    ## F2S4_160622_001_A01      L5 PT
    ## F2S4_160622_001_B01      L5 IT
    ## F2S4_160622_001_C01      L5 IT
    ## F2S4_160622_001_D01      L5 IT
    ## F2S4_160622_001_E01      L5 IT
    ## F2S4_160622_001_F01      L5 PT
    ## F2S4_160622_001_G01      L5 IT
    ## F2S4_160622_001_H01      L5 PT
    ## F2S4_160622_002_A01      L5 IT
    ## F2S4_160622_002_B01      L5 PT
    ## F2S4_160622_002_C01      L5 IT
    ## F2S4_160622_002_D01      L5 IT
    ## F2S4_160622_002_E01      L5 IT
    ## F2S4_160622_002_F01      L5 PT
    ## F2S4_160622_002_G01      L5 PT
    ## F2S4_160622_002_H01      L6 IT
    ## F2S4_160622_003_A01      L5 IT
    ## F2S4_160622_003_B01      L5 IT
    ## F2S4_160622_003_C01      L5 IT
    ## F2S4_160622_003_D01      L5 IT
    ## F2S4_160622_003_E01      L5 PT
    ## F2S4_160622_003_F01      L5 IT
    ## F2S4_160622_003_G01      L5 PT
    ## F2S4_160622_003_H01      L5 IT
    ## F2S4_160622_004_A01      L5 PT
    ## F2S4_160622_004_B01      L6 IT
    ## F2S4_160622_004_C01      L5 PT
    ## F2S4_160622_004_D01      L5 PT
    ## F2S4_160622_004_E01      L5 IT
    ## F2S4_160622_004_F01      L5 IT
    ## F2S4_160622_004_G01      L5 IT
    ## F2S4_160622_004_H01      L5 PT
    ## F2S4_160622_005_A01      L5 IT
    ## F2S4_160622_005_B01      L6 IT
    ## F2S4_160622_005_C01      L5 PT
    ## F2S4_160622_005_D01      L5 PT
    ## F2S4_160622_005_E01      L5 PT
    ## F2S4_160622_005_G01      L5 IT
    ## F2S4_160622_005_H01      L6 IT
    ## F2S4_160622_006_A01      L5 IT
    ## F2S4_160622_006_B01      L5 IT
    ## F2S4_160622_006_C01         L4
    ## F2S4_160622_006_D01      L5 IT
    ## F2S4_160622_006_E01      L5 IT
    ## F2S4_160622_006_F01      L5 IT
    ## F2S4_160622_006_G01      L5 IT
    ## F2S4_160622_006_H01      L5 IT
    ## F2S4_160622_007_A01      L5 PT
    ## F2S4_160622_007_B01      L5 PT
    ## F2S4_160622_007_C01      L5 PT
    ## F2S4_160622_007_E01      L5 IT
    ## F2S4_160622_007_F01      L5 PT
    ## F2S4_160622_007_G01      L5 IT
    ## F2S4_160622_007_H01      L5 IT
    ## F2S4_160622_008_A01      L5 PT
    ## F2S4_160622_008_B01      L5 PT
    ## F2S4_160622_008_C01      L5 IT
    ## F2S4_160622_008_D01      L5 IT
    ## F2S4_160622_008_E01      L5 IT
    ## F2S4_160622_008_F01      L5 IT
    ## F2S4_160622_008_G01      L5 IT
    ## F2S4_160622_009_B01      L5 PT
    ## F2S4_160622_009_C01      L5 IT
    ## F2S4_160622_009_D01      L5 PT
    ## F2S4_160622_009_E01      L5 PT
    ## F2S4_160622_009_F01      L5 PT
    ## F2S4_160622_009_G01      L5 PT
    ## F2S4_160622_009_H01      L5 IT
    ## F2S4_160622_010_A01      L5 IT
    ## F2S4_160622_010_B01      L5 PT
    ## F2S4_160622_010_C01      L5 IT
    ## F2S4_160622_010_D01      L5 IT
    ## F2S4_160622_010_E01      L5 PT
    ## F2S4_160622_010_F01      L5 IT
    ## F2S4_160622_010_G01      L5 IT
    ## F2S4_160622_010_H01      L5 IT
    ## F2S4_160622_011_A01      L5 IT
    ## F2S4_160622_011_B01      L5 IT
    ## F2S4_160622_011_C01      L5 PT
    ## F2S4_160622_011_D01      L5 PT
    ## F2S4_160622_011_E01      L5 PT
    ## F2S4_160622_011_F01      L5 IT
    ## F2S4_160622_011_G01      L5 PT
    ## F2S4_160622_011_H01        Sst
    ## F2S4_160622_012_A01      L5 PT
    ## F2S4_160622_012_B01      L5 IT
    ## F2S4_160622_012_C01      L5 IT
    ## F2S4_160622_012_D01      L5 PT
    ## F2S4_160622_012_E01      L5 PT
    ## F2S4_160622_012_F01      L5 IT
    ## F2S4_160622_012_G01      L5 IT
    ## F2S4_160622_012_H01      L5 IT
    ## F2S4_160622_013_A01      L5 IT
    ## F2S4_160622_013_B01      L5 PT
    ## F2S4_160622_013_C01      L5 IT
    ## F2S4_160622_013_D01      L5 IT
    ## F2S4_160622_013_E01      L5 IT
    ## F2S4_160622_013_F01      L5 PT
    ## F2S4_160622_013_G01      L5 PT
    ## F2S4_160622_013_H01      L5 IT
    ## F2S4_160622_014_A01      L5 IT
    ## F2S4_160622_014_B01      L5 PT
    ## F2S4_160622_014_C01      L5 IT
    ## F2S4_160622_014_D01      L5 IT
    ## F2S4_160622_014_E01      L5 PT
    ## F2S4_160622_014_F01      L5 PT
    ## F2S4_160622_015_A01      L5 PT
    ## F2S4_160622_015_B01      L5 PT
    ## F2S4_160622_015_C01      L5 IT
    ## F2S4_160622_015_D01      L5 PT
    ## F2S4_160622_015_E01      L5 IT
    ## F2S4_160622_015_F01      L5 IT
    ## F2S4_160622_015_G01      L5 IT
    ## F2S4_160622_016_B01      L5 IT
    ## F2S4_160622_016_C01      L5 PT
    ## F2S4_160622_016_D01      L5 PT
    ## F2S4_160622_016_E01      L5 IT
    ## F2S4_160622_016_F01      L5 PT
    ## F2S4_160622_016_G01      L5 PT
    ## F2S4_160622_016_H01      L5 IT
    ## F2S4_160622_017_A01      L5 IT
    ## F2S4_160622_017_B01      L5 IT
    ## F2S4_160622_017_C01      L5 PT
    ## F2S4_160622_017_E01      L5 PT
    ## F2S4_160622_017_H01      L6 IT
    ## F2S4_160622_018_A01      L6 CT
    ## F2S4_160622_018_B01      L6 IT
    ## F2S4_160622_018_C01      L6 IT
    ## F2S4_160622_018_D01      L5 PT
    ## F2S4_160622_018_E01      L5 PT
    ## F2S4_160622_018_F01      L5 PT
    ## F2S4_160622_019_A01      L5 IT
    ## F2S4_160622_019_B01      L5 IT
    ## F2S4_160622_019_C01      L5 IT
    ## F2S4_160622_019_E01      L5 IT
    ## F2S4_160622_019_F01      L5 IT
    ## F2S4_160622_019_G01      L5 PT
    ## F2S4_160622_019_H01      L5 PT
    ## F2S4_160622_020_A01      L5 IT
    ## F2S4_160622_020_C01      L5 IT
    ## F2S4_160622_020_D01      L5 IT
    ## F2S4_160622_020_E01      L5 IT
    ## F2S4_160622_020_F01      L5 IT
    ## F2S4_160622_020_G01      L5 IT
    ## F2S4_160622_020_H01      L5 IT
    ## F2S4_160622_021_A01      L5 IT
    ## F2S4_160622_021_B01      L5 PT
    ## F2S4_160622_021_C01      L5 PT
    ## F2S4_160622_021_F01      L5 PT
    ## F2S4_160622_021_G01      L5 IT
    ## F2S4_160622_021_H01      L5 IT
    ## F2S4_160622_022_A01      L6 IT
    ## F2S4_160622_022_B01      L5 PT
    ## F2S4_160622_022_D01      L5 PT
    ## F2S4_160622_022_E01      L5 IT
    ## F2S4_160622_022_F01      L5 PT
    ## F2S4_160622_022_G01      L5 IT
    ## F2S4_160622_022_H01      L5 IT
    ## F2S4_160623_034_A01      L5 IT
    ## F2S4_160623_034_B01      L5 IT
    ## F2S4_160623_034_C01      L5 IT
    ## F2S4_160623_034_D01      L5 IT
    ## F2S4_160623_035_A01      L5 IT
    ## F2S4_160623_035_B01      L5 PT
    ## F2S4_160623_035_C01      L5 IT
    ## F2S4_160623_035_D01      L5 IT
    ## F2S4_160623_035_E01      L5 IT
    ## F2S4_160623_035_F01      L5 IT
    ## F2S4_160623_035_G01      L5 IT
    ## F2S4_160623_035_H01      L5 IT
    ## F2S4_160623_036_A01      L5 IT
    ## F2S4_160623_036_B01      L5 IT
    ## F2S4_160623_036_C01      L5 PT
    ## F2S4_160623_036_E01      L5 PT
    ## F2S4_160623_036_F01      L5 PT
    ## F2S4_160623_036_H01      L5 IT
    ## F2S4_160623_037_A01      L5 IT
    ## F2S4_160623_037_C01      L5 PT
    ## F2S4_160623_037_D01      L5 PT
    ## F2S4_160623_037_E01      L5 IT
    ## F2S4_160623_037_F01      L5 PT
    ## F2S4_160623_037_G01      L5 IT
    ## F2S4_160623_037_H01      L5 PT
    ## F2S4_160623_038_A01      L5 PT
    ## F2S4_160623_038_B01      L5 IT
    ## F2S4_160623_038_C01      L5 IT
    ## F2S4_160623_038_D01      L5 PT
    ## F2S4_160623_038_F01      L5 IT
    ## F2S4_160623_038_G01      L5 PT
    ## F2S4_160623_039_A01      L5 IT
    ## F2S4_160623_039_C01      L5 PT
    ## F2S4_160623_039_D01      L5 IT
    ## F2S4_160623_039_E01      L5 IT
    ## F2S4_160623_039_F01      L5 IT
    ## F2S4_160623_039_G01      L5 IT
    ## F2S4_160623_039_H01      L5 PT
    ## F2S4_160623_040_A01      L5 PT
    ## F2S4_160623_040_B01      L5 PT
    ## F2S4_160623_040_D01      L5 IT
    ## F2S4_160623_040_F01        Sst
    ## F2S4_160623_040_G01      L5 PT
    ## F2S4_160623_041_A01      L5 PT
    ## F2S4_160623_041_B01      L5 IT
    ## F2S4_160623_041_C01      L5 IT
    ## F2S4_160623_041_D01      L5 IT
    ## F2S4_160623_041_E01      L5 IT
    ## F2S4_160623_041_F01         NP
    ## F2S4_160623_041_G01      L5 PT
    ## F2S4_160623_041_H01      L5 IT
    ## F2S4_160623_042_A01        Sst
    ## F2S4_160623_042_B01      L5 IT
    ## F2S4_160623_042_C01      L5 IT
    ## F2S4_160623_042_D01      L5 IT
    ## F2S4_160623_042_E01      L5 IT
    ## F2S4_160623_042_F01      L5 IT
    ## F2S4_160623_042_G01      L5 PT
    ## F2S4_160623_043_A01      L5 PT
    ## F2S4_160623_043_B01      L5 IT
    ## F2S4_160623_043_C01      L5 PT
    ## F2S4_160623_043_D01      L5 IT
    ## F2S4_160623_043_E01      L5 IT
    ## F2S4_160623_043_F01      L5 PT
    ## F2S4_160623_043_H01      L5 IT
    ## F2S4_160623_044_B01      L5 IT
    ## F2S4_160623_044_C01      L5 PT
    ## F2S4_160623_044_D01      L5 IT
    ## F2S4_160623_044_E01      L5 IT
    ## F2S4_160623_044_F01      L5 IT
    ## F2S4_160623_044_G01         L4
    ## F2S4_160623_044_H01      L5 IT
    ## F2S4_160624_021_B01      L5 IT
    ## F2S4_160624_021_C01      L5 IT
    ## F2S4_160624_021_D01      L5 IT
    ## F2S4_160624_021_E01        Sst
    ## F2S4_160624_021_F01      L5 IT
    ## F2S4_160624_021_G01      L5 IT
    ## F2S4_160624_021_H01      L5 IT
    ## F2S4_160624_022_A01      L5 IT
    ## F2S4_160624_022_B01      L5 IT
    ## F2S4_160624_022_C01      L5 IT
    ## F2S4_160624_022_D01      L5 IT
    ## F2S4_160624_022_E01      L5 IT
    ## F2S4_160624_022_F01      L5 IT
    ## F2S4_160624_022_G01      L5 IT
    ## F2S4_160624_022_H01      L5 IT
    ## F2S4_160624_023_A01      L5 IT
    ## F2S4_160624_023_B01      L5 PT
    ## F2S4_160624_023_C01      L5 IT
    ## F2S4_160624_023_E01      L5 IT
    ## F2S4_160624_023_F01      L5 IT
    ## F2S4_160624_023_G01      L5 IT
    ## F2S4_160624_023_H01      L5 IT
    ## F2S4_160624_024_A01         NP
    ## F2S4_160624_024_B01      L5 PT
    ## F2S4_160624_024_C01        Sst
    ## F2S4_160624_024_D01      L5 PT
    ## F2S4_160624_024_E01      L5 PT
    ## F2S4_160624_024_F01      L5 IT
    ## F2S4_160624_024_G01      L5 IT
    ## F2S4_160624_024_H01        Sst
    ## F2S4_160624_025_A01      L5 IT
    ## F2S4_160624_025_C01      L5 IT
    ## F2S4_160624_025_D01      L5 PT
    ## F2S4_160624_025_E01      L5 IT
    ## F2S4_160624_025_F01      L5 IT
    ## F2S4_160624_025_G01      L5 IT
    ## F2S4_160624_025_H01      L5 IT
    ## F2S4_160624_026_A01      L5 PT
    ## F2S4_160624_026_B01      L5 PT
    ## F2S4_160624_026_C01      L5 IT
    ## F2S4_160624_026_D01         NP
    ## F2S4_160624_026_E01      L5 IT
    ## F2S4_160624_026_F01        Sst
    ## F2S4_160624_026_G01      L5 IT
    ## F2S4_160624_027_A01      L5 IT
    ## F2S4_160624_027_B01      L5 PT
    ## F2S4_160624_027_C01      L5 PT
    ## F2S4_160624_027_E01      L5 PT
    ## F2S4_160624_027_F01      L5 PT
    ## F2S4_160624_027_G01      L5 PT
    ## F2S4_160624_027_H01      L5 IT
    ## F2S4_160624_028_A01      L5 IT
    ## F2S4_160624_028_B01      L5 IT
    ## F2S4_160624_028_C01      L5 IT
    ## F2S4_160624_028_D01      L5 PT
    ## F2S4_160624_028_E01      L5 IT
    ## F2S4_160624_028_F01      L5 PT
    ## F2S4_160624_028_G01      L5 IT
    ## F2S4_160624_028_H01      L5 IT
    ## F2S4_160624_029_A01      L5 IT
    ## F2S4_160624_029_B01      L5 IT
    ## F2S4_160624_029_C01      L5 IT
    ## F2S4_160624_029_D01      L5 PT
    ## F2S4_160624_029_E01      L6 IT
    ## F2S4_160624_029_F01      L5 IT
    ## F2S4_160624_029_G01      L5 PT
    ## F2S4_160624_029_H01      L5 IT
    ## F2S4_160624_030_A01      L5 IT
    ## F2S4_160624_030_B01      L5 PT
    ## F2S4_160624_030_C01      L5 IT
    ## F2S4_160624_030_D01      L5 IT
    ## F2S4_160624_030_E01      L5 IT
    ## F2S4_160624_030_F01      L5 IT
    ## F2S4_160624_030_G01      L5 IT
    ## F2S4_160624_030_H01      L5 IT
    ## F2S4_160624_031_A01      L5 IT
    ## F2S4_160624_031_B01      L5 IT
    ## F2S4_160624_031_C01      L5 IT
    ## F2S4_160624_031_D01      L5 IT
    ## F2S4_160624_031_E01      L5 IT
    ## F2S4_160624_031_F01      L5 IT
    ## F2S4_160624_031_G01      L5 IT
    ## F2S4_160624_031_H01      L6 CT
    ## F2S4_160624_032_A01      L6 CT
    ## F2S4_160624_032_B01        Sst
    ## F2S4_160624_032_C01        Sst
    ## F2S4_160624_032_D01        L6b
    ## F2S4_160627_013_F01      Pvalb
    ## F2S4_160627_013_G01      Pvalb
    ## F2S4_160627_013_H01      Pvalb
    ## F2S4_160627_014_A01      Pvalb
    ## F2S4_160627_014_B01      Pvalb
    ## F2S4_160627_014_C01      Pvalb
    ## F2S4_160627_014_D01      Pvalb
    ## F2S4_160627_014_E01      Pvalb
    ## F2S4_160627_014_F01      Pvalb
    ## F2S4_160627_014_G01      Pvalb
    ## F2S4_160627_014_H01      Pvalb
    ## F2S4_160627_015_A01      Pvalb
    ## F2S4_160627_015_B01      Pvalb
    ## F2S4_160627_015_C01      Pvalb
    ## F2S4_160627_015_D01      Pvalb
    ## F2S4_160627_015_E01      Pvalb
    ## F2S4_160627_015_F01      Pvalb
    ## F2S4_160627_015_G01      Pvalb
    ## F2S4_160627_015_H01      Pvalb
    ## F2S4_160627_016_A01      Pvalb
    ## F2S4_160627_016_B01      Pvalb
    ## F2S4_160627_016_C01      Pvalb
    ## F2S4_160627_016_D01      Pvalb
    ## F2S4_160627_016_E01      Pvalb
    ## F2S4_160627_016_F01      Pvalb
    ## F2S4_160627_016_G01      Pvalb
    ## F2S4_160627_016_H01      Pvalb
    ## F2S4_160627_017_A01      Pvalb
    ## F2S4_160627_017_B01      Pvalb
    ## F2S4_160627_017_C01      Pvalb
    ## F2S4_160627_017_D01      Pvalb
    ## F2S4_160627_017_E01      Pvalb
    ## F2S4_160627_017_F01      Pvalb
    ## F2S4_160627_017_G01      Pvalb
    ## F2S4_160627_017_H01      Pvalb
    ## F2S4_160701_001_A01    L2/3 IT
    ## F2S4_160701_001_B01    L2/3 IT
    ## F2S4_160701_001_C01         L4
    ## F2S4_160701_001_D01    L2/3 IT
    ## F2S4_160701_001_F01    L2/3 IT
    ## F2S4_160701_001_G01    L2/3 IT
    ## F2S4_160701_001_H01    L2/3 IT
    ## F2S4_160701_002_A01    L2/3 IT
    ## F2S4_160701_002_B01    L2/3 IT
    ## F2S4_160701_002_C01    L2/3 IT
    ## F2S4_160701_002_D01    L2/3 IT
    ## F2S4_160701_002_E01         L4
    ## F2S4_160701_002_F01         L4
    ## F2S4_160701_002_G01    L2/3 IT
    ## F2S4_160701_002_H01    L2/3 IT
    ## F2S4_160701_003_A01    L2/3 IT
    ## F2S4_160701_003_B01    L2/3 IT
    ## F2S4_160701_003_C01    L2/3 IT
    ## F2S4_160701_003_D01    L2/3 IT
    ## F2S4_160701_003_E01    L2/3 IT
    ## F2S4_160701_003_F01    L2/3 IT
    ## F2S4_160701_003_G01    L2/3 IT
    ## F2S4_160701_003_H01    L2/3 IT
    ## F2S4_160701_004_A01    L2/3 IT
    ## F2S4_160701_004_B01    L2/3 IT
    ## F2S4_160701_004_C01    L2/3 IT
    ## F2S4_160701_004_D01    L2/3 IT
    ## F2S4_160701_004_E01    L2/3 IT
    ## F2S4_160701_004_F01         L4
    ## F2S4_160701_004_G01    L2/3 IT
    ## F2S4_160701_004_H01    L2/3 IT
    ## F2S4_160701_005_A01    L2/3 IT
    ## F2S4_160701_005_B01    L2/3 IT
    ## F2S4_160701_005_C01    L2/3 IT
    ## F2S4_160701_005_D01    L2/3 IT
    ## F2S4_160701_005_E01         L4
    ## F2S4_160701_005_F01    L2/3 IT
    ## F2S4_160701_005_G01    L2/3 IT
    ## F2S4_160701_005_H01    L2/3 IT
    ## F2S4_160701_006_A01    L2/3 IT
    ## F2S4_160701_006_B01    L2/3 IT
    ## F2S4_160701_006_C01    L2/3 IT
    ## F2S4_160701_006_D01    L2/3 IT
    ## F2S4_160701_006_E01    L2/3 IT
    ## F2S4_160701_006_F01         L4
    ## F2S4_160701_006_G01    L2/3 IT
    ## F2S4_160701_006_H01    L2/3 IT
    ## F2S4_160701_007_A01         L4
    ## F2S4_160701_007_B01    L2/3 IT
    ## F2S4_160701_007_C01    L2/3 IT
    ## F2S4_160701_007_D01    L2/3 IT
    ## F2S4_160701_007_E01         L4
    ## F2S4_160701_007_F01         L4
    ## F2S4_160701_007_G01    L2/3 IT
    ## F2S4_160701_007_H01    L2/3 IT
    ## F2S4_160701_008_A01    L2/3 IT
    ## F2S4_160701_008_B01    L2/3 IT
    ## F2S4_160701_008_C01    L2/3 IT
    ## F2S4_160701_008_D01         L4
    ## F2S4_160701_008_E01    L2/3 IT
    ## F2S4_160701_008_F01    L2/3 IT
    ## F2S4_160701_008_G01    L2/3 IT
    ## F2S4_160701_008_H01    L2/3 IT
    ## F2S4_160701_009_A01    L2/3 IT
    ## F2S4_160701_009_B01    L2/3 IT
    ## F2S4_160701_009_C01    L2/3 IT
    ## F2S4_160701_009_D01    L2/3 IT
    ## F2S4_160701_009_E01    L2/3 IT
    ## F2S4_160701_009_F01    L2/3 IT
    ## F2S4_160701_009_G01    L2/3 IT
    ## F2S4_160701_009_H01    L2/3 IT
    ## F2S4_160701_010_A01         L4
    ## F2S4_160701_010_B01         L4
    ## F2S4_160701_010_C01         L4
    ## F2S4_160701_010_D01         L4
    ## F2S4_160701_010_E01         L4
    ## F2S4_160701_010_F01         L4
    ## F2S4_160701_010_G01         L4
    ## F2S4_160701_010_H01         L4
    ## F2S4_160701_011_A01         L4
    ## F2S4_160701_011_B01         L4
    ## F2S4_160701_011_C01         L4
    ## F2S4_160701_011_D01      L5 IT
    ## F2S4_160701_011_E01         L4
    ## F2S4_160701_011_F01         L4
    ## F2S4_160701_011_G01         L4
    ## F2S4_160701_011_H01      L5 IT
    ## F2S4_160701_012_A01         L4
    ## F2S4_160701_012_B01         L4
    ## F2S4_160701_012_C01         L4
    ## F2S4_160701_012_D01         L4
    ## F2S4_160701_012_E01         L4
    ## F2S4_160701_012_F01         L4
    ## F2S4_160701_012_G01         L4
    ## F2S4_160701_012_H01         L4
    ## F2S4_160701_013_A01      L5 IT
    ## F2S4_160701_013_B01         L4
    ## F2S4_160701_013_C01         L4
    ## F2S4_160701_013_D01         L4
    ## F2S4_160701_013_E01         L4
    ## F2S4_160701_013_F01         L4
    ## F2S4_160701_013_G01         L4
    ## F2S4_160701_013_H01         L4
    ## F2S4_160701_014_A01      L5 IT
    ## F2S4_160701_014_B01         L4
    ## F2S4_160701_014_C01         L4
    ## F2S4_160701_014_D01         L4
    ## F2S4_160701_014_E01         L4
    ## F2S4_160701_014_F01         L4
    ## F2S4_160701_014_G01      L5 IT
    ## F2S4_160701_015_A01         L4
    ## F2S4_160701_015_B01         L4
    ## F2S4_160701_015_C01         L4
    ## F2S4_160701_015_D01         L4
    ## F2S4_160701_015_E01         L4
    ## F2S4_160701_015_F01         L4
    ## F2S4_160701_015_G01      L5 IT
    ## F2S4_160701_015_H01         L4
    ## F2S4_160701_016_A01         L4
    ## F2S4_160701_016_B01         L4
    ## F2S4_160701_016_C01      L5 IT
    ## F2S4_160701_016_D01         L4
    ## F2S4_160701_016_E01      L5 IT
    ## F2S4_160701_016_F01         L4
    ## F2S4_160701_016_G01         L4
    ## F2S4_160701_016_H01         L4
    ## F2S4_160701_017_A01      L5 IT
    ## F2S4_160701_017_B01         L4
    ## F2S4_160701_017_C01      L5 IT
    ## F2S4_160701_017_D01         L4
    ## F2S4_160701_017_E01         L4
    ## F2S4_160701_017_F01         L4
    ## F2S4_160701_017_H01         L4
    ## F2S4_160701_018_A01      L6 IT
    ## F2S4_160701_018_B01      L6 IT
    ## F2S4_160701_018_C01      L6 IT
    ## F2S4_160701_018_D01      L6 CT
    ## F2S4_160701_018_E01      L6 IT
    ## F2S4_160701_018_F01      L6 IT
    ## F2S4_160701_018_G01      L6 IT
    ## F2S4_160701_018_H01      L6 CT
    ## F2S4_160701_019_A01      L6 CT
    ## F2S4_160701_019_B01      L6 CT
    ## F2S4_160701_019_C01      L6 IT
    ## F2S4_160701_019_D01      L6 CT
    ## F2S4_160701_019_E01      L6 IT
    ## F2S4_160701_019_F01      L5 IT
    ## F2S4_160701_019_G01      L6 IT
    ## F2S4_160701_019_H01      L6 IT
    ## F2S4_160701_020_A01      L6 CT
    ## F2S4_160701_020_B01      L6 CT
    ## F2S4_160701_020_C01      L6 CT
    ## F2S4_160701_020_D01      L6 IT
    ## F2S4_160701_020_E01      L6 IT
    ## F2S4_160701_020_F01      L6 IT
    ## F2S4_160701_020_G01      L6 IT
    ## F2S4_160701_020_H01      L6 IT
    ## F2S4_160701_021_A01      L6 CT
    ## F2S4_160701_021_B01      L6 IT
    ## F2S4_160701_021_C01      L6 CT
    ## F2S4_160701_021_D01      L6 IT
    ## F2S4_160701_021_E01      L6 CT
    ## F2S4_160701_021_F01      L6 IT
    ## F2S4_160701_021_G01      L6 IT
    ## F2S4_160701_021_H01      L6 IT
    ## F2S4_160701_022_A01      L6 CT
    ## F2S4_160701_022_B01      L6 IT
    ## F2S4_160701_022_C01      L6 IT
    ## F2S4_160701_022_D01      L6 IT
    ## F2S4_160701_022_E01      L6 IT
    ## F2S4_160701_022_F01      L6 IT
    ## F2S4_160701_022_G01      L6 IT
    ## F2S4_160701_022_H01      L6 IT
    ## F2S4_160701_023_A01      L6 IT
    ## F2S4_160701_023_B01      L6 IT
    ## F2S4_160701_023_C01      L6 CT
    ## F2S4_160701_023_D01      L6 CT
    ## F2S4_160701_023_E01      L6 CT
    ## F2S4_160701_023_F01      L6 CT
    ## F2S4_160701_023_G01      L6 CT
    ## F2S4_160701_023_H01      L6 CT
    ## F2S4_160701_024_A01      L6 IT
    ## F2S4_160701_024_B01      L6 CT
    ## F2S4_160701_024_C01      L6 CT
    ## F2S4_160701_024_E01      L6 IT
    ## F2S4_160701_024_G01      L6 IT
    ## F2S4_160701_024_H01      L6 IT
    ## F2S4_160701_025_A01      L6 CT
    ## F2S4_160701_025_B01      L6 IT
    ## F2S4_160701_025_C01      L6 IT
    ## F2S4_160701_025_D01      L6 CT
    ## F2S4_160701_025_E01      L6 IT
    ## F2S4_160701_025_F01      L6 IT
    ## F2S4_160701_025_G01      L6 CT
    ## F2S4_160701_025_H01      L6 CT
    ## F2S4_160701_026_A01      L6 IT
    ## F2S4_160701_026_B01      L6 CT
    ## F2S4_160701_026_C01      L6 IT
    ## F2S4_160701_026_D01      L6 IT
    ## F2S4_160701_026_E01      L6 CT
    ## F2S4_160701_026_H01      L6 CT
    ## F2S4_160701_027_A01      L6 IT
    ## F2S4_160701_027_B01      L6 IT
    ## F2S4_160701_027_C01      L6 IT
    ## F2S4_160701_027_D01      L6 IT
    ## F2S4_160701_027_E01      L6 IT
    ## F2S4_160701_027_F01      L6 CT
    ## F2S4_160701_027_G01      L6 IT
    ## F2S4_160701_027_H01      L6 CT
    ## F2S4_160701_028_A01      L6 IT
    ## F2S4_160701_028_B01      L6 CT
    ## F2S4_160701_028_C01      L6 IT
    ## F2S4_160701_028_D01      L6 CT
    ## F2S4_160701_028_E01      L6 IT
    ## F2S4_160701_028_F01      L6 IT
    ## F2S4_160701_028_G01      L6 IT
    ## F2S4_160701_028_H01      L6 IT
    ## F2S4_160701_029_A01      L6 IT
    ## F2S4_160701_029_B01      L6 IT
    ## F2S4_160701_029_C01      L6 IT
    ## F2S4_160701_029_D01         L4
    ## F2S4_160701_029_E01      L6 IT
    ## F2S4_160701_029_F01      L6 IT
    ## F2S4_160701_029_G01      L6 CT
    ## F2S4_160701_029_H01      L6 CT
    ## F2S4_160701_030_A01      L6 IT
    ## F2S4_160701_030_B01      L6 CT
    ## F2S4_160701_030_C01      L6 CT
    ## F2S4_160701_030_D01      L6 CT
    ## F2S4_160701_030_E01      L6 CT
    ## F2S4_160701_030_F01      L6 IT
    ## F2S4_160701_030_G01      L6 IT
    ## F2S4_160701_030_H01      L6 IT
    ## F2S4_160701_031_A01      L6 IT
    ## F2S4_160701_031_B01      L6 CT
    ## F2S4_160701_031_C01      L6 IT
    ## F2S4_160701_031_D01      L6 CT
    ## F2S4_160701_031_E01      L6 IT
    ## F2S4_160701_031_F01      L6 CT
    ## F2S4_160701_031_G01      L6 CT
    ## F2S4_160701_031_H01      L6 IT
    ## F2S4_160701_032_A01         L4
    ## F2S4_160701_032_B01         L4
    ## F2S4_160701_032_D01         L4
    ## F2S4_160701_032_E01         L4
    ## F2S4_160701_032_G01         L4
    ## F2S4_160701_032_H01         L4
    ## F2S4_160701_033_A01         L4
    ## F2S4_160701_033_B01         L4
    ## F2S4_160701_033_C01         L4
    ## F2S4_160701_033_D01         L4
    ## F2S4_160701_033_E01         L4
    ## F2S4_160701_033_F01         L4
    ## F2S4_160701_033_G01         L4
    ## F2S4_160701_033_H01         L4
    ## F2S4_160701_034_A01         L4
    ## F2S4_160701_034_B01         L4
    ## F2S4_160701_034_C01         L4
    ## F2S4_160701_034_D01         L4
    ## F2S4_160701_034_E01      L5 PT
    ## F2S4_160701_034_F01         L4
    ## F2S4_160701_034_G01         L4
    ## F2S4_160701_034_H01         L4
    ## F2S4_160701_035_A01         L4
    ## F2S4_160701_035_C01         L4
    ## F2S4_160701_035_D01         L4
    ## F2S4_160701_035_E01         L4
    ## F2S4_160701_035_F01         L4
    ## F2S4_160701_035_G01         L4
    ## F2S4_160701_035_H01         L4
    ## F2S4_160701_036_B01         L4
    ## F2S4_160701_036_C01         L4
    ## F2S4_160701_036_D01         L4
    ## F2S4_160701_036_G01         L4
    ## F2S4_160701_036_H01         L4
    ## F2S4_160701_037_A01         L4
    ## F2S4_160701_037_B01         L4
    ## F2S4_160701_037_C01         L4
    ## F2S4_160701_037_D01      L5 IT
    ## F2S4_160701_037_E01      L5 IT
    ## F2S4_160701_037_F01         L4
    ## F2S4_160701_037_H01         L4
    ## F2S4_160701_038_A01         L4
    ## F2S4_160701_038_B01         L4
    ## F2S4_160701_038_C01         L4
    ## F2S4_160701_038_D01         L4
    ## F2S4_160701_038_E01         L4
    ## F2S4_160701_038_G01         L4
    ## F2S4_160701_038_H01         L4
    ## F2S4_160701_039_A01         L4
    ## F2S4_160701_039_B01         L4
    ## F2S4_160701_039_C01         L4
    ## F2S4_160701_039_D01         L4
    ## F2S4_160701_039_E01         L4
    ## F2S4_160701_039_F01         L4
    ## F2S4_160701_039_G01         L4
    ## F2S4_160701_039_H01         L4
    ## F2S4_160701_040_A01    L2/3 IT
    ## F2S4_160701_040_B01    L2/3 IT
    ## F2S4_160701_040_C01    L2/3 IT
    ## F2S4_160701_040_D01    L2/3 IT
    ## F2S4_160701_040_E01    L2/3 IT
    ## F2S4_160701_040_F01    L2/3 IT
    ## F2S4_160701_040_G01    L2/3 IT
    ## F2S4_160701_040_H01    L2/3 IT
    ## F2S4_160701_041_A01    L2/3 IT
    ## F2S4_160701_041_B01    L2/3 IT
    ## F2S4_160701_041_C01    L2/3 IT
    ## F2S4_160701_041_D01    L2/3 IT
    ## F2S4_160701_041_E01    L2/3 IT
    ## F2S4_160701_041_F01    L2/3 IT
    ## F2S4_160701_041_G01    L2/3 IT
    ## F2S4_160701_041_H01    L2/3 IT
    ## F2S4_160701_042_A01    L2/3 IT
    ## F2S4_160701_042_B01    L2/3 IT
    ## F2S4_160701_042_C01    L2/3 IT
    ## F2S4_160701_042_D01    L2/3 IT
    ## F2S4_160701_042_E01    L2/3 IT
    ## F2S4_160701_042_F01    L2/3 IT
    ## F2S4_160701_042_G01    L2/3 IT
    ## F2S4_160701_042_H01    L2/3 IT
    ## F2S4_160701_043_A01    L2/3 IT
    ## F2S4_160701_043_B01    L2/3 IT
    ## F2S4_160701_043_C01    L2/3 IT
    ## F2S4_160701_043_D01    L2/3 IT
    ## F2S4_160701_043_E01    L2/3 IT
    ## F2S4_160701_043_F01    L2/3 IT
    ## F2S4_160701_043_G01    L2/3 IT
    ## F2S4_160701_043_H01    L2/3 IT
    ## F2S4_160701_044_B01    L2/3 IT
    ## F2S4_160701_044_C01    L2/3 IT
    ## F2S4_160701_044_D01    L2/3 IT
    ## F2S4_160701_044_E01    L2/3 IT
    ## F2S4_160701_044_F01    L2/3 IT
    ## F2S4_160701_044_G01    L2/3 IT
    ## F2S4_160701_044_H01    L2/3 IT
    ## F2S4_160701_045_A01    L2/3 IT
    ## F2S4_160701_045_B01    L2/3 IT
    ## F2S4_160701_045_C01    L2/3 IT
    ## F2S4_160701_045_D01    L2/3 IT
    ## F2S4_160701_045_E01    L2/3 IT
    ## F2S4_160701_045_F01    L2/3 IT
    ## F2S4_160701_045_G01         L4
    ## F2S4_160701_045_H01    L2/3 IT
    ## F2S4_160701_046_A01    L2/3 IT
    ## F2S4_160701_046_B01    L2/3 IT
    ## F2S4_160701_046_C01    L2/3 IT
    ## F2S4_160701_046_D01    L2/3 IT
    ## F2S4_160701_046_E01    L2/3 IT
    ## F2S4_160701_046_F01    L2/3 IT
    ## F2S4_160701_046_G01    L2/3 IT
    ## F2S4_160701_046_H01    L2/3 IT
    ## F2S4_160701_047_A01    L2/3 IT
    ## F2S4_160701_047_B01    L2/3 IT
    ## F2S4_160701_047_C01    L2/3 IT
    ## F2S4_160701_047_D01    L2/3 IT
    ## F2S4_160701_047_E01    L2/3 IT
    ## F2S4_160701_047_F01    L2/3 IT
    ## F2S4_160701_047_G01    L2/3 IT
    ## F2S4_160701_047_H01    L2/3 IT
    ## F2S4_160705_001_A01        L6b
    ## F2S4_160705_001_B01        L6b
    ## F2S4_160705_001_C01        L6b
    ## F2S4_160705_001_F01        L6b
    ## F2S4_160705_001_H01      L6 IT
    ## F2S4_160705_002_B01       VLMC
    ## F2S4_160705_002_C01        L6b
    ## F2S4_160706_001_A01      L6 IT
    ## F2S4_160706_001_B01      L6 CT
    ## F2S4_160706_001_C01      L6 IT
    ## F2S4_160706_001_D01        L6b
    ## F2S4_160706_001_E01      L6 IT
    ## F2S4_160706_001_F01      L6 IT
    ## F2S4_160706_001_G01      L6 IT
    ## F2S4_160706_001_H01      L6 IT
    ## F2S4_160706_002_A01      L6 IT
    ## F2S4_160706_002_B01      L6 IT
    ## F2S4_160706_002_C01      L6 IT
    ## F2S4_160706_002_D01      L6 IT
    ## F2S4_160706_002_E01      L6 IT
    ## F2S4_160706_002_F01      L6 IT
    ## F2S4_160706_002_G01      L6 IT
    ## F2S4_160706_002_H01      L6 CT
    ## F2S4_160706_003_A01      L6 CT
    ## F2S4_160706_003_B01      L6 IT
    ## F2S4_160706_003_C01      L6 CT
    ## F2S4_160706_003_D01      L6 IT
    ## F2S4_160706_003_E01      L6 IT
    ## F2S4_160706_003_F01      L6 CT
    ## F2S4_160706_003_G01      L6 IT
    ## F2S4_160706_003_H01      L6 IT
    ## F2S4_160706_004_A01      L6 IT
    ## F2S4_160706_004_B01      L6 IT
    ## F2S4_160706_004_C01      L6 IT
    ## F2S4_160706_004_D01      L6 CT
    ## F2S4_160706_004_E01      L6 CT
    ## F2S4_160706_004_F01      L6 CT
    ## F2S4_160706_004_G01      L6 IT
    ## F2S4_160706_004_H01      L6 IT
    ## F2S4_160706_005_A01      L6 CT
    ## F2S4_160706_005_B01      L6 IT
    ## F2S4_160706_005_D01      L6 IT
    ## F2S4_160706_005_E01      L6 IT
    ## F2S4_160706_005_F01      L6 IT
    ## F2S4_160706_005_G01      L6 CT
    ## F2S4_160706_005_H01      L6 IT
    ## F2S4_160706_006_A01      L6 IT
    ## F2S4_160706_006_B01      L6 IT
    ## F2S4_160706_006_C01      L6 IT
    ## F2S4_160706_006_D01      L6 IT
    ## F2S4_160706_006_E01      L6 IT
    ## F2S4_160706_006_F01      L6 IT
    ## F2S4_160706_006_G01      L6 IT
    ## F2S4_160706_006_H01      L6 IT
    ## F2S4_160706_007_A01      L6 IT
    ## F2S4_160706_007_B01      L6 IT
    ## F2S4_160706_007_C01      L6 IT
    ## F2S4_160706_007_D01      L6 CT
    ## F2S4_160706_007_E01      L6 CT
    ## F2S4_160706_007_F01      L6 IT
    ## F2S4_160706_007_G01      L6 CT
    ## F2S4_160706_007_H01      L6 IT
    ## F2S4_160706_008_A01      L6 IT
    ## F2S4_160706_008_B01      L6 CT
    ## F2S4_160706_008_C01      L6 IT
    ## F2S4_160706_008_D01      L6 CT
    ## F2S4_160706_008_E01      L6 IT
    ## F2S4_160706_008_F01      L6 IT
    ## F2S4_160706_008_G01      L6 CT
    ## F2S4_160706_008_H01      L6 CT
    ## F2S4_160706_009_A01      L6 IT
    ## F2S4_160706_009_B01      L6 IT
    ## F2S4_160706_009_C01      L6 IT
    ## F2S4_160706_009_D01      L6 CT
    ## F2S4_160706_009_E01      L6 IT
    ## F2S4_160706_009_F01      L6 IT
    ## F2S4_160706_009_G01      L6 CT
    ## F2S4_160706_009_H01      L6 IT
    ## F2S4_160706_010_A01      L6 IT
    ## F2S4_160706_010_B01      L6 IT
    ## F2S4_160706_010_C01      L6 IT
    ## F2S4_160706_010_D01      L6 IT
    ## F2S4_160706_010_E01      L6 IT
    ## F2S4_160706_010_F01      L6 IT
    ## F2S4_160706_010_G01      L6 IT
    ## F2S4_160706_010_H01      L6 IT
    ## F2S4_160706_011_A01      L6 CT
    ## F2S4_160706_011_B01      L6 IT
    ## F2S4_160706_011_D01      L6 IT
    ## F2S4_160706_011_E01      L6 IT
    ## F2S4_160706_011_F01      L6 IT
    ## F2S4_160706_011_G01      L6 IT
    ## F2S4_160706_011_H01      L6 IT
    ## F2S4_160706_012_A01      L6 IT
    ## F2S4_160706_012_B01      L6 IT
    ## F2S4_160706_012_C01      L6 IT
    ## F2S4_160706_012_D01      L6 IT
    ## F2S4_160706_012_E01      L6 IT
    ## F2S4_160706_012_F01      L6 IT
    ## F2S4_160706_012_G01      L6 IT
    ## F2S4_160706_013_A01      L6 CT
    ## F2S4_160706_013_B01      L6 IT
    ## F2S4_160706_013_C01      L6 IT
    ## F2S4_160706_013_D01      L6 IT
    ## F2S4_160706_013_E01      L6 CT
    ## F2S4_160706_013_F01      L6 CT
    ## F2S4_160706_013_G01      L6 IT
    ## F2S4_160706_013_H01      L6 IT
    ## F2S4_160706_014_A01      L6 CT
    ## F2S4_160706_014_B01      L6 IT
    ## F2S4_160706_014_C01      L6 CT
    ## F2S4_160706_014_E01      L6 IT
    ## F2S4_160706_014_F01      L6 IT
    ## F2S4_160706_014_G01      L6 IT
    ## F2S4_160706_014_H01      L6 IT
    ## F2S4_160706_015_A01      L6 IT
    ## F2S4_160706_015_B01      L6 IT
    ## F2S4_160706_015_C01      L6 CT
    ## F2S4_160706_015_D01      L6 IT
    ## F2S4_160706_015_E01      L6 CT
    ## F2S4_160706_015_F01      L6 IT
    ## F2S4_160706_015_G01      L6 IT
    ## F2S4_160706_015_H01      L6 IT
    ## F2S4_160706_016_A01      L6 IT
    ## F2S4_160706_016_B01      L6 IT
    ## F2S4_160706_016_D01      L6 IT
    ## F2S4_160706_016_E01      L6 IT
    ## F2S4_160706_016_F01      L6 IT
    ## F2S4_160706_016_G01      L6 IT
    ## F2S4_160706_016_H01      L6 IT
    ## F2S4_160706_017_A01      L6 CT
    ## F2S4_160706_017_B01      L6 IT
    ## F2S4_160706_017_C01      L6 IT
    ## F2S4_160706_017_D01      L6 CT
    ## F2S4_160706_017_E01      L6 CT
    ## F2S4_160706_017_F01      L6 IT
    ## F2S4_160706_017_G01      L6 IT
    ## F2S4_160706_017_H01      L6 IT
    ## F2S4_160706_018_A01      L6 IT
    ## F2S4_160706_018_B01      L6 IT
    ## F2S4_160706_018_C01      L6 IT
    ## F2S4_160706_018_D01      L6 IT
    ## F2S4_160706_018_E01      L6 IT
    ## F2S4_160706_018_F01      L6 IT
    ## F2S4_160706_018_G01      L6 IT
    ## F2S4_160706_018_H01      L6 IT
    ## F2S4_160706_019_A01      L6 CT
    ## F2S4_160706_019_B01      L6 IT
    ## F2S4_160706_019_C01      L6 CT
    ## F2S4_160706_019_D01      L6 IT
    ## F2S4_160706_019_E01      L6 IT
    ## F2S4_160706_019_F01      L6 IT
    ## F2S4_160706_019_G01      L6 IT
    ## F2S4_160706_019_H01      L6 IT
    ## F2S4_160706_020_A01      L6 CT
    ## F2S4_160706_020_B01      L6 CT
    ## F2S4_160706_020_C01      L6 IT
    ## F2S4_160706_020_D01      L6 IT
    ## F2S4_160706_020_E01      L6 IT
    ## F2S4_160706_020_F01      L6 IT
    ## F2S4_160706_020_G01      L6 IT
    ## F2S4_160706_020_H01      L6 CT
    ## F2S4_160706_021_A01      L6 IT
    ## F2S4_160706_021_B01      L6 CT
    ## F2S4_160706_021_C01      L6 CT
    ## F2S4_160706_021_D01      L6 IT
    ## F2S4_160706_021_E01      L6 IT
    ## F2S4_160706_021_F01      L6 IT
    ## F2S4_160706_021_G01      L6 IT
    ## F2S4_160706_021_H01      L6 IT
    ## F2S4_160706_022_A01      L6 CT
    ## F2S4_160706_022_B01      L6 CT
    ## F2S4_160706_022_C01      L6 IT
    ## F2S4_160706_022_D01      L6 CT
    ## F2S4_160706_022_E01      L6 CT
    ## F2S4_160706_022_F01      L6 CT
    ## F2S4_160706_022_G01      L6 IT
    ## F2S4_160706_022_H01      L6 CT
    ## F2S4_160706_023_A01      L6 IT
    ## F2S4_160706_023_B01      L6 CT
    ## F2S4_160706_023_C01      L6 IT
    ## F2S4_160706_023_D01      L6 CT
    ## F2S4_160706_023_E01      L6 CT
    ## F2S4_160706_023_F01      L6 IT
    ## F2S4_160706_023_G01      L6 IT
    ## F2S4_160706_023_H01      L6 IT
    ## F2S4_160706_024_A01      L6 IT
    ## F2S4_160706_024_B01      L6 IT
    ## F2S4_160706_024_C01      L6 IT
    ## F2S4_160706_024_D01      L6 IT
    ## F2S4_160706_024_E01      L6 IT
    ## F2S4_160706_024_F01      L6 IT
    ## F2S4_160706_024_G01      L6 IT
    ## F2S4_160706_024_H01      L6 CT
    ## F2S4_160706_025_A01      L6 CT
    ## F2S4_160706_025_B01      L6 CT
    ## F2S4_160706_025_C01      L6 IT
    ## F2S4_160706_025_D01      L6 IT
    ## F2S4_160706_025_E01      L6 CT
    ## F2S4_160706_025_F01      L6 IT
    ## F2S4_160706_025_G01      L6 CT
    ## F2S4_160706_025_H01      L6 IT
    ## F2S4_160706_026_A01        L6b
    ## F2S4_160706_026_B01      L6 IT
    ## F2S4_160706_026_C01      L6 CT
    ## F2S4_160706_026_D01      L6 CT
    ## F2S4_160706_026_E01      L6 IT
    ## F2S4_160706_026_F01      L6 IT
    ## F2S4_160706_026_G01      L6 CT
    ## F2S4_160706_026_H01      L6 IT
    ## F2S4_160706_027_A01      L6 CT
    ## F2S4_160706_027_B01      L6 IT
    ## F2S4_160706_027_C01      L6 CT
    ## F2S4_160706_027_D01      L6 CT
    ## F2S4_160706_027_E01      L6 CT
    ## F2S4_160706_027_F01      L6 IT
    ## F2S4_160706_027_G01      L6 CT
    ## F2S4_160706_027_H01      L6 IT
    ## F2S4_160706_028_A01      L6 IT
    ## F2S4_160706_028_B01      L6 IT
    ## F2S4_160706_028_C01      L6 IT
    ## F2S4_160706_028_D01      L6 IT
    ## F2S4_160706_028_E01      L6 IT
    ## F2S4_160706_028_F01      L6 IT
    ## F2S4_160706_028_G01      L6 CT
    ## F2S4_160706_028_H01      L6 IT
    ## F2S4_160706_029_A01      L6 IT
    ## F2S4_160706_029_B01      L6 CT
    ## F2S4_160706_029_C01      L6 IT
    ## F2S4_160706_029_D01      L6 IT
    ## F2S4_160706_029_E01      L6 IT
    ## F2S4_160706_029_F01      L6 IT
    ## F2S4_160706_029_G01      L6 IT
    ## F2S4_160706_029_H01      L6 IT
    ## F2S4_160706_030_A01      L6 CT
    ## F2S4_160706_030_B01      L6 CT
    ## F2S4_160706_030_C01      L6 IT
    ## F2S4_160706_030_D01      L6 IT
    ## F2S4_160706_030_E01      L6 IT
    ## F2S4_160706_030_F01      L6 IT
    ## F2S4_160706_030_G01      L6 IT
    ## F2S4_160706_030_H01      L6 IT
    ## F2S4_160706_031_B01      L6 CT
    ## F2S4_160706_031_C01      L6 IT
    ## F2S4_160706_031_D01      L6 IT
    ## F2S4_160706_031_E01      L6 IT
    ## F2S4_160706_031_F01      L6 CT
    ## F2S4_160706_031_G01      L6 IT
    ## F2S4_160706_031_H01      L6 CT
    ## F2S4_160706_032_A01      L6 IT
    ## F2S4_160706_032_B01      L6 CT
    ## F2S4_160706_032_C01      L6 IT
    ## F2S4_160706_032_D01      L6 IT
    ## F2S4_160706_032_E01      L6 IT
    ## F2S4_160706_032_F01      L6 IT
    ## F2S4_160706_032_G01      L6 IT
    ## F2S4_160706_032_H01      L6 IT
    ## F2S4_160706_033_A01      L6 IT
    ## F2S4_160706_033_B01      L6 IT
    ## F2S4_160706_033_C01      L6 IT
    ## F2S4_160706_033_D01      L6 IT
    ## F2S4_160706_033_E01      L6 IT
    ## F2S4_160706_033_F01      L6 CT
    ## F2S4_160706_033_G01      L6 CT
    ## F2S4_160706_033_H01      L6 CT
    ## F2S4_160706_034_A01      L6 CT
    ## F2S4_160706_034_B01      L6 CT
    ## F2S4_160706_034_E01      L6 IT
    ## F2S4_160706_034_F01      L6 IT
    ## F2S4_160706_034_G01      L6 CT
    ## F2S4_160706_034_H01      L6 IT
    ## F2S4_160706_035_A01      L6 IT
    ## F2S4_160706_035_B01      L6 CT
    ## F2S4_160706_035_C01      L6 IT
    ## F2S4_160706_035_D01      L6 IT
    ## F2S4_160706_035_E01      L6 IT
    ## F2S4_160706_035_F01      L6 CT
    ## F2S4_160706_035_G01      L6 IT
    ## F2S4_160706_035_H01      L6 IT
    ## F2S4_160706_036_A01      L6 IT
    ## F2S4_160706_036_B01      L6 IT
    ## F2S4_160706_036_C01      L6 IT
    ## F2S4_160706_036_D01      L6 IT
    ## F2S4_160706_036_E01      L6 IT
    ## F2S4_160706_036_F01      L6 CT
    ## F2S4_160706_036_G01      L6 CT
    ## F2S4_160706_036_H01      L6 IT
    ## F2S4_160707_007_F01        L6b
    ## F2S4_160707_007_G01        L6b
    ## F2S4_160707_007_H01        L6b
    ## F2S4_160707_008_A01        L6b
    ## F2S4_160707_008_B01        L6b
    ## F2S4_160707_008_C01        L6b
    ## F2S4_160707_008_D01        L6b
    ## F2S4_160707_008_E01        L6b
    ## F2S4_160707_008_F01        L6b
    ## F2S4_160707_008_H01        L6b
    ## F2S4_160707_009_A01        L6b
    ## F2S4_160707_009_B01       VLMC
    ## F2S4_160707_009_C01        L6b
    ## F2S4_160707_009_D01        L6b
    ## F2S4_160707_009_E01        L6b
    ## F2S4_160707_009_F01        L6b
    ## F2S4_160707_009_G01        L6b
    ## F2S4_160707_009_H01       Endo
    ## F2S4_160707_010_A01        L6b
    ## F2S4_160707_010_B01      L6 CT
    ## F2S4_160707_010_C01        L6b
    ## F2S4_160707_010_D01        L6b
    ## F2S4_160707_010_E01        L6b
    ## F2S4_160707_010_F01       Endo
    ## F2S4_160707_010_G01        L6b
    ## F2S4_160707_010_H01       Endo
    ## F2S4_160707_011_E01        L6b
    ## F2S4_160707_011_F01        L6b
    ## F2S4_160707_011_G01        L6b
    ## F2S4_160707_011_H01        L6b
    ## F2S4_160708_006_G01      Pvalb
    ## F2S4_160708_006_H01      Pvalb
    ## F2S4_160712_001_A01      Pvalb
    ## F2S4_160712_001_B01      Pvalb
    ## F2S4_160712_001_C01      Pvalb
    ## F2S4_160712_001_D01      Pvalb
    ## F2S4_160712_001_E01      Pvalb
    ## F2S4_160712_001_F01      Pvalb
    ## F2S4_160712_001_G01      Pvalb
    ## F2S4_160712_001_H01      Pvalb
    ## F2S4_160712_002_A01      Pvalb
    ## F2S4_160712_002_B01      Pvalb
    ## F2S4_160712_002_C01      Pvalb
    ## F2S4_160712_002_D01      Pvalb
    ## F2S4_160712_002_E01      Pvalb
    ## F2S4_160712_002_F01      Pvalb
    ## F2S4_160712_002_G01      Pvalb
    ## F2S4_160712_002_H01      Pvalb
    ## F2S4_160713_001_A01         L4
    ## F2S4_160713_001_B01         L4
    ## F2S4_160713_001_C01         L4
    ## F2S4_160713_001_D01         L4
    ## F2S4_160713_001_E01         L4
    ## F2S4_160713_001_F01         L4
    ## F2S4_160713_001_G01         L4
    ## F2S4_160713_001_H01         L4
    ## F2S4_160713_002_A01         L4
    ## F2S4_160713_002_B01         L4
    ## F2S4_160713_002_C01         L4
    ## F2S4_160713_002_D01         L4
    ## F2S4_160713_002_E01         L4
    ## F2S4_160713_002_F01         L4
    ## F2S4_160713_002_G01         L4
    ## F2S4_160713_002_H01         L4
    ## F2S4_160713_003_A01         L4
    ## F2S4_160713_003_B01         L4
    ## F2S4_160713_003_C01         L4
    ## F2S4_160713_003_D01         L4
    ## F2S4_160713_003_E01         L4
    ## F2S4_160713_003_F01         L4
    ## F2S4_160713_003_G01         L4
    ## F2S4_160713_003_H01         L4
    ## F2S4_160713_004_A01         L4
    ## F2S4_160713_004_C01         L4
    ## F2S4_160713_004_D01         L4
    ## F2S4_160713_004_E01         L4
    ## F2S4_160713_004_F01         L4
    ## F2S4_160713_004_G01         L4
    ## F2S4_160713_004_H01         L4
    ## F2S4_160713_005_A01         L4
    ## F2S4_160713_005_B01         L4
    ## F2S4_160713_005_C01         L4
    ## F2S4_160713_005_D01         L4
    ## F2S4_160713_005_E01         L4
    ## F2S4_160713_005_F01         L4
    ## F2S4_160713_005_G01         L4
    ## F2S4_160713_005_H01         L4
    ## F2S4_160713_006_A01      L5 IT
    ## F2S4_160713_006_B01         L4
    ## F2S4_160713_006_C01         L4
    ## F2S4_160713_006_D01         L4
    ## F2S4_160713_006_E01         L4
    ## F2S4_160713_006_F01         L4
    ## F2S4_160713_006_G01         L4
    ## F2S4_160713_006_H01         L4
    ## F2S4_160713_007_A01       VLMC
    ## F2S4_160713_007_B01        L6b
    ## F2S4_160713_007_C01        L6b
    ## F2S4_160713_007_D01        L6b
    ## F2S4_160713_007_E01        L6b
    ## F2S4_160713_007_F01        L6b
    ## F2S4_160713_007_G01       VLMC
    ## F2S4_160713_007_H01         L4
    ## F2S4_160713_008_A01        L6b
    ## F2S4_160713_008_B01       Endo
    ## F2S4_160713_008_C01        L6b
    ## F2S4_160713_008_D01        L6b
    ## F2S4_160714_001_A01        L6b
    ## F2S4_160714_001_B01        L6b
    ## F2S4_160714_001_C01        L6b
    ## F2S4_160714_001_D01        L6b
    ## F2S4_160714_001_E01        L6b
    ## F2S4_160714_001_F01       Endo
    ## F2S4_160714_001_G01        L6b
    ## F2S4_160714_001_H01       VLMC
    ## F2S4_160714_002_A01        L6b
    ## F2S4_160714_002_B01        L6b
    ## F2S4_160714_002_C01        L6b
    ## F2S4_160714_002_D01        L6b
    ## F2S4_160714_002_E01        L6b
    ## F2S4_160714_002_F01        Vip
    ## F2S4_160714_002_G01        Vip
    ## F2S4_160714_002_H01        Vip
    ## F2S4_160714_003_A01        Vip
    ## F2S4_160714_003_B01        Vip
    ## F2S4_160714_003_D01        Vip
    ## F2S4_160714_003_E01        Vip
    ## F2S4_160715_001_B01        L6b
    ## F2S4_160715_001_C01        L6b
    ## F2S4_160715_001_D01       Endo
    ## F2S4_160715_001_F01       Endo
    ## F2S4_160715_001_G01       Endo
    ## F2S4_160715_001_H01       Endo
    ## F2S4_160715_002_A01        L6b
    ## F2S4_160715_002_B01       Endo
    ## F2S4_160715_002_C01       Endo
    ## F2S4_160715_002_D01       Endo
    ## F2S4_160715_002_E01       Endo
    ## F2S4_160715_002_F01        L6b
    ## F2S4_160715_002_G01       Endo
    ## F2S4_160715_002_H01       VLMC
    ## F2S4_160715_003_B01        L6b
    ## F2S4_160715_003_D01       Endo
    ## F2S4_160715_003_E01        L6b
    ## F2S4_160715_003_F01        L6b
    ## F2S4_160715_003_G01        L6b
    ## F2S4_160715_003_H01        L6b
    ## F2S4_160715_004_A01       Endo
    ## F2S4_160715_004_B01       Endo
    ## F2S4_160715_004_C01       Endo
    ## F2S4_160715_004_D01        L6b
    ## F2S4_160715_004_F01       Endo
    ## F2S4_160715_004_G01       Endo
    ## F2S4_160715_004_H01       Endo
    ## F2S4_160715_005_A01        L6b
    ## F2S4_160715_005_B01        L6b
    ## F2S4_160715_005_D01       Endo
    ## F2S4_160715_005_E01        L6b
    ## F2S4_160715_005_F01       Endo
    ## F2S4_160715_005_G01       Endo
    ## F2S4_160715_005_H01       Endo
    ## F2S4_160715_006_A01        L6b
    ## F2S4_160715_006_C01        L6b
    ## F2S4_160715_006_D01        L6b
    ## F2S4_160715_006_E01        L6b
    ## F2S4_160715_006_F01        L6b
    ## F2S4_160715_006_G01       Endo
    ## F2S4_160715_006_H01        L6b
    ## F2S4_160715_007_A01       Endo
    ## F2S4_160715_007_B01        L6b
    ## F2S4_160715_007_D01        L6b
    ## F2S4_160715_007_E01       Endo
    ## F2S4_160715_007_G01        L6b
    ## F2S4_160715_007_H01       Endo
    ## F2S4_160715_008_A01       Endo
    ## F2S4_160715_008_B01       Endo
    ## F2S4_160715_008_C01        L6b
    ## F2S4_160715_008_D01       Endo
    ## F2S4_160715_008_E01       Endo
    ## F2S4_160715_008_F01       Endo
    ## F2S4_160715_008_G01       Endo
    ## F2S4_160715_008_H01       Endo
    ## F2S4_160715_009_A01       Endo
    ## F2S4_160715_009_B01        L6b
    ## F2S4_160715_009_C01       Endo
    ## F2S4_160715_009_D01        Vip
    ## F2S4_160715_009_E01        Vip
    ## F2S4_160715_009_F01        Vip
    ## F2S4_160715_009_G01        Vip
    ## F2S4_160715_009_H01        Vip
    ## F2S4_160715_010_A01        Vip
    ## F2S4_160715_010_B01        Vip
    ## F2S4_160715_010_C01        Vip
    ## F2S4_160715_010_D01        Vip
    ## F2S4_160715_010_F01        Vip
    ## F2S4_160715_010_G01        Vip
    ## F2S4_160718_004_B01        L6b
    ## F2S4_160718_004_C01        L6b
    ## F2S4_160718_004_D01        L6b
    ## F2S4_160718_004_E01        L6b
    ## F2S4_160718_004_F01        L6b
    ## F2S4_160718_004_G01        L6b
    ## F2S4_160718_004_H01        L6b
    ## F2S4_160718_005_A01        L6b
    ## F2S4_160718_005_B01        L6b
    ## F2S4_160718_005_C01        L6b
    ## F2S4_160718_005_E01        L6b
    ## F2S4_160718_005_F01        L6b
    ## F2S4_160718_005_G01        L6b
    ## F2S4_160718_005_H01        L6b
    ## F2S4_160718_006_A01        L6b
    ## F2S4_160718_006_B01       Endo
    ## F2S4_160718_006_C01       Endo
    ## F2S4_160718_006_D01        L6b
    ## F2S4_160718_006_E01        L6b
    ## F2S4_160718_006_F01        L6b
    ## F2S4_160718_006_G01        L6b
    ## F2S4_160718_006_H01       Endo
    ## F2S4_160718_007_A01       Endo
    ## F2S4_160718_007_B01        L6b
    ## F2S4_160718_007_C01       Endo
    ## F2S4_160718_007_D01        L6b
    ## F2S4_160718_007_E01        L6b
    ## F2S4_160718_007_F01        L6b
    ## F2S4_160718_007_G01       Endo
    ## F2S4_160718_007_H01        L6b
    ## F2S4_160718_008_A01        L6b
    ## F2S4_160718_008_B01        L6b
    ## F2S4_160718_008_D01        L6b
    ## F2S4_160718_008_E01        L6b
    ## F2S4_160718_008_F01        L6b
    ## F2S4_160718_008_G01        L6b
    ## F2S4_160719_001_B01         L4
    ## F2S4_160719_001_C01         L4
    ## F2S4_160719_001_D01         L4
    ## F2S4_160719_001_E01         L4
    ## F2S4_160719_001_F01         L4
    ## F2S4_160719_001_G01         L4
    ## F2S4_160719_001_H01         L4
    ## F2S4_160719_002_A01         L4
    ## F2S4_160719_002_B01         L4
    ## F2S4_160719_002_C01         L4
    ## F2S4_160719_002_D01         L4
    ## F2S4_160719_002_E01         L4
    ## F2S4_160719_002_F01         L4
    ## F2S4_160719_002_G01         L4
    ## F2S4_160719_002_H01         L4
    ## F2S4_160719_003_A01         L4
    ## F2S4_160719_003_B01         L4
    ## F2S4_160719_003_C01         L4
    ## F2S4_160719_003_D01         L4
    ## F2S4_160719_003_E01         L4
    ## F2S4_160719_003_F01         L4
    ## F2S4_160719_003_G01         L4
    ## F2S4_160719_003_H01         L4
    ## F2S4_160719_004_A01         L4
    ## F2S4_160719_004_B01         L4
    ## F2S4_160719_004_C01         L4
    ## F2S4_160719_004_D01         L4
    ## F2S4_160719_004_E01         L4
    ## F2S4_160719_004_F01         L4
    ## F2S4_160719_004_G01         L4
    ## F2S4_160719_004_H01         L4
    ## F2S4_160719_005_A01         L4
    ## F2S4_160719_005_B01         L4
    ## F2S4_160719_005_C01         L4
    ## F2S4_160719_005_D01         L4
    ## F2S4_160719_005_E01         L4
    ## F2S4_160719_005_F01         L4
    ## F2S4_160719_005_G01         L4
    ## F2S4_160719_005_H01         L4
    ## F2S4_160719_006_E01         L4
    ## F2S4_160719_006_F01         L4
    ## F2S4_160719_006_G01         L4
    ## F2S4_160719_006_H01         L4
    ## F2S4_160719_007_A01         L4
    ## F2S4_160719_007_B01         L4
    ## F2S4_160719_007_C01         L4
    ## F2S4_160719_007_E01         L4
    ## F2S4_160719_007_F01         L4
    ## F2S4_160719_007_G01         L4
    ## F2S4_160719_007_H01         L4
    ## F2S4_160719_008_A01         L4
    ## F2S4_160719_008_B01      L5 IT
    ## F2S4_160719_008_C01         L4
    ## F2S4_160719_008_D01         L4
    ## F2S4_160719_008_E01         L4
    ## F2S4_160719_008_F01         L4
    ## F2S4_160719_008_G01         L4
    ## F2S4_160719_008_H01         L4
    ## F2S4_160719_009_A01         L4
    ## F2S4_160719_009_B01         L4
    ## F2S4_160719_009_C01         L4
    ## F2S4_160719_009_D01         L4
    ## F2S4_160719_009_E01         L4
    ## F2S4_160719_009_F01         L4
    ## F2S4_160719_009_G01         L4
    ## F2S4_160719_009_H01         L4
    ## F2S4_160719_010_A01         L4
    ## F2S4_160719_010_B01         L4
    ## F2S4_160719_010_D01         L4
    ## F2S4_160719_010_E01         L4
    ## F2S4_160719_010_F01         L4
    ## F2S4_160719_010_G01         L4
    ## F2S4_160719_010_H01         L4
    ## F2S4_160719_011_A01         L4
    ## F2S4_160719_011_B01         L4
    ## F2S4_160719_011_C01         L4
    ## F2S4_160719_011_D01         L4
    ## F2S4_160719_011_E01         L4
    ## F2S4_160719_011_F01         L4
    ## F2S4_160719_011_G01         L4
    ## F2S4_160719_011_H01         L4
    ## F2S4_160719_012_A01         L4
    ## F2S4_160719_012_B01         L4
    ## F2S4_160719_012_C01         L4
    ## F2S4_160719_012_D01         L4
    ## F2S4_160719_012_E01         L4
    ## F2S4_160719_012_F01         L4
    ## F2S4_160719_012_G01         L4
    ## F2S4_160719_012_H01         L4
    ## F2S4_160719_013_A01         L4
    ## F2S4_160719_013_B01         L4
    ## F2S4_160719_013_C01         L4
    ## F2S4_160719_013_D01         L4
    ## F2S4_160719_013_E01         L4
    ## F2S4_160719_013_F01      L5 IT
    ## F2S4_160719_013_G01         L4
    ## F2S4_160719_013_H01         L4
    ## F2S4_160720_003_A01       Endo
    ## F2S4_160720_003_B01        L6b
    ## F2S4_160720_003_C01       Endo
    ## F2S4_160720_003_E01        L6b
    ## F2S4_160720_003_F01        L6b
    ## F2S4_160720_003_G01        L6b
    ## F2S4_160720_003_H01        L6b
    ## F2S4_160725_002_G01      Pvalb
    ## F2S4_160725_002_H01        Sst
    ## F2S4_160725_003_A01        Sst
    ## F2S4_160725_003_E01      Pvalb
    ## F2S4_160725_003_F01      Pvalb
    ## F2S4_160725_004_A01      Pvalb
    ## F2S4_160727_002_A01      Pvalb
    ## F2S4_160727_002_B01        Sst
    ## F2S4_160727_002_D01        Sst
    ## F2S4_160727_002_E01        Sst
    ## F2S4_160727_002_F01      Pvalb
    ## F2S4_160727_002_G01      Pvalb
    ## F2S4_160727_003_A01        Sst
    ## F2S4_160727_003_B01        Sst
    ## F2S4_160727_003_C01        Sst
    ## F2S4_160727_003_D01        Sst
    ## F2S4_160727_003_E01        Sst
    ## F2S4_160727_003_F01        Sst
    ## F2S4_160727_003_G01        Sst
    ## F2S4_160727_004_A01        Sst
    ## F2S4_160727_004_B01      Pvalb
    ## F2S4_160727_004_C01        Sst
    ## F2S4_160727_004_D01      Pvalb
    ## F2S4_160727_004_E01        Sst
    ## F2S4_160727_004_F01        Sst
    ## F2S4_160727_004_G01        Sst
    ## F2S4_160727_004_H01        Sst
    ## F2S4_160727_005_A01      Pvalb
    ## F2S4_160802_001_A01      Lamp5
    ## F2S4_160802_001_B01      Lamp5
    ## F2S4_160802_001_C01        Vip
    ## F2S4_160802_001_D01        Sst
    ## F2S4_160802_001_E01        Vip
    ## F2S4_160802_001_F01      Lamp5
    ## F2S4_160802_001_G01      Lamp5
    ## F2S4_160802_001_H01      Lamp5
    ## F2S4_160802_002_A01      Lamp5
    ## F2S4_160802_002_B01      Lamp5
    ## F2S4_160802_002_C01      Lamp5
    ## F2S4_160802_002_D01   Serpinf1
    ## F2S4_160802_002_E01      Lamp5
    ## F2S4_160802_002_F01      Lamp5
    ## F2S4_160802_002_G01      Lamp5
    ## F2S4_160802_002_H01        Vip
    ## F2S4_160802_003_A01        Vip
    ## F2S4_160802_003_B01        Sst
    ## F2S4_160802_003_C01        Sst
    ## F2S4_160802_003_D01        Sst
    ## F2S4_160802_003_E01        Sst
    ## F2S4_160802_003_F01        Sst
    ## F2S4_160802_003_G01        Sst
    ## F2S4_160802_003_H01      L5 IT
    ## F2S4_160802_004_A01        Sst
    ## F2S4_160802_004_B01        Sst
    ## F2S4_160802_004_C01      L5 IT
    ## F2S4_160802_004_D01        Sst
    ## F2S4_160802_004_E01      Pvalb
    ## F2S4_160802_004_F01        Sst
    ## F2S4_160802_004_G01        Sst
    ## F2S4_160802_004_H01      Pvalb
    ## F2S4_160802_005_A01        Sst
    ## F2S4_160802_005_B01      Lamp5
    ## F2S4_160802_005_C01        Sst
    ## F2S4_160802_005_D01        Sst
    ## F2S4_160802_005_F01        Sst
    ## F2S4_160802_005_G01        Vip
    ## F2S4_160802_005_H01        Sst
    ## F2S4_160802_006_A01      Lamp5
    ## F2S4_160802_006_C01      Lamp5
    ## F2S4_160802_006_D01      Lamp5
    ## F2S4_160802_006_E01      Lamp5
    ## F2S4_160802_006_F01      Lamp5
    ## F2S4_160802_006_G01      Lamp5
    ## F2S4_160802_006_H01        Sst
    ## F2S4_160802_007_A01      Lamp5
    ## F2S4_160802_007_B01        Vip
    ## F2S4_160802_007_C01      Lamp5
    ## F2S4_160802_007_D01      Lamp5
    ## F2S4_160802_007_E01      Lamp5
    ## F2S4_160802_007_F01      Lamp5
    ## F2S4_160802_007_G01      Lamp5
    ## F2S4_160802_007_H01        Vip
    ## F2S4_160802_008_A01        Sst
    ## F2S4_160802_008_B01        Sst
    ## F2S4_160802_008_C01        Sst
    ## F2S4_160802_008_D01      Pvalb
    ## F2S4_160802_008_E01        Sst
    ## F2S4_160802_008_F01      L6 IT
    ## F2S4_160802_008_G01        Sst
    ## F2S4_160802_008_H01        Sst
    ## F2S4_160802_009_A01        Sst
    ## F2S4_160802_009_B01        Sst
    ## F2S4_160802_009_C01        Sst
    ## F2S4_160802_009_D01        Sst
    ## F2S4_160802_009_E01        Sst
    ## F2S4_160802_009_F01        Sst
    ## F2S4_160802_009_G01        Sst
    ## F2S4_160802_009_H01        Sst
    ## F2S4_160803_038_A01        Sst
    ## F2S4_160803_038_B01        Sst
    ## F2S4_160803_038_C01        Sst
    ## F2S4_160803_038_D01        Sst
    ## F2S4_160803_038_E01        Sst
    ## F2S4_160803_038_F01        Sst
    ## F2S4_160803_038_G01        Sst
    ## F2S4_160803_038_H01        Sst
    ## F2S4_160803_039_A01        Sst
    ## F2S4_160803_039_B01        Sst
    ## F2S4_160803_039_C01        Sst
    ## F2S4_160803_039_D01        Sst
    ## F2S4_160803_039_E01        Sst
    ## F2S4_160803_039_F01      Pvalb
    ## F2S4_160803_039_G01        Sst
    ## F2S4_160803_039_H01        Sst
    ## F2S4_160803_040_A01        Sst
    ## F2S4_160803_040_B01        Sst
    ## F2S4_160804_001_A01      Pvalb
    ## F2S4_160804_001_B01      Pvalb
    ## F2S4_160804_001_C01      Pvalb
    ## F2S4_160804_001_D01      Pvalb
    ## F2S4_160804_001_E01      Pvalb
    ## F2S4_160804_001_F01      Pvalb
    ## F2S4_160804_001_G01      Pvalb
    ## F2S4_160804_001_H01      Pvalb
    ## F2S4_160804_002_A01      Pvalb
    ## F2S4_160804_002_B01      Pvalb
    ## F2S4_160804_002_C01      Pvalb
    ## F2S4_160804_002_D01      Pvalb
    ## F2S4_160804_002_E01      Pvalb
    ## F2S4_160804_002_F01      Pvalb
    ## F2S4_160804_002_G01        Sst
    ## F2S4_160804_002_H01      Pvalb
    ## F2S4_160804_003_A01      Pvalb
    ## F2S4_160804_003_B01      Pvalb
    ## F2S4_160804_003_C01      Pvalb
    ## F2S4_160804_003_D01      Pvalb
    ## F2S4_160804_003_E01      Pvalb
    ## F2S4_160804_003_F01      Pvalb
    ## F2S4_160804_003_G01      Pvalb
    ## F2S4_160804_004_A01      Pvalb
    ## F2S4_160804_004_B01      Pvalb
    ## F2S4_160804_004_C01      Pvalb
    ## F2S4_160804_004_D01      Pvalb
    ## F2S4_160804_004_E01      Pvalb
    ## F2S4_160804_004_F01      Pvalb
    ## F2S4_160804_004_G01      Pvalb
    ## F2S4_160804_004_H01      Pvalb
    ## F2S4_160805_001_A01      Pvalb
    ## F2S4_160805_001_B01      Pvalb
    ## F2S4_160805_001_C01      Pvalb
    ## F2S4_160805_001_D01      Pvalb
    ## F2S4_160805_001_E01      Pvalb
    ## F2S4_160805_001_F01      Pvalb
    ## F2S4_160805_001_G01      Pvalb
    ## F2S4_160805_001_H01      Pvalb
    ## F2S4_160805_002_A01      Pvalb
    ## F2S4_160805_002_B01      Pvalb
    ## F2S4_160805_002_C01      Pvalb
    ## F2S4_160805_002_D01      Pvalb
    ## F2S4_160805_002_E01      Pvalb
    ## F2S4_160805_002_F01      Pvalb
    ## F2S4_160805_002_G01      Pvalb
    ## F2S4_160805_002_H01      Pvalb
    ## F2S4_160805_003_A01      Pvalb
    ## F2S4_160805_003_B01      Pvalb
    ## F2S4_160805_003_C01      Pvalb
    ## F2S4_160805_003_D01      Pvalb
    ## F2S4_160805_003_E01      Pvalb
    ## F2S4_160805_003_F01      Pvalb
    ## F2S4_160805_003_G01      Pvalb
    ## F2S4_160805_003_H01      Pvalb
    ## F2S4_160805_004_A01      Pvalb
    ## F2S4_160805_004_B01      Pvalb
    ## F2S4_160805_004_C01      Pvalb
    ## F2S4_160805_004_D01      Pvalb
    ## F2S4_160805_004_E01      Pvalb
    ## F2S4_160805_004_F01      Pvalb
    ## F2S4_160805_004_G01      Pvalb
    ## F2S4_160805_004_H01      Pvalb
    ## F2S4_160805_005_A01      Pvalb
    ## F2S4_160805_005_B01      Pvalb
    ## F2S4_160805_005_C01      Pvalb
    ## F2S4_160805_005_D01      Pvalb
    ## F2S4_160805_005_E01      Pvalb
    ## F2S4_160805_005_F01      Pvalb
    ## F2S4_160805_005_G01      Pvalb
    ## F2S4_160805_005_H01      Pvalb
    ## F2S4_160805_006_A01      Pvalb
    ## F2S4_160805_006_B01      Pvalb
    ## F2S4_160805_006_C01      Pvalb
    ## F2S4_160805_006_D01      Pvalb
    ## F2S4_160805_006_E01      Pvalb
    ## F2S4_160805_006_F01      Pvalb
    ## F2S4_160805_006_G01      Pvalb
    ## F2S4_160805_006_H01      Pvalb
    ## F2S4_160805_007_G01      Pvalb
    ## F2S4_160805_007_H01      Pvalb
    ## F2S4_160805_008_B01      Pvalb
    ## F2S4_160805_008_C01      Pvalb
    ## F2S4_160805_008_D01      Pvalb
    ## F2S4_160805_008_E01      Pvalb
    ## F2S4_160805_008_F01      Pvalb
    ## F2S4_160805_008_G01      Pvalb
    ## F2S4_160805_008_H01      Pvalb
    ## F2S4_160805_009_A01      Pvalb
    ## F2S4_160805_009_B01      Pvalb
    ## F2S4_160805_009_C01      Pvalb
    ## F2S4_160805_009_D01      Pvalb
    ## F2S4_160805_009_E01      Pvalb
    ## F2S4_160805_009_F01      Pvalb
    ## F2S4_160805_009_G01      Pvalb
    ## F2S4_160805_009_H01      Pvalb
    ## F2S4_160805_010_G01      Pvalb
    ## F2S4_160805_010_H01      Pvalb
    ## F2S4_160808_032_B01        Vip
    ## F2S4_160808_032_C01        Vip
    ## F2S4_160808_032_E01        Vip
    ## F2S4_160808_032_F01        Vip
    ## F2S4_160808_032_G01        Vip
    ## F2S4_160808_032_H01        Vip
    ## F2S4_160811_001_A01      Pvalb
    ## F2S4_160811_001_B01      Pvalb
    ## F2S4_160811_001_C01      Pvalb
    ## F2S4_160811_001_D01      Pvalb
    ## F2S4_160811_001_E01      Pvalb
    ## F2S4_160811_001_F01        Sst
    ## F2S4_160811_001_G01      Pvalb
    ## F2S4_160811_001_H01      Pvalb
    ## F2S4_160811_002_A01      Pvalb
    ## F2S4_160811_002_B01      Pvalb
    ## F2S4_160811_002_C01      Pvalb
    ## F2S4_160811_002_D01      Pvalb
    ## F2S4_160811_002_E01      Pvalb
    ## F2S4_160811_002_F01      Pvalb
    ## F2S4_160811_002_G01      Pvalb
    ## F2S4_160811_002_H01      Pvalb
    ## F2S4_160811_003_A01      Pvalb
    ## F2S4_160811_003_B01      Pvalb
    ## F2S4_160811_003_C01      Pvalb
    ## F2S4_160811_003_D01      Pvalb
    ## F2S4_160811_003_E01      Pvalb
    ## F2S4_160811_003_F01      Pvalb
    ## F2S4_160811_003_G01      Pvalb
    ## F2S4_160811_003_H01      Pvalb
    ## F2S4_160811_004_F01      Pvalb
    ## F2S4_160811_004_G01        Sst
    ## F2S4_160811_004_H01        Sst
    ## F2S4_160812_001_A01        Sst
    ## F2S4_160812_001_B01      Pvalb
    ## F2S4_160812_001_C01        Sst
    ## F2S4_160812_001_D01      Pvalb
    ## F2S4_160812_001_E01      Pvalb
    ## F2S4_160812_001_F01      Pvalb
    ## F2S4_160812_001_G01        Sst
    ## F2S4_160812_001_H01      Pvalb
    ## F2S4_160812_002_A01      Pvalb
    ## F2S4_160812_002_B01      Pvalb
    ## F2S4_160812_002_C01      Pvalb
    ## F2S4_160812_002_D01      Pvalb
    ## F2S4_160812_002_E01        Sst
    ## F2S4_160812_002_F01      Pvalb
    ## F2S4_160812_002_G01      Pvalb
    ## F2S4_160812_002_H01      Pvalb
    ## F2S4_160812_003_A01      Pvalb
    ## F2S4_160812_003_B01      Pvalb
    ## F2S4_160812_003_C01        Sst
    ## F2S4_160812_003_E01        Sst
    ## F2S4_160812_003_F01      Pvalb
    ## F2S4_160812_003_G01      Pvalb
    ## F2S4_160812_003_H01      Pvalb
    ## F2S4_160818_001_A01      L5 PT
    ## F2S4_160818_001_B01      L5 PT
    ## F2S4_160818_001_C01      L5 PT
    ## F2S4_160818_001_D01      L5 IT
    ## F2S4_160818_001_E01      L5 IT
    ## F2S4_160818_001_F01      L5 IT
    ## F2S4_160818_001_G01      L5 IT
    ## F2S4_160818_001_H01      L5 IT
    ## F2S4_160818_002_A01      L5 IT
    ## F2S4_160818_002_B01      L5 IT
    ## F2S4_160818_002_C01      L5 PT
    ## F2S4_160818_002_D01      L6 IT
    ## F2S4_160818_002_E01      L5 PT
    ## F2S4_160818_002_F01      L5 IT
    ## F2S4_160818_002_G01      L5 PT
    ## F2S4_160818_002_H01      L5 IT
    ## F2S4_160818_003_A01      L5 IT
    ## F2S4_160818_003_B01      L5 IT
    ## F2S4_160818_003_C01      L5 PT
    ## F2S4_160818_003_D01      L5 IT
    ## F2S4_160818_003_E01      L5 IT
    ## F2S4_160818_003_F01      L5 PT
    ## F2S4_160818_003_G01      L5 IT
    ## F2S4_160818_003_H01      L5 PT
    ## F2S4_160818_004_A01      L5 IT
    ## F2S4_160818_004_B01      L5 IT
    ## F2S4_160818_004_C01      L5 PT
    ## F2S4_160818_004_D01      L5 IT
    ## F2S4_160818_004_E01      L5 IT
    ## F2S4_160818_004_F01      L6 IT
    ## F2S4_160818_004_G01      L5 IT
    ## F2S4_160818_004_H01      L5 IT
    ## F2S4_160818_005_A01      L5 PT
    ## F2S4_160818_005_B01        Sst
    ## F2S4_160818_005_C01      L5 IT
    ## F2S4_160818_005_D01      L5 IT
    ## F2S4_160818_005_E01      L5 IT
    ## F2S4_160818_005_F01      L6 IT
    ## F2S4_160818_005_G01      L5 IT
    ## F2S4_160818_005_H01      L5 IT
    ## F2S4_160818_006_A01      L5 IT
    ## F2S4_160818_006_B01      L5 IT
    ## F2S4_160818_006_C01      L6 IT
    ## F2S4_160818_006_D01      L5 PT
    ## F2S4_160818_006_E01      L5 PT
    ## F2S4_160818_006_F01      L5 PT
    ## F2S4_160818_006_G01      L5 PT
    ## F2S4_160818_006_H01      L5 IT
    ## F2S4_160818_007_A01      L5 IT
    ## F2S4_160818_007_B01      L5 IT
    ## F2S4_160818_007_C01      L6 CT
    ## F2S4_160818_007_D01      L6 IT
    ## F2S4_160818_007_E01      L5 IT
    ## F2S4_160818_007_F01      L5 PT
    ## F2S4_160818_007_G01      L5 IT
    ## F2S4_160818_007_H01      L5 PT
    ## F2S4_160818_008_A01      L5 IT
    ## F2S4_160818_008_B01      L5 IT
    ## F2S4_160818_008_C01      L5 IT
    ## F2S4_160818_008_D01      L5 IT
    ## F2S4_160818_008_E01        Sst
    ## F2S4_160818_008_F01      L5 IT
    ## F2S4_160818_008_G01      L5 IT
    ## F2S4_160818_008_H01      L5 IT
    ## F2S4_160818_009_A01      L5 IT
    ## F2S4_160818_009_B01      L5 PT
    ## F2S4_160818_009_C01        Sst
    ## F2S4_160818_009_D01      L5 IT
    ## F2S4_160818_009_E01      L6 CT
    ## F2S4_160818_009_F01      L6 CT
    ## F2S4_160818_009_H01        Sst
    ## F2S4_160818_010_A01        Sst
    ## F2S4_160824_049_A01      Pvalb
    ## F2S4_160824_049_B01      Pvalb
    ## F2S4_160824_049_C01      Pvalb
    ## F2S4_160824_049_D01      Pvalb
    ## F2S4_160824_049_E01      Pvalb
    ## F2S4_160824_049_F01      Pvalb
    ## F2S4_160824_049_G01      Pvalb
    ## F2S4_160824_049_H01      Pvalb
    ## F2S4_160824_050_A01      Pvalb
    ## F2S4_160824_050_B01      Pvalb
    ## F2S4_160824_050_C01      Pvalb
    ## F2S4_160824_050_D01      Pvalb
    ## F2S4_160824_050_E01      Pvalb
    ## F2S4_160824_050_F01      Pvalb
    ## F2S4_160824_050_G01      Pvalb
    ## F2S4_160824_050_H01      Pvalb
    ## F2S4_160824_051_A01      Pvalb
    ## F2S4_160824_051_B01      Pvalb
    ## F2S4_160824_051_C01      Pvalb
    ## F2S4_160824_051_D01      Pvalb
    ## F2S4_160824_051_E01      Pvalb
    ## F2S4_160824_051_F01      Pvalb
    ## F2S4_160824_051_G01      Pvalb
    ## F2S4_160824_051_H01      Pvalb
    ## F2S4_160824_052_A01      Pvalb
    ## F2S4_160824_052_B01      Pvalb
    ## F2S4_160824_052_C01      Pvalb
    ## F2S4_160824_052_D01      Pvalb
    ## F2S4_160824_052_E01      Pvalb
    ## F2S4_160824_052_F01      Pvalb
    ## F2S4_160824_052_G01      Pvalb
    ## F2S4_160824_052_H01      Pvalb
    ## F2S4_160824_053_A01      Pvalb
    ## F2S4_160824_053_B01      Pvalb
    ## F2S4_160824_053_C01      Pvalb
    ## F2S4_160824_053_D01      Pvalb
    ## F2S4_160824_053_E01      Pvalb
    ## F2S4_160824_053_F01      Pvalb
    ## F2S4_160824_053_G01      Pvalb
    ## F2S4_160824_053_H01      Pvalb
    ## F2S4_160824_054_A01      Pvalb
    ## F2S4_160824_054_B01      Pvalb
    ## F2S4_160824_054_C01 Macrophage
    ## F2S4_160824_054_D01      Pvalb
    ## F2S4_160824_054_E01      Pvalb
    ## F2S4_160824_054_F01      Pvalb
    ## F2S4_160825_055_A01      Pvalb
    ## F2S4_160825_055_B01      Pvalb
    ## F2S4_160825_055_C01      Pvalb
    ## F2S4_160825_055_D01      Pvalb
    ## F2S4_160825_055_E01      Pvalb
    ## F2S4_160825_055_F01      Pvalb
    ## F2S4_160825_055_G01      Pvalb
    ## F2S4_160825_055_H01      Pvalb
    ## F2S4_160825_056_A01      Pvalb
    ## F2S4_160825_056_B01      Pvalb
    ## F2S4_160825_056_C01      Pvalb
    ## F2S4_160825_056_D01      Pvalb
    ## F2S4_160825_056_E01      Pvalb
    ## F2S4_160825_056_F01      Pvalb
    ## F2S4_160825_056_G01      Pvalb
    ## F2S4_160825_056_H01      Pvalb
    ## F2S4_160825_057_A01        Sst
    ## F2S4_160825_057_B01      Pvalb
    ## F2S4_160825_057_C01      Pvalb
    ## F2S4_160825_057_D01      Pvalb
    ## F2S4_160825_057_E01      Pvalb
    ## F2S4_160825_057_F01      Pvalb
    ## F2S4_160825_057_G01      Pvalb
    ## F2S4_160825_057_H01      Pvalb
    ## F2S4_160825_058_A01      Pvalb
    ## F2S4_160825_058_B01      Pvalb
    ## F2S4_160825_058_C01      Pvalb
    ## F2S4_160825_058_D01      Pvalb
    ## F2S4_160825_058_E01      Pvalb
    ## F2S4_160825_058_F01      Pvalb
    ## F2S4_160825_058_G01      Pvalb
    ## F2S4_160825_058_H01      Pvalb
    ## F2S4_160825_059_A01      Pvalb
    ## F2S4_160825_059_B01      Pvalb
    ## F2S4_160825_059_C01      Pvalb
    ## F2S4_160825_059_D01      Pvalb
    ## F2S4_160825_059_E01      Pvalb
    ## F2S4_160825_059_F01      Pvalb
    ## F2S4_160825_059_G01      Pvalb
    ## F2S4_160825_059_H01      Pvalb
    ## F2S4_160825_060_A01      Pvalb
    ## F2S4_160825_060_B01      Pvalb
    ## F2S4_160825_060_C01      Pvalb
    ## F2S4_160825_060_D01      Pvalb
    ## F2S4_160825_060_E01      Pvalb
    ## F2S4_160825_060_F01      Pvalb
    ## F2S4_160825_060_G01      Pvalb
    ## F2S4_160825_060_H01      Pvalb
    ## F2S4_160825_061_A01      Pvalb
    ## F2S4_160825_061_B01      Pvalb
    ## F2S4_160825_061_C01      Pvalb
    ## F2S4_160825_061_D01      Pvalb
    ## F2S4_160825_061_E01      Pvalb
    ## F2S4_160825_061_F01      Pvalb
    ## F2S4_160825_061_G01      Pvalb
    ## F2S4_160825_061_H01      Pvalb
    ## F2S4_160826_001_A01      Pvalb
    ## F2S4_160826_001_B01      Pvalb
    ## F2S4_160826_001_C01      Lamp5
    ## F2S4_160826_001_D01        Sst
    ## F2S4_160826_001_E01      Pvalb
    ## F2S4_160826_001_F01      Pvalb
    ## F2S4_160826_001_G01      Pvalb
    ## F2S4_160826_001_H01      Pvalb
    ## F2S4_160826_002_A01      Pvalb
    ## F2S4_160826_002_B01      Pvalb
    ## F2S4_160826_002_C01      Pvalb
    ## F2S4_160826_002_D01      Pvalb
    ## F2S4_160826_002_E01      Pvalb
    ## F2S4_160826_002_F01      Pvalb
    ## F2S4_160826_002_G01      Pvalb
    ## F2S4_160826_002_H01      Pvalb
    ## F2S4_160826_003_A01        Sst
    ## F2S4_160826_003_B01      Pvalb
    ## F2S4_160826_003_C01      Pvalb
    ## F2S4_160826_003_D01        Sst
    ## F2S4_160826_003_E01      Pvalb
    ## F2S4_160826_003_F01      Pvalb
    ## F2S4_160826_003_G01      Pvalb
    ## F2S4_160826_003_H01      Pvalb
    ## F2S4_160826_004_A01      Pvalb
    ## F2S4_160826_004_B01      Pvalb
    ## F2S4_160826_004_C01      Pvalb
    ## F2S4_160826_004_D01      Pvalb
    ## F2S4_160826_004_E01      Pvalb
    ## F2S4_160826_004_F01      Pvalb
    ## F2S4_160826_004_G01      Pvalb
    ## F2S4_160826_004_H01      Pvalb
    ## F2S4_160829_001_A01      Pvalb
    ## F2S4_160829_001_B01      Pvalb
    ## F2S4_160829_001_C01      Pvalb
    ## F2S4_160829_001_D01      Pvalb
    ## F2S4_160829_001_E01      Pvalb
    ## F2S4_160829_001_F01      Pvalb
    ## F2S4_160829_001_G01      Pvalb
    ## F2S4_160829_001_H01      Pvalb
    ## F2S4_160829_002_A01      Pvalb
    ## F2S4_160829_002_B01      Pvalb
    ## F2S4_160829_002_C01      Pvalb
    ## F2S4_160829_002_D01      Pvalb
    ## F2S4_160829_002_E01      Pvalb
    ## F2S4_160829_002_F01      Pvalb
    ## F2S4_160829_002_G01      Pvalb
    ## F2S4_160829_002_H01      Pvalb
    ## F2S4_160829_003_A01      Pvalb
    ## F2S4_160829_003_B01      Pvalb
    ## F2S4_160829_003_C01      Pvalb
    ## F2S4_160829_003_D01      Pvalb
    ## F2S4_160829_003_E01      Pvalb
    ## F2S4_160829_003_F01      Pvalb
    ## F2S4_160829_003_G01      Pvalb
    ## F2S4_160829_003_H01      Pvalb
    ## F2S4_160829_004_A01      Pvalb
    ## F2S4_160829_004_B01      Pvalb
    ## F2S4_160829_004_C01      Pvalb
    ## F2S4_160829_004_D01      Pvalb
    ## F2S4_160829_004_E01      Pvalb
    ## F2S4_160829_004_F01      Pvalb
    ## F2S4_160829_004_G01      Pvalb
    ## F2S4_160829_004_H01      Pvalb
    ## F2S4_160829_005_B01      Pvalb
    ## F2S4_160829_005_C01      Pvalb
    ## F2S4_160829_005_D01      Pvalb
    ## F2S4_160829_005_E01      Pvalb
    ## F2S4_160829_005_F01      Pvalb
    ## F2S4_160829_005_G01      Pvalb
    ## F2S4_160829_005_H01      Pvalb
    ## F2S4_160829_006_B01      Pvalb
    ## F2S4_160829_006_C01      Pvalb
    ## F2S4_160829_006_D01      Pvalb
    ## F2S4_160829_006_E01      Pvalb
    ## F2S4_160829_006_F01      Pvalb
    ## F2S4_160829_006_G01      Pvalb
    ## F2S4_160829_006_H01      Pvalb
    ## F2S4_160901_001_A01      Astro
    ## F2S4_160901_001_B01      Astro
    ## F2S4_160901_001_C01      Pvalb
    ## F2S4_160901_001_D01        SMC
    ## F2S4_160901_001_E01      Astro
    ## F2S4_160901_001_F01       VLMC
    ## F2S4_160901_001_G01        SMC
    ## F2S4_160901_001_H01      Astro
    ## F2S4_160901_002_C01      L5 IT
    ## F2S4_160901_002_D01        SMC
    ## F2S4_160901_002_E01         L4
    ## F2S4_160901_002_F01      Astro
    ## F2S4_160901_002_G01        Vip
    ## F2S4_160901_002_H01      L5 IT
    ## F2S4_160901_003_A01      Astro
    ## F2S4_160901_003_B01         L4
    ## F2S4_160901_003_D01      Astro
    ## F2S4_160901_003_E01      Astro
    ## F2S4_160901_003_F01      Astro
    ## F2S4_160901_003_G01        SMC
    ## F2S4_160901_004_B01      Astro
    ## F2S4_160901_004_C01      Astro
    ## F2S4_160901_004_D01         L4
    ## F2S4_160901_004_E01      Astro
    ## F2S4_160901_004_F01      L6 IT
    ## F2S4_160901_004_G01      Astro
    ## F2S4_160901_004_H01      L6 IT
    ## F2S4_160901_005_A01        SMC
    ## F2S4_160901_005_B01      Astro
    ## F2S4_160901_005_C01      Astro
    ## F2S4_160901_005_D01         NP
    ## F2S4_160901_005_E01      Astro
    ## F2S4_160901_005_G01      L6 IT
    ## F2S4_160901_005_H01      Astro
    ## F2S4_160901_006_A01      Astro
    ## F2S4_160901_006_B01      Astro
    ## F2S4_160901_006_C01        Sst
    ## F2S4_160901_006_D01      L5 PT
    ## F2S4_160901_006_E01         NP
    ## F2S4_160901_006_F01    L2/3 IT
    ## F2S4_160901_006_G01      Astro
    ## F2S4_160901_006_H01        SMC
    ## F2S4_160901_007_A01      Astro
    ## F2S4_160901_007_C01      Astro
    ## F2S4_160901_007_D01         L4
    ## F2S4_160901_007_E01         L4
    ## F2S4_160901_007_F01      Lamp5
    ## F2S4_160901_007_H01      Astro
    ## F2S4_160901_008_A01      Pvalb
    ## F2S4_160901_008_B01        SMC
    ## F2S4_160901_008_C01      Astro
    ## F2S4_160901_008_D01        Vip
    ## F2S4_160901_008_E01      Astro
    ## F2S4_160901_008_F01         L4
    ## F2S4_160901_008_G01      Astro
    ## F2S4_160901_008_H01      Astro
    ## F2S4_160901_009_A01        SMC
    ## F2S4_160901_009_B01        SMC
    ## F2S4_160901_009_C01      Meis2
    ## F2S4_160901_009_D01        SMC
    ## F2S4_160901_009_E01        SMC
    ## F2S4_160901_009_F01        SMC
    ## F2S4_160901_009_G01      L5 IT
    ## F2S4_160901_009_H01        Sst
    ## F2S4_160901_010_A01      L5 PT
    ## F2S4_160901_010_B01         L4
    ## F2S4_160901_010_C01         L4
    ## F2S4_160901_010_D01        SMC
    ## F2S4_160901_010_E01        SMC
    ## F2S4_160901_010_F01        SMC
    ## F2S4_160901_010_G01      Astro
    ## F2S4_160901_010_H01        SMC
    ## F2S4_160901_011_A01        Vip
    ## F2S4_160901_011_B01        SMC
    ## F2S4_160901_011_D01        SMC
    ## F2S4_160901_011_E01        SMC
    ## F2S4_160901_011_F01        SMC
    ## F2S4_160901_012_A01        SMC
    ## F2S4_160901_012_C01       VLMC
    ## F2S4_160901_012_D01        Vip
    ## F2S4_160901_012_E01    L2/3 IT
    ## F2S4_160901_012_F01        SMC
    ## F2S4_160901_012_G01        SMC
    ## F2S4_160901_012_H01      Pvalb
    ## F2S4_160901_013_A01      Lamp5
    ## F2S4_160901_013_B01        SMC
    ## F2S4_160901_013_C01      Astro
    ## F2S4_160901_013_D01         L4
    ## F2S4_160901_013_G01        SMC
    ## F2S4_160901_013_H01        Vip
    ## F2S4_160901_014_A01      Astro
    ## F2S4_160901_014_B01        SMC
    ## F2S4_160901_014_C01        SMC
    ## F2S4_160901_014_D01      Pvalb
    ## F2S4_160901_014_F01        SMC
    ## F2S4_160901_014_G01      Astro
    ## F2S4_160901_014_H01      Astro
    ## F2S4_160901_015_A01      Astro
    ## F2S4_160901_015_B01        Vip
    ## F2S4_160901_015_C01        Vip
    ## F2S4_160901_015_D01         L4
    ## F2S4_160901_015_E01         NP
    ## F2S4_160901_015_F01        SMC
    ## F2S4_160901_015_G01      Pvalb
    ## F2S4_160901_015_H01        Sst
    ## F2S4_160901_016_A01      Astro
    ## F2S4_160901_016_B01      L6 CT
    ## F2S4_160901_016_C01        SMC
    ## F2S4_160901_016_D01        SMC
    ## F2S4_160901_016_E01        SMC
    ## F2S4_160901_016_G01        SMC
    ## F2S4_160901_016_H01      Pvalb
    ## F2S4_161003_001_A01         NP
    ## F2S4_161003_001_B01         NP
    ## F2S4_161003_001_C01         NP
    ## F2S4_161003_001_D01      L5 PT
    ## F2S4_161003_001_E01         L4
    ## F2S4_161003_001_F01         L4
    ## F2S4_161003_001_G01         NP
    ## F2S4_161003_002_A01      L5 PT
    ## F2S4_161003_002_B01         NP
    ## F2S4_161003_002_C01      L5 PT
    ## F2S4_161003_002_D01      L5 PT
    ## F2S4_161003_002_E01      L5 PT
    ## F2S4_161003_002_G01      L5 PT
    ## F2S4_161003_002_H01      L6 IT
    ## F2S4_161003_003_A01      L5 PT
    ## F2S4_161003_003_B01         NP
    ## F2S4_161003_003_C01      L5 PT
    ## F2S4_161003_003_E01      L5 PT
    ## F2S4_161003_003_F01         NP
    ## F2S4_161003_003_G01         NP
    ## F2S4_161004_001_A01         L4
    ## F2S4_161004_001_B01         L4
    ## F2S4_161004_001_C01         L4
    ## F2S4_161004_001_D01         L4
    ## F2S4_161004_001_E01         L4
    ## F2S4_161004_001_F01         L4
    ## F2S4_161004_001_G01         L4
    ## F2S4_161004_001_H01         L4
    ## F2S4_161004_002_A01         L4
    ## F2S4_161004_002_B01         L4
    ## F2S4_161004_002_C01         L4
    ## F2S4_161004_002_D01         L4
    ## F2S4_161004_002_E01         L4
    ## F2S4_161004_002_F01         L4
    ## F2S4_161004_002_G01         L4
    ## F2S4_161004_002_H01         L4
    ## F2S4_161004_003_A01         L4
    ## F2S4_161004_003_B01         L4
    ## F2S4_161004_003_C01         L4
    ## F2S4_161004_003_D01         L4
    ## F2S4_161004_003_E01         L4
    ## F2S4_161004_003_F01         L4
    ## F2S4_161004_003_G01      L5 PT
    ## F2S4_161004_003_H01      L5 IT
    ## F2S4_161004_004_A01         L4
    ## F2S4_161004_004_B01         L4
    ## F2S4_161004_004_C01         L4
    ## F2S4_161004_004_D01         L4
    ## F2S4_161004_004_E01      L5 IT
    ## F2S4_161004_004_F01      L5 IT
    ## F2S4_161004_004_G01         L4
    ## F2S4_161004_004_H01         L4
    ## F2S4_161004_005_A01         L4
    ## F2S4_161004_005_B01         L4
    ## F2S4_161004_005_C01         L4
    ## F2S4_161004_005_D01         L4
    ## F2S4_161004_005_E01         L4
    ## F2S4_161004_005_F01         L4
    ## F2S4_161004_005_G01      L5 IT
    ## F2S4_161004_005_H01         L4
    ## F2S4_161004_008_A01      L5 IT
    ## F2S4_161004_008_B01         L4
    ## F2S4_161004_008_C01      L5 IT
    ## F2S4_161004_008_D01         L4
    ## F2S4_161004_008_E01         L4
    ## F2S4_161004_008_F01         L4
    ## F2S4_161004_008_G01         L4
    ## F2S4_161004_008_H01         L4
    ## F2S4_161004_009_B01         L4
    ## F2S4_161004_009_C01         L4
    ## F2S4_161004_009_D01         L4
    ## F2S4_161004_009_E01         L4
    ## F2S4_161004_009_F01         L4
    ## F2S4_161004_009_G01         L4
    ## F2S4_161004_009_H01         L4
    ## F2S4_161004_010_A01    L2/3 IT
    ## F2S4_161004_010_B01      Lamp5
    ## F2S4_161004_010_C01    L2/3 IT
    ## F2S4_161004_010_D01        Vip
    ## F2S4_161004_010_E01    L2/3 IT
    ## F2S4_161004_010_F01    L2/3 IT
    ## F2S4_161004_010_G01    L2/3 IT
    ## F2S4_161004_010_H01    L2/3 IT
    ## F2S4_161004_011_A01    L2/3 IT
    ## F2S4_161004_011_B01    L2/3 IT
    ## F2S4_161004_011_C01    L2/3 IT
    ## F2S4_161004_011_D01        Vip
    ## F2S4_161004_011_E01      Lamp5
    ## F2S4_161004_011_F01         L4
    ## F2S4_161004_011_G01        Vip
    ## F2S4_161004_011_H01    L2/3 IT
    ## F2S4_161004_012_A01    L2/3 IT
    ## F2S4_161004_012_B01    L2/3 IT
    ## F2S4_161004_012_C01    L2/3 IT
    ## F2S4_161004_012_D01    L2/3 IT
    ## F2S4_161004_012_E01    L2/3 IT
    ## F2S4_161004_012_F01    L2/3 IT
    ## F2S4_161004_012_G01    L2/3 IT
    ## F2S4_161004_012_H01    L2/3 IT
    ## F2S4_161004_013_A01    L2/3 IT
    ## F2S4_161004_013_B01    L2/3 IT
    ## F2S4_161004_013_C01    L2/3 IT
    ## F2S4_161004_013_D01      Lamp5
    ## F2S4_161004_013_E01    L2/3 IT
    ## F2S4_161004_013_F01    L2/3 IT
    ## F2S4_161004_013_G01         L4
    ## F2S4_161004_013_H01    L2/3 IT
    ## F2S4_161004_014_A01    L2/3 IT
    ## F2S4_161004_014_B01    L2/3 IT
    ## F2S4_161004_014_C01    L2/3 IT
    ## F2S4_161004_014_D01    L2/3 IT
    ## F2S4_161004_014_E01    L2/3 IT
    ## F2S4_161004_014_F01    L2/3 IT
    ## F2S4_161004_014_G01    L2/3 IT
    ## F2S4_161004_014_H01    L2/3 IT
    ## F2S4_161004_015_A01    L2/3 IT
    ## F2S4_161004_015_B01    L2/3 IT
    ## F2S4_161004_015_C01    L2/3 IT
    ## F2S4_161004_015_D01    L2/3 IT
    ## F2S4_161004_015_E01    L2/3 IT
    ## F2S4_161004_015_F01    L2/3 IT
    ## F2S4_161004_015_G01    L2/3 IT
    ## F2S4_161004_015_H01    L2/3 IT
    ## F2S4_161004_016_A01    L2/3 IT
    ## F2S4_161004_016_B01    L2/3 IT
    ## F2S4_161004_016_C01    L2/3 IT
    ## F2S4_161004_016_E01    L2/3 IT
    ## F2S4_161004_016_F01    L2/3 IT
    ## F2S4_161004_016_G01      Lamp5
    ## F2S4_161004_016_H01    L2/3 IT
    ## F2S4_161004_017_A01    L2/3 IT
    ## F2S4_161004_017_B01      Lamp5
    ## F2S4_161004_017_C01    L2/3 IT
    ## F2S4_161004_017_D01    L2/3 IT
    ## F2S4_161004_017_E01    L2/3 IT
    ## F2S4_161004_017_F01    L2/3 IT
    ## F2S4_161004_017_G01    L2/3 IT
    ## F2S4_161004_017_H01    L2/3 IT
    ## F2S4_161004_018_A01    L2/3 IT
    ## F2S4_161004_018_B01    L2/3 IT
    ## F2S4_161004_018_C01    L2/3 IT
    ## F2S4_161004_018_D01    L2/3 IT
    ## F2S4_161004_018_E01         L4
    ## F2S4_161004_018_F01        Vip
    ## F2S4_161004_018_G01    L2/3 IT
    ## F2S4_161004_018_H01    L2/3 IT
    ## F2S4_161004_019_A01        L6b
    ## F2S4_161004_019_B01      L6 CT
    ## F2S4_161004_019_C01         L4
    ## F2S4_161004_019_D01        L6b
    ## F2S4_161004_019_E01        L6b
    ## F2S4_161004_019_F01        L6b
    ## F2S4_161004_019_G01      L6 CT
    ## F2S4_161004_019_H01        L6b
    ## F2S4_161004_020_A01      L6 CT
    ## F2S4_161004_020_B01      L5 IT
    ## F2S4_161004_020_E01        Vip
    ## F2S4_161004_020_F01      L6 IT
    ## F2S4_161004_020_G01         L4
    ## F2S4_161004_020_H01       Sncg
    ## F2S4_161004_021_A01         L4
    ## F2S4_161004_021_B01        Vip
    ## F2S4_161004_021_C01         L4
    ## F2S4_161004_021_D01         L4
    ## F2S4_161004_021_F01         L4
    ## F2S4_161004_021_G01         L4
    ## F2S4_161004_021_H01         L4
    ## F2S4_161004_022_A01         L4
    ## F2S4_161004_022_B01         L4
    ## F2S4_161004_022_C01         L4
    ## F2S4_161004_022_D01         L4
    ## F2S4_161004_022_E01         L4
    ## F2S4_161004_022_F01         L4
    ## F2S4_161004_022_G01         L4
    ## F2S4_161004_023_A01         L4
    ## F2S4_161004_023_B01         L4
    ## F2S4_161004_023_C01         L4
    ## F2S4_161004_023_D01         L4
    ## F2S4_161004_023_E01         L4
    ## F2S4_161004_023_F01         L4
    ## F2S4_161004_023_G01         L4
    ## F2S4_161004_023_H01         L4
    ## F2S4_161007_001_A01         L4
    ## F2S4_161007_001_B01      L5 IT
    ## F2S4_161007_001_C01         L4
    ## F2S4_161007_001_D01         L4
    ## F2S4_161007_001_E01         L4
    ## F2S4_161007_001_F01         L4
    ## F2S4_161007_001_G01      L5 IT
    ## F2S4_161007_001_H01      L5 IT
    ## F2S4_161007_002_A01         L4
    ## F2S4_161007_002_B01      L5 IT
    ## F2S4_161007_002_C01         L4
    ## F2S4_161007_002_D01      L5 IT
    ## F2S4_161007_002_E01      L5 IT
    ## F2S4_161007_002_F01         L4
    ## F2S4_161007_002_G01         L4
    ## F2S4_161007_002_H01      L5 IT
    ## F2S4_161007_003_A01         L4
    ## F2S4_161007_003_C01         L4
    ## F2S4_161007_003_D01      L5 IT
    ## F2S4_161007_003_F01         L4
    ## F2S4_161007_003_G01         L4
    ## F2S4_161007_003_H01      L5 IT
    ## F2S4_161007_004_A01      L5 IT
    ## F2S4_161007_004_B01      L5 IT
    ## F2S4_161007_004_D01      L5 IT
    ## F2S4_161007_004_E01      L5 IT
    ## F2S4_161007_004_F01         L4
    ## F2S4_161007_004_G01         L4
    ## F2S4_161007_004_H01      L5 IT
    ## F2S4_161007_005_A01         L4
    ## F2S4_161007_005_B01         L4
    ## F2S4_161007_005_C01      L5 IT
    ## F2S4_161007_005_D01         L4
    ## F2S4_161007_005_G01         L4
    ## F2S4_161007_005_H01      L5 IT
    ## F2S4_161007_006_A01      L5 IT
    ## F2S4_161007_006_B01         L4
    ## F2S4_161007_006_C01      L5 IT
    ## F2S4_161007_006_D01         L4
    ## F2S4_161007_006_F01      L5 IT
    ## F2S4_161007_006_G01         L4
    ## F2S4_161007_006_H01         L4
    ## F2S4_161007_010_C01      L5 IT
    ## F2S4_161007_010_D01         L4
    ## F2S4_161007_010_E01      L5 IT
    ## F2S4_161007_010_F01      L5 IT
    ## F2S4_161007_010_G01         L4
    ## F2S4_161007_010_H01      L5 IT
    ## F2S4_161007_011_B01      L5 IT
    ## F2S4_161007_011_C01      L5 IT
    ## F2S4_161007_011_D01      L5 IT
    ## F2S4_161007_011_E01         L4
    ## F2S4_161007_011_F01      L5 IT
    ## F2S4_161007_011_G01         L4
    ## F2S4_161007_011_H01      L5 IT
    ## F2S4_161007_012_A01         L4
    ## F2S4_161007_012_B01      L5 IT
    ## F2S4_161007_012_C01         L4
    ## F2S4_161007_012_D01         L4
    ## F2S4_161007_012_E01      L5 IT
    ## F2S4_161007_012_F01      L5 IT
    ## F2S4_161007_012_G01      L5 IT
    ## F2S4_161007_012_H01      L5 IT
    ## F2S4_161010_001_A01        Sst
    ## F2S4_161010_001_B01        Sst
    ## F2S4_161010_001_C01        Sst
    ## F2S4_161010_001_D01        Sst
    ## F2S4_161010_001_E01        Sst
    ## F2S4_161010_001_F01        Sst
    ## F2S4_161010_001_G01        Sst
    ## F2S4_161010_001_H01        Sst
    ## F2S4_161010_002_A01        Sst
    ## F2S4_161010_002_B01        Sst
    ## F2S4_161010_002_C01        Sst
    ## F2S4_161010_002_D01        Sst
    ## F2S4_161010_002_E01      Pvalb
    ## F2S4_161010_002_F01        Sst
    ## F2S4_161010_002_G01        Sst
    ## F2S4_161010_002_H01        Sst
    ## F2S4_161010_003_A01        Sst
    ## F2S4_161010_003_B01        Sst
    ## F2S4_161010_003_C01      Pvalb
    ## F2S4_161010_003_D01        Sst
    ## F2S4_161010_003_E01        Sst
    ## F2S4_161010_003_F01        Sst
    ## F2S4_161010_003_G01        Sst
    ## F2S4_161010_003_H01        Sst
    ## F2S4_161010_004_A01        Sst
    ## F2S4_161010_004_B01        Sst
    ## F2S4_161010_004_C01        Sst
    ## F2S4_161010_004_D01        Sst
    ## F2S4_161010_004_F01        Sst
    ## F2S4_161010_004_G01        Sst
    ## F2S4_161010_004_H01        Sst
    ## F2S4_161010_005_B01        Sst
    ## F2S4_161010_005_C01        Sst
    ## F2S4_161010_005_D01        Sst
    ## F2S4_161010_005_E01        Sst
    ## F2S4_161010_005_F01        Sst
    ## F2S4_161010_005_G01        Sst
    ## F2S4_161010_005_H01        Sst
    ## F2S4_161010_006_A01        Sst
    ## F2S4_161010_006_B01        Sst
    ## F2S4_161010_006_C01        Sst
    ## F2S4_161010_006_D01        Sst
    ## F2S4_161010_006_E01        Sst
    ## F2S4_161010_006_F01        Sst
    ## F2S4_161010_006_G01        Sst
    ## F2S4_161010_006_H01        Sst
    ## F2S4_161010_007_A01        Sst
    ## F2S4_161010_007_B01        Sst
    ## F2S4_161010_007_C01        Sst
    ## F2S4_161010_007_D01        Sst
    ## F2S4_161010_007_E01      Pvalb
    ## F2S4_161010_007_F01        Sst
    ## F2S4_161010_007_G01        Sst
    ## F2S4_161010_007_H01        Sst
    ## F2S4_161010_008_A01        Sst
    ## F2S4_161010_008_B01        Sst
    ## F2S4_161010_008_C01        Sst
    ## F2S4_161010_008_D01        Sst
    ## F2S4_161010_008_E01        Sst
    ## F2S4_161010_008_F01        Sst
    ## F2S4_161010_008_G01        Sst
    ## F2S4_161010_008_H01        Sst
    ## F2S4_161010_009_A01        Sst
    ## F2S4_161010_009_B01        Sst
    ## F2S4_161010_009_C01        Sst
    ## F2S4_161010_009_D01        Sst
    ## F2S4_161010_009_E01        Sst
    ## F2S4_161010_009_F01        Sst
    ## F2S4_161010_009_G01        Sst
    ## F2S4_161010_009_H01        Sst
    ## F2S4_161010_010_A01        Sst
    ## F2S4_161010_010_B01        Sst
    ## F2S4_161010_010_C01        Sst
    ## F2S4_161010_010_D01        Sst
    ## F2S4_161010_010_E01      Pvalb
    ## F2S4_161010_010_F01        Sst
    ## F2S4_161010_010_G01        Sst
    ## F2S4_161010_010_H01        Sst
    ## F2S4_161010_011_A01        Sst
    ## F2S4_161010_011_B01        Sst
    ## F2S4_161010_011_C01      Pvalb
    ## F2S4_161010_011_D01        Sst
    ## F2S4_161010_011_E01        Sst
    ## F2S4_161010_011_F01        Sst
    ## F2S4_161010_011_G01        Sst
    ## F2S4_161010_011_H01        Sst
    ## F2S4_161010_012_A01        Sst
    ## F2S4_161010_012_B01        Sst
    ## F2S4_161010_012_C01        Sst
    ## F2S4_161010_012_D01        Sst
    ## F2S4_161010_012_E01        Sst
    ## F2S4_161010_012_G01        Sst
    ## F2S4_161010_012_H01        Sst
    ## F2S4_161010_013_A01      Pvalb
    ## F2S4_161010_013_C01      Pvalb
    ## F2S4_161010_013_D01        Sst
    ## F2S4_161010_013_E01        Sst
    ## F2S4_161010_013_F01        Sst
    ## F2S4_161010_013_H01        Sst
    ## F2S4_161010_014_A01        Sst
    ## F2S4_161010_014_B01        Sst
    ## F2S4_161010_014_C01        Sst
    ## F2S4_161010_014_D01        Sst
    ## F2S4_161010_014_E01        Sst
    ## F2S4_161010_014_F01        Sst
    ## F2S4_161010_014_G01        Sst
    ## F2S4_161010_014_H01        Sst
    ## F2S4_161010_015_A01        Sst
    ## F2S4_161010_015_B01        Sst
    ## F2S4_161010_015_C01      Pvalb
    ## F2S4_161010_015_D01        Sst
    ## F2S4_161010_015_E01        Sst
    ## F2S4_161010_015_F01        Sst
    ## F2S4_161010_015_G01        Sst
    ## F2S4_161010_015_H01        Sst
    ## F2S4_161010_016_A01      Pvalb
    ## F2S4_161010_016_B01        Sst
    ## F2S4_161010_016_C01        Sst
    ## F2S4_161010_016_D01      Pvalb
    ## F2S4_161010_016_E01        Sst
    ## F2S4_161010_016_F01        Sst
    ## F2S4_161010_016_G01        Sst
    ## F2S4_161010_016_H01      Pvalb
    ## F2S4_161010_017_A01      Pvalb
    ## F2S4_161010_017_B01        Sst
    ## F2S4_161010_017_C01        Sst
    ## F2S4_161010_017_D01        Sst
    ## F2S4_161010_017_E01        Sst
    ## F2S4_161010_017_F01        Sst
    ## F2S4_161010_017_G01        Sst
    ## F2S4_161010_017_H01        Sst
    ## F2S4_161010_018_A01        Sst
    ## F2S4_161010_018_B01        Sst
    ## F2S4_161010_018_C01        Sst
    ## F2S4_161010_018_D01        Sst
    ## F2S4_161010_018_F01      Pvalb
    ## F2S4_161010_018_G01        Sst
    ## F2S4_161010_018_H01        Sst
    ## F2S4_161010_019_A01        Sst
    ## F2S4_161010_019_B01        Sst
    ## F2S4_161010_019_C01      Pvalb
    ## F2S4_161010_019_D01        Sst
    ## F2S4_161010_019_E01        Sst
    ## F2S4_161010_019_F01      Pvalb
    ## F2S4_161010_019_G01      Pvalb
    ## F2S4_161010_019_H01      Pvalb
    ## F2S4_161010_020_A01        Sst
    ## F2S4_161010_020_B01        Sst
    ## F2S4_161010_020_C01      Pvalb
    ## F2S4_161010_020_D01        Sst
    ## F2S4_161010_020_E01      Pvalb
    ## F2S4_161010_020_F01        Sst
    ## F2S4_161010_020_G01        Sst
    ## F2S4_161010_020_H01        Sst
    ## F2S4_161010_021_A01        Sst
    ## F2S4_161010_021_B01        Sst
    ## F2S4_161010_021_C01        Sst
    ## F2S4_161010_021_D01        Sst
    ## F2S4_161010_021_E01        Sst
    ## F2S4_161010_021_F01        Sst
    ## F2S4_161010_021_G01        Sst
    ## F2S4_161010_021_H01        Sst
    ## F2S4_161010_022_A01        Sst
    ## F2S4_161010_022_B01        Sst
    ## F2S4_161010_022_C01        Sst
    ## F2S4_161010_022_D01        Sst
    ## F2S4_161010_022_E01        Sst
    ## F2S4_161010_022_F01        Sst
    ## F2S4_161010_022_G01      Pvalb
    ## F2S4_161010_022_H01        Sst
    ## F2S4_161010_023_A01        Sst
    ## F2S4_161010_023_B01        Sst
    ## F2S4_161010_023_C01        Sst
    ## F2S4_161010_023_D01        Sst
    ## F2S4_161010_023_E01        Sst
    ## F2S4_161010_023_F01        Sst
    ## F2S4_161010_023_G01        Sst
    ## F2S4_161010_023_H01        Sst
    ## F2S4_161010_024_B01        Sst
    ## F2S4_161010_024_C01        Sst
    ## F2S4_161010_024_D01        Sst
    ## F2S4_161010_024_E01        Sst
    ## F2S4_161010_024_F01        Sst
    ## F2S4_161010_024_G01        Sst
    ## F2S4_161010_024_H01        Sst
    ## F2S4_161010_025_A01        Sst
    ## F2S4_161010_025_B01        Sst
    ## F2S4_161010_025_C01        Sst
    ## F2S4_161010_025_D01        Sst
    ## F2S4_161010_025_E01        Sst
    ## F2S4_161010_025_F01        Sst
    ## F2S4_161010_025_G01        Sst
    ## F2S4_161010_025_H01        Sst
    ## F2S4_161010_026_A01        Sst
    ## F2S4_161010_026_B01        Sst
    ## F2S4_161010_026_C01        Sst
    ## F2S4_161010_026_D01        Sst
    ## F2S4_161010_026_E01        Sst
    ## F2S4_161010_026_F01        Sst
    ## F2S4_161010_026_G01        Sst
    ## F2S4_161010_026_H01        Sst
    ## F2S4_161010_027_A01        Sst
    ## F2S4_161010_027_B01        Sst
    ## F2S4_161010_027_C01        Sst
    ## F2S4_161010_027_D01        Sst
    ## F2S4_161010_027_E01        Sst
    ## F2S4_161010_027_F01        Sst
    ## F2S4_161010_027_G01        Sst
    ## F2S4_161010_027_H01        Sst
    ## F2S4_161010_028_A01        Sst
    ## F2S4_161010_028_B01        Sst
    ## F2S4_161010_028_C01        Sst
    ## F2S4_161010_028_D01      Pvalb
    ## F2S4_161010_028_E01        Sst
    ## F2S4_161010_028_F01        Sst
    ## F2S4_161010_028_G01        Sst
    ## F2S4_161010_028_H01        Sst
    ## F2S4_161010_029_A01        Sst
    ## F2S4_161010_029_B01        Sst
    ## F2S4_161010_029_C01        Sst
    ## F2S4_161010_029_D01        Sst
    ## F2S4_161010_029_E01        Sst
    ## F2S4_161010_029_F01      Pvalb
    ## F2S4_161010_029_G01        Sst
    ## F2S4_161010_029_H01        Sst
    ## F2S4_161010_030_A01        Sst
    ## F2S4_161010_030_B01        Sst
    ## F2S4_161010_030_C01        Sst
    ## F2S4_161010_030_D01        Sst
    ## F2S4_161010_030_E01        Sst
    ## F2S4_161010_030_F01        Sst
    ## F2S4_161010_030_G01        Sst
    ## F2S4_161010_030_H01      Pvalb
    ## F2S4_161010_031_A01        Sst
    ## F2S4_161010_031_B01      Pvalb
    ## F2S4_161010_031_C01        Sst
    ## F2S4_161010_031_D01        Sst
    ## F2S4_161010_031_E01        Sst
    ## F2S4_161010_031_F01        Sst
    ## F2S4_161010_031_G01        Sst
    ## F2S4_161010_031_H01        Sst
    ## F2S4_161010_032_A01        Sst
    ## F2S4_161010_032_B01        Sst
    ## F2S4_161010_032_C01        Sst
    ## F2S4_161010_032_D01        Sst
    ## F2S4_161010_032_E01        Sst
    ## F2S4_161010_032_F01        Sst
    ## F2S4_161010_032_G01        Vip
    ## F2S4_161010_032_H01        Sst
    ## F2S4_161010_033_B01        Sst
    ## F2S4_161010_033_C01        Sst
    ## F2S4_161010_033_D01        Sst
    ## F2S4_161010_033_E01        Sst
    ## F2S4_161010_033_F01        Sst
    ## F2S4_161010_033_G01        Sst
    ## F2S4_161010_033_H01        Sst
    ## F2S4_161010_034_G01        Sst
    ## F2S4_161010_034_H01        Sst
    ## F2S4_161017_001_A01        Vip
    ## F2S4_161017_001_B01        Vip
    ## F2S4_161017_001_C01        Vip
    ## F2S4_161017_001_D01        Vip
    ## F2S4_161017_001_E01        Vip
    ## F2S4_161017_001_F01        Vip
    ## F2S4_161017_001_G01        Vip
    ## F2S4_161017_001_H01        Vip
    ## F2S4_161017_002_A01        Vip
    ## F2S4_161017_002_B01        Vip
    ## F2S4_161017_002_C01        Vip
    ## F2S4_161017_002_E01        Vip
    ## F2S4_161017_002_F01        Vip
    ## F2S4_161017_002_G01        Vip
    ## F2S4_161017_002_H01        Vip
    ## F2S4_161017_003_A01        Vip
    ## F2S4_161017_003_B01       Sncg
    ## F2S4_161017_003_C01        Vip
    ## F2S4_161017_003_D01        Vip
    ## F2S4_161017_003_E01        Vip
    ## F2S4_161017_003_F01        Vip
    ## F2S4_161017_003_G01        Vip
    ## F2S4_161017_003_H01        Vip
    ## F2S4_161017_004_A01        Vip
    ## F2S4_161017_004_C01        Vip
    ## F2S4_161017_004_D01        Vip
    ## F2S4_161017_004_E01       Sncg
    ## F2S4_161017_004_F01        Vip
    ## F2S4_161017_004_G01        Vip
    ## F2S4_161017_004_H01        Vip
    ## F2S4_161017_005_A01        Vip
    ## F2S4_161017_005_B01        Vip
    ## F2S4_161017_005_C01        Vip
    ## F2S4_161017_005_D01        Vip
    ## F2S4_161017_005_E01        Vip
    ## F2S4_161017_005_F01        Vip
    ## F2S4_161017_005_G01        Vip
    ## F2S4_161017_005_H01        Vip
    ## F2S4_161017_006_A01        Vip
    ## F2S4_161017_006_B01        Vip
    ## F2S4_161017_006_C01        Vip
    ## F2S4_161017_006_D01        Vip
    ## F2S4_161017_006_E01        Vip
    ## F2S4_161017_006_F01        Vip
    ## F2S4_161017_006_G01        Vip
    ## F2S4_161017_006_H01        Vip
    ## F2S4_161017_007_A01        Vip
    ## F2S4_161017_007_B01        Vip
    ## F2S4_161017_007_C01        Vip
    ## F2S4_161017_007_D01        Vip
    ## F2S4_161017_007_E01        Vip
    ## F2S4_161017_007_F01        Vip
    ## F2S4_161017_007_G01        Vip
    ## F2S4_161017_007_H01       Sncg
    ## F2S4_161017_008_A01        Vip
    ## F2S4_161017_008_B01       Sncg
    ## F2S4_161017_008_C01        Vip
    ## F2S4_161017_008_D01        Vip
    ## F2S4_161017_008_E01        Vip
    ## F2S4_161017_008_F01        Vip
    ## F2S4_161017_008_G01        Vip
    ## F2S4_161017_008_H01        Vip
    ## F2S4_161017_009_A01        Vip
    ## F2S4_161017_009_B01        Vip
    ## F2S4_161017_009_C01        Vip
    ## F2S4_161017_009_D01        Vip
    ## F2S4_161017_009_E01        Vip
    ## F2S4_161017_009_F01        Vip
    ## F2S4_161017_009_G01        Vip
    ## F2S4_161017_009_H01        Vip
    ## F2S4_161017_010_A01        Vip
    ## F2S4_161017_010_B01        Vip
    ## F2S4_161017_010_C01        Vip
    ## F2S4_161017_010_D01        Vip
    ## F2S4_161017_010_E01        Vip
    ## F2S4_161017_010_F01        Vip
    ## F2S4_161017_010_G01        Vip
    ## F2S4_161017_010_H01        Vip
    ## F2S4_161017_011_A01        Vip
    ## F2S4_161017_011_B01        Vip
    ## F2S4_161017_011_C01        Vip
    ## F2S4_161017_011_D01        Vip
    ## F2S4_161017_011_E01        Vip
    ## F2S4_161017_011_F01        Vip
    ## F2S4_161017_011_G01        Vip
    ## F2S4_161017_011_H01        Vip
    ## F2S4_161017_012_A01        Vip
    ## F2S4_161017_012_B01        Vip
    ## F2S4_161017_012_C01        Vip
    ## F2S4_161017_012_D01        Vip
    ## F2S4_161017_012_E01        Vip
    ## F2S4_161017_012_F01       Sncg
    ## F2S4_161017_012_G01        Vip
    ## F2S4_161017_012_H01        Vip
    ## F2S4_161017_013_A01        Vip
    ## F2S4_161017_013_B01        Vip
    ## F2S4_161017_013_C01        Vip
    ## F2S4_161017_013_D01        Vip
    ## F2S4_161017_013_E01        Vip
    ## F2S4_161017_013_F01        Vip
    ## F2S4_161017_013_G01        Vip
    ## F2S4_161017_013_H01        Vip
    ## F2S4_161017_014_A01        Vip
    ## F2S4_161017_014_B01        Vip
    ## F2S4_161017_014_C01        Vip
    ## F2S4_161017_014_D01        Vip
    ## F2S4_161017_014_E01        Vip
    ## F2S4_161017_014_F01        Vip
    ## F2S4_161017_014_G01        Vip
    ## F2S4_161017_014_H01        Vip
    ## F2S4_161017_015_A01        Vip
    ## F2S4_161017_015_B01        Vip
    ## F2S4_161017_015_C01        Vip
    ## F2S4_161017_015_D01        Vip
    ## F2S4_161017_015_E01        Vip
    ## F2S4_161017_015_F01        Vip
    ## F2S4_161017_015_G01        Vip
    ## F2S4_161017_015_H01        Vip
    ## F2S4_161017_016_A01        Vip
    ## F2S4_161017_016_B01        Vip
    ## F2S4_161017_016_C01        Vip
    ## F2S4_161017_016_D01        Vip
    ## F2S4_161017_016_E01        Vip
    ## F2S4_161017_016_F01        Vip
    ## F2S4_161017_016_G01        Vip
    ## F2S4_161017_016_H01        Vip
    ## F2S4_161017_017_A01        Vip
    ## F2S4_161017_017_B01        Vip
    ## F2S4_161017_017_C01        Vip
    ## F2S4_161017_017_D01        Vip
    ## F2S4_161017_017_E01        Vip
    ## F2S4_161017_017_F01        Vip
    ## F2S4_161017_017_G01        Vip
    ## F2S4_161017_017_H01        Vip
    ## F2S4_161017_018_A01        Vip
    ## F2S4_161017_018_B01        Vip
    ## F2S4_161017_018_C01        Vip
    ## F2S4_161017_018_D01        Vip
    ## F2S4_161017_018_E01        Vip
    ## F2S4_161017_018_F01        Vip
    ## F2S4_161017_018_G01        Vip
    ## F2S4_161017_018_H01        Vip
    ## F2S4_161017_019_A01        Vip
    ## F2S4_161017_019_B01        Vip
    ## F2S4_161017_019_C01        Vip
    ## F2S4_161017_019_D01        Vip
    ## F2S4_161017_019_E01        Vip
    ## F2S4_161017_019_F01        Vip
    ## F2S4_161017_019_G01        Vip
    ## F2S4_161017_019_H01        Vip
    ## F2S4_161017_020_A01        Vip
    ## F2S4_161017_020_B01        Vip
    ## F2S4_161017_020_C01        Vip
    ## F2S4_161017_020_D01        Vip
    ## F2S4_161017_020_E01        Vip
    ## F2S4_161017_020_F01        Vip
    ## F2S4_161017_020_G01        Vip
    ## F2S4_161017_020_H01        Vip
    ## F2S4_161017_021_A01        Vip
    ## F2S4_161017_021_B01        Vip
    ## F2S4_161017_021_C01        Vip
    ## F2S4_161017_021_D01        Vip
    ## F2S4_161017_021_E01        Vip
    ## F2S4_161017_021_F01        Vip
    ## F2S4_161017_021_G01        Vip
    ## F2S4_161017_021_H01   Serpinf1
    ## F2S4_161017_022_A01        Vip
    ## F2S4_161017_022_B01        Vip
    ## F2S4_161017_022_C01        Vip
    ## F2S4_161017_022_D01        Vip
    ## F2S4_161017_022_E01        Vip
    ## F2S4_161017_022_F01        Vip
    ## F2S4_161017_022_G01       Sncg
    ## F2S4_161017_023_A01        Vip
    ## F2S4_161017_023_B01        Vip
    ## F2S4_161017_023_C01        Vip
    ## F2S4_161017_023_D01        Vip
    ## F2S4_161017_023_E01        Vip
    ## F2S4_161017_023_F01        Vip
    ## F2S4_161017_023_G01        Vip
    ## F2S4_161017_023_H01        Vip
    ## F2S4_161017_024_A01        Vip
    ## F2S4_161017_024_B01        Vip
    ## F2S4_161017_024_C01       Sncg
    ## F2S4_161017_024_D01        Vip
    ## F2S4_161017_024_E01        Vip
    ## F2S4_161017_024_F01        Vip
    ## F2S4_161017_024_G01        Vip
    ## F2S4_161017_024_H01        Vip
    ## F2S4_161017_025_A01        Vip
    ## F2S4_161017_025_B01        Vip
    ## F2S4_161017_025_C01        Vip
    ## F2S4_161017_025_D01        Vip
    ## F2S4_161017_025_E01   Serpinf1
    ## F2S4_161017_025_F01        Vip
    ## F2S4_161017_025_H01        Vip
    ## F2S4_161017_026_A01        Vip
    ## F2S4_161017_026_B01        Vip
    ## F2S4_161017_026_C01        Vip
    ## F2S4_161017_026_D01        Vip
    ## F2S4_161017_026_E01        Vip
    ## F2S4_161017_026_F01        Vip
    ## F2S4_161017_026_G01        Vip
    ## F2S4_161017_026_H01        Vip
    ## F2S4_161018_001_A01         L4
    ## F2S4_161018_001_C01         NP
    ## F2S4_161018_001_D01      L5 PT
    ## F2S4_161018_001_E01         NP
    ## F2S4_161018_001_F01      L5 IT
    ## F2S4_161018_001_G01      L5 PT
    ## F2S4_161018_001_H01      L5 IT
    ## F2S4_161018_002_B01         NP
    ## F2S4_161018_002_D01         NP
    ## F2S4_161018_002_F01      L5 PT
    ## F2S4_161018_002_G01         L4
    ## F2S4_161018_003_D01         NP
    ## F2S4_161018_003_E01         NP
    ## F2S4_161018_003_F01      L5 PT
    ## F2S4_161018_003_H01      L5 PT
    ## F2S4_161020_009_A01        Vip
    ## F2S4_161020_009_B01        Vip
    ## F2S4_161020_009_C01        Vip
    ## F2S4_161020_009_D01        Vip
    ## F2S4_161020_009_E01        Vip
    ## F2S4_161020_009_F01        Vip
    ## F2S4_161020_009_H01        Vip
    ## F2S4_161020_010_A01        Vip
    ## F2S4_161020_010_B01        Vip
    ## F2S4_161020_010_C01        Vip
    ## F2S4_161020_010_D01        Vip
    ## F2S4_161020_010_E01        Vip
    ## F2S4_161020_010_F01      Pvalb
    ## F2S4_161020_010_G01      Pvalb
    ## F2S4_161020_011_A01      Pvalb
    ## F2S4_161020_011_B01      Pvalb
    ## F2S4_161020_011_D01        Vip
    ## F2S4_161020_011_E01        Vip
    ## F2S4_161020_011_F01        Vip
    ## F2S4_161020_011_G01        Vip
    ## F2S4_161020_011_H01        Vip
    ## F2S4_161020_012_B01        Vip
    ## F2S4_161020_012_C01        Vip
    ## F2S4_161020_012_D01        Vip
    ## F2S4_161020_012_E01        Vip
    ## F2S4_161020_012_F01        Vip
    ## F2S4_161020_012_G01      Pvalb
    ## F2S4_161020_012_H01      Pvalb
    ## F2S4_161020_013_A01        Sst
    ## F2S4_161028_001_A01      Pvalb
    ## F2S4_161028_001_B01      Pvalb
    ## F2S4_161028_001_C01      Pvalb
    ## F2S4_161028_001_D01        Sst
    ## F2S4_161028_001_E01      Pvalb
    ## F2S4_161028_001_F01        Sst
    ## F2S4_161028_001_G01        Sst
    ## F2S4_161028_001_H01        Sst
    ## F2S4_161028_002_A01        Sst
    ## F2S4_161028_002_B01        Sst
    ## F2S4_161028_002_C01        Sst
    ## F2S4_161028_002_D01        Sst
    ## F2S4_161028_002_E01      Pvalb
    ## F2S4_161028_002_F01        Sst
    ## F2S4_161028_002_G01        Sst
    ## F2S4_161028_002_H01      Pvalb
    ## F2S4_161028_003_A01        Sst
    ## F2S4_161028_003_B01        Sst
    ## F2S4_161028_003_C01        Sst
    ## F2S4_161028_003_D01        Sst
    ## F2S4_161028_003_E01        Sst
    ## F2S4_161028_003_F01        Sst
    ## F2S4_161028_003_G01        Sst
    ## F2S4_161028_003_H01      Pvalb
    ## F2S4_161028_004_A01        Sst
    ## F2S4_161028_004_B01        Sst
    ## F2S4_161028_004_C01        Sst
    ## F2S4_161028_004_D01      Pvalb
    ## F2S4_161028_004_E01        Sst
    ## F2S4_161028_004_F01        Sst
    ## F2S4_161028_004_G01      Pvalb
    ## F2S4_161028_004_H01        Sst
    ## F2S4_161028_005_A01        Sst
    ## F2S4_161028_005_B01      Pvalb
    ## F2S4_161028_005_C01        Sst
    ## F2S4_161028_005_D01        Sst
    ## F2S4_161028_005_E01        Sst
    ## F2S4_161028_005_F01        Sst
    ## F2S4_161028_005_G01        Sst
    ## F2S4_161028_005_H01        Sst
    ## F2S4_161028_006_A01        Sst
    ## F2S4_161028_006_C01        Sst
    ## F2S4_161028_006_D01        Sst
    ## F2S4_161028_006_E01        Sst
    ## F2S4_161028_006_F01        Sst
    ## F2S4_161028_006_G01        Sst
    ## F2S4_161028_006_H01        Sst
    ## F2S4_161028_007_B01        Sst
    ## F2S4_161028_007_C01        Sst
    ## F2S4_161028_007_D01        Sst
    ## F2S4_161028_007_E01        Sst
    ## F2S4_161028_007_F01        Sst
    ## F2S4_161028_007_G01        Sst
    ## F2S4_161028_007_H01        Sst
    ## F2S4_161028_008_A01        Sst
    ## F2S4_161028_008_B01        Sst
    ## F2S4_161028_008_C01        Sst
    ## F2S4_161028_008_D01        Sst
    ## F2S4_161028_008_E01        Sst
    ## F2S4_161028_008_F01        Sst
    ## F2S4_161028_008_G01        Sst
    ## F2S4_161028_009_A01        Sst
    ## F2S4_161028_009_B01        Sst
    ## F2S4_161028_009_C01        Sst
    ## F2S4_161028_009_D01        Sst
    ## F2S4_161028_009_E01        Sst
    ## F2S4_161028_009_F01        Sst
    ## F2S4_161028_009_G01        Sst
    ## F2S4_161028_009_H01        Sst
    ## F2S4_161028_010_A01        Sst
    ## F2S4_161028_010_B01        Sst
    ## F2S4_161028_010_C01        Sst
    ## F2S4_161028_010_D01        Sst
    ## F2S4_161028_010_E01        Sst
    ## F2S4_161028_010_F01        Sst
    ## F2S4_161028_010_G01        Sst
    ## F2S4_161028_010_H01        Sst
    ## F2S4_161028_011_A01        Sst
    ## F2S4_161028_011_B01        Sst
    ## F2S4_161028_011_C01        Sst
    ## F2S4_161028_011_D01        Sst
    ## F2S4_161028_011_E01        Sst
    ## F2S4_161028_011_F01        Sst
    ## F2S4_161028_011_H01        Sst
    ## F2S4_161028_012_A01        Sst
    ## F2S4_161028_012_B01        Sst
    ## F2S4_161028_012_C01        Sst
    ## F2S4_161028_012_D01        Sst
    ## F2S4_161028_012_E01        Sst
    ## F2S4_161028_012_F01        Sst
    ## F2S4_161028_012_G01        Sst
    ## F2S4_161028_012_H01        Sst
    ## F2S4_161028_013_A01        Sst
    ## F2S4_161028_013_B01        Sst
    ## F2S4_161028_013_C01        Sst
    ## F2S4_161028_013_D01        Sst
    ## F2S4_161028_013_E01        Sst
    ## F2S4_161028_013_G01        Sst
    ## F2S4_161028_014_A01         L4
    ## F2S4_161028_014_B01         L4
    ## F2S4_161028_014_C01      L5 IT
    ## F2S4_161028_014_D01      L5 IT
    ## F2S4_161028_014_E01      L5 IT
    ## F2S4_161028_014_F01      L5 IT
    ## F2S4_161028_014_G01         L4
    ## F2S4_161028_014_H01         L4
    ## F2S4_161028_015_A01         L4
    ## F2S4_161028_015_C01      L5 IT
    ## F2S4_161028_015_D01      L5 IT
    ## F2S4_161028_015_E01      L5 IT
    ## F2S4_161028_015_F01         L4
    ## F2S4_161028_015_G01      L5 IT
    ## F2S4_161028_016_A01         L4
    ## F2S4_161028_016_B01         L4
    ## F2S4_161028_016_C01      L5 IT
    ## F2S4_161028_016_D01      L5 IT
    ## F2S4_161028_016_E01      L5 IT
    ## F2S4_161028_016_F01      L5 IT
    ## F2S4_161028_016_G01         L4
    ## F2S4_161028_016_H01      L5 IT
    ## F2S4_161028_017_A01      L5 IT
    ## F2S4_161028_017_B01      L5 IT
    ## F2S4_161028_017_D01      L5 IT
    ## F2S4_161028_017_E01      L5 IT
    ## F2S4_161028_017_F01      L5 IT
    ## F2S4_161028_017_G01         L4
    ## F2S4_161028_017_H01      L5 IT
    ## F2S4_161028_018_A01      L5 IT
    ## F2S4_161028_018_C01      L5 IT
    ## F2S4_161028_018_D01      L5 IT
    ## F2S4_161028_018_E01         L4
    ## F2S4_161028_018_F01      L5 IT
    ## F2S4_161028_018_H01      L5 IT
    ## F2S4_161028_019_A01      L5 IT
    ## F2S4_161028_019_B01         L4
    ## F2S4_161028_019_C01         L4
    ## F2S4_161028_019_D01      L5 IT
    ## F2S4_161028_019_G01      L5 IT
    ## F2S4_161028_019_H01         L4
    ## F2S4_161028_020_A01      L5 IT
    ## F2S4_161028_020_B01      L5 IT
    ## F2S4_161028_020_E01         L4
    ## F2S4_161028_020_G01         L4
    ## F2S4_161028_020_H01      L5 IT
    ## F2S4_161028_021_A01      L5 IT
    ## F2S4_161028_021_B01      L5 IT
    ## F2S4_161028_021_C01         L4
    ## F2S4_161028_021_F01      L5 IT
    ## F2S4_161028_021_G01      L5 IT
    ## F2S4_161028_021_H01         L4
    ## F2S4_161103_001_A01      L5 IT
    ## F2S4_161103_001_B01         NP
    ## F2S4_161103_001_C01      L6 CT
    ## F2S4_161103_001_D01      L5 IT
    ## F2S4_161103_001_E01      L5 IT
    ## F2S4_161103_001_F01      L5 PT
    ## F2S4_161103_001_G01      L5 IT
    ## F2S4_161103_001_H01         NP
    ## F2S4_161103_002_A01         NP
    ## F2S4_161103_002_B01         L4
    ## F2S4_161103_002_D01         NP
    ## F2S4_161103_002_F01         L4
    ## F2S4_161103_002_G01         L4
    ## F2S4_161103_002_H01      L5 IT
    ## F2S4_161103_003_A01      L5 IT
    ## F2S4_161103_003_D01      L5 IT
    ## F2S4_161103_003_E01         NP
    ## F2S4_161103_003_H01         L4
    ## F2S4_161103_004_A01      L5 PT
    ## F2S4_161103_004_B01         NP
    ## F2S4_161103_004_C01      L5 PT
    ## F2S4_161103_004_D01         L4
    ## F2S4_161103_004_E01      L5 IT
    ## F2S4_161103_004_H01      L5 PT
    ## F2S4_161103_005_B01         NP
    ## F2S4_161103_005_C01         L4
    ## F2S4_161103_005_G01         L4
    ## F2S4_161103_006_A01         NP
    ## F2S4_161103_006_B01         L4
    ## F2S4_161103_006_C01         NP
    ## F2S4_161103_006_D01      L5 IT
    ## F2S4_161103_006_E01         NP
    ## F2S4_161103_007_B01      L6 CT
    ## F2S4_161103_007_C01      L5 IT
    ## F2S4_161103_007_D01      L5 IT
    ## F2S4_161103_007_E01      L5 PT
    ## F2S4_161103_007_F01         NP
    ## F2S4_161103_007_G01      L6 CT
    ## F2S4_161103_007_H01      L5 PT
    ## F2S4_161103_008_A01        Sst
    ## F2S4_161103_008_B01      Lamp5
    ## F2S4_161103_008_D01      Pvalb
    ## F2S4_161103_008_E01      Pvalb
    ## F2S4_161103_008_F01      Pvalb
    ## F2S4_161103_008_G01      Pvalb
    ## F2S4_161107_001_A01        Vip
    ## F2S4_161107_001_B01        Vip
    ## F2S4_161107_001_C01        Vip
    ## F2S4_161107_001_D01        Vip
    ## F2S4_161107_001_F01        Vip
    ## F2S4_161107_001_G01        Sst
    ## F2S4_161107_001_H01        Vip
    ## F2S4_161107_002_A01        Vip
    ## F2S4_161107_002_B01        Vip
    ## F2S4_161107_002_C01        Vip
    ## F2S4_161107_002_D01        Vip
    ## F2S4_161107_002_E01        Vip
    ## F2S4_161107_002_F01        Vip
    ## F2S4_161107_002_G01        Vip
    ## F2S4_161107_002_H01        Vip
    ## F2S4_161107_003_A01        Vip
    ## F2S4_161107_003_B01        Vip
    ## F2S4_161107_003_C01        Vip
    ## F2S4_161107_003_D01        Vip
    ## F2S4_161107_003_E01        Vip
    ## F2S4_161107_003_F01        Vip
    ## F2S4_161107_003_G01        Vip
    ## F2S4_161107_003_H01        Vip
    ## F2S4_161107_004_A01        Vip
    ## F2S4_161107_004_B01        Vip
    ## F2S4_161107_004_C01        Vip
    ## F2S4_161107_004_D01        Vip
    ## F2S4_161107_004_E01        Vip
    ## F2S4_161107_004_F01        Vip
    ## F2S4_161107_004_G01       Sncg
    ## F2S4_161107_004_H01        Vip
    ## F2S4_161107_005_A01        Vip
    ## F2S4_161107_005_B01        Vip
    ## F2S4_161107_005_C01        Vip
    ## F2S4_161107_005_D01        Vip
    ## F2S4_161107_005_E01        Vip
    ## F2S4_161107_005_F01        Vip
    ## F2S4_161107_005_G01        Vip
    ## F2S4_161107_005_H01        Vip
    ## F2S4_161107_006_A01        Vip
    ## F2S4_161107_006_B01        Vip
    ## F2S4_161107_006_C01        Vip
    ## F2S4_161107_006_D01        Vip
    ## F2S4_161107_006_E01       Sncg
    ## F2S4_161107_006_F01        Vip
    ## F2S4_161107_006_G01        Vip
    ## F2S4_161107_006_H01        Vip
    ## F2S4_161107_007_A01       Sncg
    ## F2S4_161107_007_B01        Vip
    ## F2S4_161107_007_C01        Vip
    ## F2S4_161107_007_D01        Vip
    ## F2S4_161107_007_E01        Vip
    ## F2S4_161107_007_F01        Vip
    ## F2S4_161107_007_G01        Vip
    ## F2S4_161107_007_H01        Vip
    ## F2S4_161107_008_A01        Vip
    ## F2S4_161107_008_B01        Vip
    ## F2S4_161107_008_C01        Vip
    ## F2S4_161107_008_D01   Serpinf1
    ## F2S4_161107_008_E01        Vip
    ## F2S4_161107_008_F01        Vip
    ## F2S4_161107_008_G01        Vip
    ## F2S4_161107_008_H01        Vip
    ## F2S4_161107_009_A01        Vip
    ## F2S4_161107_009_B01        Vip
    ## F2S4_161107_009_C01        Vip
    ## F2S4_161107_009_D01        Vip
    ## F2S4_161107_009_E01       Sncg
    ## F2S4_161107_009_F01        Vip
    ## F2S4_161107_009_G01        Vip
    ## F2S4_161107_009_H01        Vip
    ## F2S4_161107_010_A01        Vip
    ## F2S4_161107_010_B01        Vip
    ## F2S4_161107_010_C01        Vip
    ## F2S4_161107_010_D01        Vip
    ## F2S4_161107_010_E01       Sncg
    ## F2S4_161107_010_F01        Vip
    ## F2S4_161107_010_G01        Vip
    ## F2S4_161107_010_H01        Vip
    ## F2S4_161107_011_A01        Vip
    ## F2S4_161107_011_B01        Vip
    ## F2S4_161107_011_C01        Vip
    ## F2S4_161107_011_D01        Vip
    ## F2S4_161107_011_E01        Vip
    ## F2S4_161107_011_F01        Vip
    ## F2S4_161107_011_G01        Vip
    ## F2S4_161107_011_H01        Vip
    ## F2S4_161107_012_A01        Vip
    ## F2S4_161107_012_B01        Vip
    ## F2S4_161107_012_C01        Vip
    ## F2S4_161107_012_D01        Vip
    ## F2S4_161107_012_E01        Vip
    ## F2S4_161107_012_F01        Vip
    ## F2S4_161107_012_G01        Vip
    ## F2S4_161107_012_H01        Vip
    ## F2S4_161107_013_A01        Vip
    ## F2S4_161107_013_B01        Vip
    ## F2S4_161107_013_C01        Vip
    ## F2S4_161107_013_D01        Vip
    ## F2S4_161107_013_E01        Vip
    ## F2S4_161107_013_F01       Sncg
    ## F2S4_161107_013_G01        Vip
    ## F2S4_161107_013_H01        Vip
    ## F2S4_161107_014_A01        Vip
    ## F2S4_161107_014_B01        Vip
    ## F2S4_161107_014_D01        Vip
    ## F2S4_161107_014_E01        Vip
    ## F2S4_161107_014_F01        Vip
    ## F2S4_161107_014_G01        Vip
    ## F2S4_161107_014_H01        Vip
    ## F2S4_161107_015_B01        Vip
    ## F2S4_161107_015_C01        Vip
    ## F2S4_161107_015_D01        Vip
    ## F2S4_161107_015_E01        Vip
    ## F2S4_161107_015_F01        Vip
    ## F2S4_161107_015_G01        Vip
    ## F2S4_161107_015_H01        Vip
    ## F2S4_161107_016_A01        Vip
    ## F2S4_161107_016_B01        Vip
    ## F2S4_161107_016_C01        Vip
    ## F2S4_161107_016_D01        Vip
    ## F2S4_161107_016_E01        Vip
    ## F2S4_161107_016_F01        Vip
    ## F2S4_161107_016_G01        Vip
    ## F2S4_161107_016_H01        Vip
    ## F2S4_161107_017_A01        Vip
    ## F2S4_161107_017_B01        Vip
    ## F2S4_161107_017_C01        Vip
    ## F2S4_161107_017_D01        Vip
    ## F2S4_161107_017_E01        Vip
    ## F2S4_161107_017_F01        Vip
    ## F2S4_161107_017_G01        Vip
    ## F2S4_161107_017_H01        Vip
    ## F2S4_161107_018_A01        Vip
    ## F2S4_161107_018_B01        Vip
    ## F2S4_161107_018_C01        Vip
    ## F2S4_161107_018_D01        Vip
    ## F2S4_161107_018_E01        Vip
    ## F2S4_161107_018_F01        Vip
    ## F2S4_161107_018_G01        Vip
    ## F2S4_161107_018_H01        Vip
    ## F2S4_161109_001_B01      L5 IT
    ## F2S4_161109_001_C01      L5 IT
    ## F2S4_161109_001_D01        Vip
    ## F2S4_161109_001_E01        Vip
    ## F2S4_161109_001_F01         NP
    ## F2S4_161109_001_H01      Lamp5
    ## F2S4_161109_002_A01      Pvalb
    ## F2S4_161109_002_B01         NP
    ## F2S4_161109_002_C01      L5 IT
    ## F2S4_161110_001_A01        Sst
    ## F2S4_161110_001_B01      Lamp5
    ## F2S4_161110_001_C01       VLMC
    ## F2S4_161110_001_D01        L6b
    ## F2S4_161110_001_E01       VLMC
    ## F2S4_161110_001_F01       VLMC
    ## F2S4_161110_001_G01        Sst
    ## F2S4_161110_001_H01       VLMC
    ## F2S4_161110_002_A01       VLMC
    ## F2S4_161110_002_B01        Sst
    ## F2S4_161110_002_D01       VLMC
    ## F2S4_161110_002_E01       VLMC
    ## F2S4_161110_002_F01       VLMC
    ## F2S4_161110_002_G01       VLMC
    ## F2S4_161110_002_H01       VLMC
    ## F2S4_161110_003_C01       VLMC
    ## F2S4_161110_003_D01       VLMC
    ## F2S4_161110_003_E01       VLMC
    ## F2S4_161110_003_F01      L5 IT
    ## F2S4_161110_003_G01      L5 IT
    ## F2S4_161111_001_A01      L5 IT
    ## F2S4_161111_001_B01      L6 IT
    ## F2S4_161111_001_C01      L5 IT
    ## F2S4_161111_001_D01      L5 IT
    ## F2S4_161111_001_E01      L5 IT
    ## F2S4_161111_001_F01      L6 CT
    ## F2S4_161111_001_G01      L5 PT
    ## F2S4_161111_001_H01      L5 PT
    ## F2S4_161111_002_A01      L5 PT
    ## F2S4_161111_002_C01      L5 PT
    ## F2S4_161111_002_D01      L6 CT
    ## F2S4_161111_002_F01      L6 CT
    ## F2S4_161111_002_G01      L5 PT
    ## F2S4_161111_002_H01      L5 PT
    ## F2S4_161111_003_A01      L5 PT
    ## F2S4_161111_003_B01      L5 PT
    ## F2S4_161111_003_C01      L6 CT
    ## F2S4_161111_003_D01      L6 CT
    ## F2S4_161111_003_E01      L6 CT
    ## F2S4_161111_003_F01      L5 PT
    ## F2S4_161111_003_G01      L6 CT
    ## F2S4_161111_004_A01      L5 PT
    ## F2S4_161111_004_B01      L6 CT
    ## F2S4_161111_004_C01      L6 CT
    ## F2S4_161111_004_D01      L6 CT
    ## F2S4_161111_004_E01      L5 PT
    ## F2S4_161111_004_F01      L5 PT
    ## F2S4_161111_004_G01      L6 CT
    ## F2S4_161111_004_H01      L6 CT
    ## F2S4_161111_005_A01      L5 PT
    ## F2S4_161111_005_B01      L6 CT
    ## F2S4_161111_005_C01      L5 PT
    ## F2S4_161111_005_D01      L5 PT
    ## F2S4_161111_005_E01      L5 PT
    ## F2S4_161111_005_F01      L6 CT
    ## F2S4_161111_005_G01      L6 CT
    ## F2S4_161111_005_H01      L5 PT
    ## F2S4_161111_006_A01      L5 PT
    ## F2S4_161111_006_B01      L5 PT
    ## F2S4_161111_006_C01      L5 PT
    ## F2S4_161111_006_D01      L5 PT
    ## F2S4_161111_006_E01      L5 PT
    ## F2S4_161111_006_F01      L6 CT
    ## F2S4_161111_006_G01      L6 CT
    ## F2S4_161111_006_H01      L5 PT
    ## F2S4_161111_007_G01      L5 PT
    ## F2S4_161111_007_H01      L5 PT
    ## F2S4_161114_001_A01      L6 CT
    ## F2S4_161114_001_B01      L6 CT
    ## F2S4_161114_001_C01      L6 CT
    ## F2S4_161114_001_D01      L6 CT
    ## F2S4_161114_001_E01      L6 CT
    ## F2S4_161114_001_F01      L6 CT
    ## F2S4_161114_001_G01      L5 PT
    ## F2S4_161114_001_H01      L6 CT
    ## F2S4_161114_002_B01      L6 CT
    ## F2S4_161114_002_C01      L6 CT
    ## F2S4_161114_002_D01      L6 CT
    ## F2S4_161114_002_E01      L6 CT
    ## F2S4_161114_002_F01      L6 CT
    ## F2S4_161114_002_G01      L6 CT
    ## F2S4_161114_002_H01      L6 CT
    ## F2S4_161115_001_B01      L5 PT
    ## F2S4_161116_001_A01    L2/3 IT
    ## F2S4_161116_001_B01      L6 IT
    ## F2S4_161116_001_C01    L2/3 IT
    ## F2S4_161116_001_D01    L2/3 IT
    ## F2S4_161116_001_F01        L6b
    ## F2S4_161116_001_G01    L2/3 IT
    ## F2S4_161116_001_H01      L5 IT
    ## F2S4_161116_002_B01      L6 IT
    ## F2S4_161116_002_G01        L6b
    ## F2S4_161116_002_H01      L6 IT
    ## F2S4_161117_001_A01      L5 IT
    ## F2S4_161117_001_B01    L2/3 IT
    ## F2S4_161117_001_E01        L6b
    ## F2S4_161117_001_F01        L6b
    ## F2S4_161117_001_G01       VLMC
    ## F2S4_161117_002_B01      L5 IT
    ## F2S4_161117_002_C01      L5 IT
    ## F2S4_161117_002_D01    L2/3 IT
    ## F2S4_161117_002_F01        L6b
    ## F2S4_161117_002_G01    L2/3 IT
    ## F2S4_161117_002_H01       VLMC
    ## F2S4_161117_003_A01    L2/3 IT
    ## F2S4_161117_003_B01        L6b
    ## F2S4_161117_003_C01      L5 IT
    ## F2S4_161117_003_D01      L5 IT
    ## F2S4_161117_003_E01      L5 IT
    ## F2S4_161117_003_F01        Sst
    ## F2S4_161117_003_G01        L6b
    ## F2S4_161117_003_H01    L2/3 IT
    ## F2S4_161117_004_A01        L6b
    ## F2S4_161117_004_B01      L5 IT
    ## F2S4_161117_004_D01       VLMC
    ## F2S4_161117_004_E01      L5 IT
    ## F2S4_161117_004_F01        L6b
    ## F2S4_161117_004_G01    L2/3 IT
    ## F2S4_161117_004_H01       VLMC
    ## F2S4_161117_005_A01         L4
    ## F2S4_161117_005_B01      L5 IT
    ## F2S4_161117_005_C01         L4
    ## F2S4_161117_005_D01         L4
    ## F2S4_161117_005_E01    L2/3 IT
    ## F2S4_161117_005_F01    L2/3 IT
    ## F2S4_161117_005_G01        L6b
    ## F2S4_161117_005_H01      L5 IT
    ## F2S4_161117_006_A01      L5 IT
    ## F2S4_161117_006_B01      L5 IT
    ## F2S4_161117_006_C01      L6 IT
    ## F2S4_161117_006_D01    L2/3 IT
    ## F2S4_161117_006_E01        L6b
    ## F2S4_161117_006_F01        L6b
    ## F2S4_161117_006_G01      L5 IT
    ## F2S4_161117_006_H01         L4
    ## F2S4_161117_007_A01        L6b
    ## F2S4_161121_001_A01      L6 CT
    ## F2S4_161121_001_B01      L6 CT
    ## F2S4_161121_001_C01      L6 CT
    ## F2S4_161121_001_D01      L6 CT
    ## F2S4_161121_001_F01      L6 CT
    ## F2S4_161121_001_G01      L6 CT
    ## F2S4_161121_001_H01      L6 CT
    ## F2S4_161121_002_A01        L6b
    ## F2S4_161121_002_B01      L6 CT
    ## F2S4_161121_002_C01      L6 CT
    ## F2S4_161121_002_D01      L6 CT
    ## F2S4_161121_002_E01      L6 CT
    ## F2S4_161121_002_F01        L6b
    ## F2S4_161121_002_G01      L6 CT
    ## F2S4_161121_002_H01      L6 CT
    ## F2S4_161121_003_A01      L6 CT
    ## F2S4_161121_003_B01        L6b
    ## F2S4_161121_003_D01      L6 CT
    ## F2S4_161121_003_E01        L6b
    ## F2S4_161121_003_G01        L6b
    ## F2S4_161121_003_H01      L6 CT
    ## F2S4_161121_004_A01      L6 CT
    ## F2S4_161121_004_B01      L6 CT
    ## F2S4_161121_004_C01      L6 CT
    ## F2S4_161121_004_D01      L6 CT
    ## F2S4_161121_004_E01      L6 CT
    ## F2S4_161121_004_F01      L5 PT
    ## F2S4_161121_004_G01      L6 CT
    ## F2S4_161121_004_H01      L6 CT
    ## F2S4_161121_005_A01      L6 CT
    ## F2S4_161121_005_B01      L6 CT
    ## F2S4_161121_005_C01      L6 CT
    ## F2S4_161121_005_D01      L6 CT
    ## F2S4_161121_005_F01      L6 CT
    ## F2S4_161121_005_G01      L6 CT
    ## F2S4_161121_005_H01      L6 CT
    ## F2S4_161121_006_D01      L6 CT
    ## F2S4_161121_006_E01      L6 CT
    ## F2S4_161121_006_F01        L6b
    ## F2S4_161121_006_G01      L6 CT
    ## F2S4_161121_006_H01      L6 CT
    ## F2S4_161122_001_A01      L6 CT
    ## F2S4_161122_001_B01      L5 PT
    ## F2S4_161122_001_C01      L6 CT
    ## F2S4_161122_001_D01      L6 CT
    ## F2S4_161122_001_E01      L5 PT
    ## F2S4_161122_001_G01        L6b
    ## F2S4_161122_001_H01      L6 CT
    ## F2S4_161122_002_A01      L6 CT
    ## F2S4_161122_002_B01      L6 CT
    ## F2S4_161122_002_C01      L5 PT
    ## F2S4_161122_002_D01      L6 CT
    ## F2S4_161122_002_E01      L6 CT
    ## F2S4_161122_002_F01      L6 CT
    ## F2S4_161122_002_G01      L6 CT
    ## F2S4_161122_002_H01        L6b
    ## F2S4_161122_003_A01      L6 CT
    ## F2S4_161122_003_B01      L5 PT
    ## F2S4_161122_003_C01      L5 PT
    ## F2S4_161122_003_D01      L6 CT
    ## F2S4_161122_003_E01        L6b
    ## F2S4_161122_003_F01      L6 CT
    ## F2S4_161122_003_G01      L6 CT
    ## F2S4_161122_003_H01      L6 CT
    ## F2S4_161122_004_H01      L5 PT
    ## F2S4_161123_162_A01      L6 IT
    ## F2S4_161123_162_B01      L6 IT
    ## F2S4_161123_162_C01    L2/3 IT
    ## F2S4_161123_162_D01    L2/3 IT
    ## F2S4_161123_162_F01      L6 IT
    ## F2S4_161123_162_H01      L6 IT
    ## F2S4_161123_163_A01      L6 IT
    ## F2S4_161123_163_C01    L2/3 IT
    ## F2S4_161123_163_D01      L6 IT
    ## F2S4_161123_163_E01      L6 IT
    ## F2S4_161123_163_F01    L2/3 IT
    ## F2S4_161123_163_G01      L6 IT
    ## F2S4_161123_163_H01      L6 IT
    ## F2S4_161123_164_A01        L6b
    ## F2S4_161123_164_B01      L6 IT
    ## F2S4_161123_164_C01      L6 IT
    ## F2S4_161123_164_D01      L6 IT
    ## F2S4_161123_164_E01    L2/3 IT
    ## F2S4_161128_001_A01      L5 IT
    ## F2S4_161128_001_B01      L5 IT
    ## F2S4_161128_001_C01    L2/3 IT
    ## F2S4_161128_001_D01      L6 IT
    ## F2S4_161128_001_E01      L6 IT
    ## F2S4_161128_001_F01    L2/3 IT
    ## F2S4_161128_001_G01    L2/3 IT
    ## F2S4_161128_001_H01      L6 IT
    ## F2S4_161128_002_A01    L2/3 IT
    ## F2S4_161128_002_B01      L6 IT
    ## F2S4_161128_002_C01      L6 IT
    ## F2S4_161128_002_D01    L2/3 IT
    ## F2S4_161128_002_E01      L6 IT
    ## F2S4_161128_002_F01      L6 IT
    ## F2S4_161128_002_G01      L5 IT
    ## F2S4_161128_002_H01      L6 IT
    ## F2S4_161128_003_A01    L2/3 IT
    ## F2S4_161128_003_B01      L6 IT
    ## F2S4_161128_003_C01      L6 IT
    ## F2S4_161128_003_D01      L6 IT
    ## F2S4_161128_003_E01      L6 IT
    ## F2S4_161128_003_F01      L6 IT
    ## F2S4_161128_003_G01    L2/3 IT
    ## F2S4_161128_003_H01      L6 IT
    ## F2S4_161128_004_A01      L6 IT
    ## F2S4_161128_004_B01      L6 IT
    ## F2S4_161128_004_C01    L2/3 IT
    ## F2S4_161128_004_D01      L6 IT
    ## F2S4_161128_004_E01      L6 IT
    ## F2S4_161128_004_F01    L2/3 IT
    ## F2S4_161128_004_G01      L6 IT
    ## F2S4_161128_004_H01      L6 IT
    ## F2S4_161128_005_A01      L6 IT
    ## F2S4_161128_005_B01      L6 IT
    ## F2S4_161128_005_C01    L2/3 IT
    ## F2S4_161128_005_D01      L5 IT
    ## F2S4_161128_005_E01      L6 IT
    ## F2S4_161128_005_F01      L6 IT
    ## F2S4_161128_005_G01      L6 IT
    ## F2S4_161128_005_H01    L2/3 IT
    ## F2S4_161128_006_A01      L6 IT
    ## F2S4_161128_006_B01      L5 IT
    ## F2S4_161128_006_C01      L5 PT
    ## F2S4_161128_006_D01      L6 IT
    ## F2S4_161128_006_E01      L6 IT
    ## F2S4_161128_006_F01      L6 IT
    ## F2S4_161128_006_G01      L6 IT
    ## F2S4_161128_006_H01      L5 IT
    ## F2S4_161128_007_A01    L2/3 IT
    ## F2S4_161128_007_B01      L6 IT
    ## F2S4_161128_007_C01    L2/3 IT
    ## F2S4_161128_007_D01    L2/3 IT
    ## F2S4_161128_007_E01      L6 IT
    ## F2S4_161128_007_F01      L6 IT
    ## F2S4_161128_007_G01      L6 IT
    ## F2S4_161128_007_H01      L5 IT
    ## F2S4_161128_008_B01      L6 IT
    ## F2S4_161128_008_C01      L6 IT
    ## F2S4_161128_008_E01    L2/3 IT
    ## F2S4_161128_008_F01      L6 IT
    ## F2S4_161128_008_G01      L6 IT
    ## F2S4_161128_008_H01      L6 IT
    ## F2S4_161128_009_A01      L5 IT
    ## F2S4_161128_009_B01      L6 IT
    ## F2S4_161128_009_C01      L6 IT
    ## F2S4_161128_009_D01      L6 IT
    ## F2S4_161128_009_E01      L6 IT
    ## F2S4_161128_009_F01    L2/3 IT
    ## F2S4_161128_009_G01      L6 IT
    ## F2S4_161128_009_H01    L2/3 IT
    ## F2S4_161128_010_A01      L6 IT
    ## F2S4_161128_010_B01      L6 IT
    ## F2S4_161128_010_C01      L6 IT
    ## F2S4_161128_010_D01    L2/3 IT
    ## F2S4_161128_010_E01      L6 IT
    ## F2S4_161128_010_F01      L6 IT
    ## F2S4_161128_010_G01      L6 IT
    ## F2S4_161128_010_H01    L2/3 IT
    ## F2S4_161128_011_A01      L6 IT
    ## F2S4_161128_011_B01      L6 IT
    ## F2S4_161128_011_C01      L6 IT
    ## F2S4_161128_011_D01      L6 IT
    ## F2S4_161128_011_F01    L2/3 IT
    ## F2S4_161128_011_G01    L2/3 IT
    ## F2S4_161128_011_H01      L5 IT
    ## F2S4_161128_012_A01      L6 IT
    ## F2S4_161128_012_B01      L6 IT
    ## F2S4_161128_012_C01      L6 IT
    ## F2S4_161128_012_D01    L2/3 IT
    ## F2S4_161128_012_E01      L6 IT
    ## F2S4_161128_012_F01      L6 IT
    ## F2S4_161128_012_G01      L6 IT
    ## F2S4_161128_012_H01      L6 IT
    ## F2S4_161128_013_A01      L6 IT
    ## F2S4_161128_013_B01      L5 IT
    ## F2S4_161128_013_C01      L6 IT
    ## F2S4_161128_013_D01      L6 IT
    ## F2S4_161128_013_E01      L6 IT
    ## F2S4_161128_013_F01    L2/3 IT
    ## F2S4_161128_013_G01      L6 IT
    ## F2S4_161128_014_A01      L5 IT
    ## F2S4_161128_014_B01      L5 IT
    ## F2S4_161128_014_C01      L6 IT
    ## F2S4_161128_014_D01      L5 IT
    ## F2S4_161128_014_E01      L6 IT
    ## F2S4_161128_014_F01    L2/3 IT
    ## F2S4_161128_014_G01      L5 IT
    ## F2S4_161128_014_H01      L6 IT
    ## F2S4_161128_015_A01      L6 IT
    ## F2S4_161128_015_C01      L5 IT
    ## F2S4_161128_015_D01      L6 IT
    ## F2S4_161128_015_E01      L6 IT
    ## F2S4_161128_015_F01      L6 IT
    ## F2S4_161128_015_G01      L6 IT
    ## F2S4_161128_015_H01    L2/3 IT
    ## F2S4_161128_016_A01      L6 IT
    ## F2S4_161128_016_B01      L6 IT
    ## F2S4_161128_016_C01    L2/3 IT
    ## F2S4_161128_016_D01      L6 IT
    ## F2S4_161128_016_E01    L2/3 IT
    ## F2S4_161128_016_F01      L6 IT
    ## F2S4_161128_016_G01      L6 IT
    ## F2S4_161128_016_H01      L5 IT
    ## F2S4_161128_017_A01      L6 IT
    ## F2S4_161128_017_B01      L6 IT
    ## F2S4_161128_017_D01      L6 IT
    ## F2S4_161128_017_E01      L6 IT
    ## F2S4_161128_017_F01      L6 IT
    ## F2S4_161128_017_G01      L6 IT
    ## F2S4_161128_017_H01      L6 IT
    ## F2S4_161128_018_A01    L2/3 IT
    ## F2S4_161128_018_B01      L6 IT
    ## F2S4_161128_018_D01      L6 IT
    ## F2S4_161128_018_E01      L5 IT
    ## F2S4_161128_018_F01    L2/3 IT
    ## F2S4_161128_018_G01      L5 IT
    ## F2S4_161128_018_H01      L6 IT
    ## F2S4_161129_001_A01    L2/3 IT
    ## F2S4_161129_001_B01    L2/3 IT
    ## F2S4_161129_001_C01      L5 IT
    ## F2S4_161129_001_D01    L2/3 IT
    ## F2S4_161129_001_E01    L2/3 IT
    ## F2S4_161129_001_F01    L2/3 IT
    ## F2S4_161129_001_G01    L2/3 IT
    ## F2S4_161129_001_H01    L2/3 IT
    ## F2S4_161129_002_A01    L2/3 IT
    ## F2S4_161129_002_B01      L6 IT
    ## F2S4_161129_002_C01      L6 IT
    ## F2S4_161129_002_D01      L5 PT
    ## F2S4_161129_002_E01    L2/3 IT
    ## F2S4_161129_002_F01    L2/3 IT
    ## F2S4_161129_002_G01    L2/3 IT
    ## F2S4_161129_002_H01    L2/3 IT
    ## F2S4_161129_003_A01    L2/3 IT
    ## F2S4_161129_003_B01    L2/3 IT
    ## F2S4_161129_003_C01      L5 PT
    ## F2S4_161129_003_D01    L2/3 IT
    ## F2S4_161129_003_F01    L2/3 IT
    ## F2S4_161129_003_G01    L2/3 IT
    ## F2S4_161129_003_H01    L2/3 IT
    ## F2S4_161129_004_A01    L2/3 IT
    ## F2S4_161129_004_B01    L2/3 IT
    ## F2S4_161129_004_C01    L2/3 IT
    ## F2S4_161129_004_D01    L2/3 IT
    ## F2S4_161129_004_E01    L2/3 IT
    ## F2S4_161129_004_F01      L6 IT
    ## F2S4_161129_004_G01    L2/3 IT
    ## F2S4_161129_004_H01    L2/3 IT
    ## F2S4_161129_005_B01    L2/3 IT
    ## F2S4_161129_005_C01    L2/3 IT
    ## F2S4_161129_005_E01    L2/3 IT
    ## F2S4_161129_005_F01    L2/3 IT
    ## F2S4_161129_005_H01    L2/3 IT
    ## F2S4_161129_006_A01    L2/3 IT
    ## F2S4_161129_006_B01    L2/3 IT
    ## F2S4_161129_006_C01    L2/3 IT
    ## F2S4_161129_006_D01    L2/3 IT
    ## F2S4_161129_006_E01      L5 IT
    ## F2S4_161129_006_F01    L2/3 IT
    ## F2S4_161129_006_G01    L2/3 IT
    ## F2S4_161129_006_H01    L2/3 IT
    ## F2S4_161129_007_A01      L6 IT
    ## F2S4_161129_007_B01    L2/3 IT
    ## F2S4_161129_007_C01    L2/3 IT
    ## F2S4_161129_007_D01    L2/3 IT
    ## F2S4_161129_007_E01    L2/3 IT
    ## F2S4_161129_007_F01      L5 PT
    ## F2S4_161129_007_G01    L2/3 IT
    ## F2S4_161129_007_H01    L2/3 IT
    ## F2S4_161129_008_A01    L2/3 IT
    ## F2S4_161129_008_B01    L2/3 IT
    ## F2S4_161129_008_C01    L2/3 IT
    ## F2S4_161129_008_D01    L2/3 IT
    ## F2S4_161129_008_E01    L2/3 IT
    ## F2S4_161129_008_F01    L2/3 IT
    ## F2S4_161129_008_G01    L2/3 IT
    ## F2S4_161129_008_H01    L2/3 IT
    ## F2S4_161129_009_B01    L2/3 IT
    ## F2S4_161206_001_A01      L6 CT
    ## F2S4_161206_001_B01      L6 CT
    ## F2S4_161206_001_C01      L6 CT
    ## F2S4_161206_001_D01      L6 CT
    ## F2S4_161206_001_E01      L6 CT
    ## F2S4_161206_001_F01      L6 CT
    ## F2S4_161206_001_G01      L6 CT
    ## F2S4_161206_001_H01      L6 CT
    ## F2S4_161206_002_A01      L6 CT
    ## F2S4_161206_002_B01      L6 CT
    ## F2S4_161206_002_C01      L6 CT
    ## F2S4_161206_002_D01      L6 CT
    ## F2S4_161206_002_E01      L6 CT
    ## F2S4_161206_002_F01      L6 CT
    ## F2S4_161206_002_G01      L6 CT
    ## F2S4_161206_002_H01        L6b
    ## F2S4_161206_003_A01      L6 CT
    ## F2S4_161206_003_B01      L6 CT
    ## F2S4_161206_003_C01      L6 CT
    ## F2S4_161206_003_D01      L6 CT
    ## F2S4_161206_003_E01      L6 CT
    ## F2S4_161206_003_F01      L6 CT
    ## F2S4_161206_003_G01      L6 CT
    ## F2S4_161206_004_D01      L6 CT
    ## F2S4_161206_004_F01      L6 CT
    ## F2S4_161206_004_G01      L6 CT
    ## F2S4_161206_004_H01      L6 CT
    ## F2S4_161207_001_A01      L6 CT
    ## F2S4_161207_001_B01      L6 CT
    ## F2S4_161207_001_C01      L6 CT
    ## F2S4_161207_001_D01      L6 CT
    ## F2S4_161207_001_E01      L6 CT
    ## F2S4_161207_001_F01      L6 CT
    ## F2S4_161207_001_G01      L6 CT
    ## F2S4_161207_001_H01      L6 CT
    ## F2S4_161207_002_A01      L6 CT
    ## F2S4_161207_002_B01      L6 CT
    ## F2S4_161207_002_C01      L6 CT
    ## F2S4_161207_002_D01      L6 CT
    ## F2S4_161207_002_E01      L6 CT
    ## F2S4_161207_002_F01      L6 CT
    ## F2S4_161207_002_G01      L6 CT
    ## F2S4_161207_002_H01      L6 CT
    ## F2S4_161207_003_B01      L6 CT
    ## F2S4_161207_003_C01      L6 CT
    ## F2S4_161207_003_D01      L6 CT
    ## F2S4_161207_003_E01      L6 CT
    ## F2S4_161207_003_F01      L6 CT
    ## F2S4_161207_003_G01      L6 CT
    ## F2S4_161207_003_H01      L6 CT
    ## F2S4_161207_004_A01      L6 CT
    ## F2S4_161207_004_B01      L6 CT
    ## F2S4_161207_004_C01      L6 CT
    ## F2S4_161207_004_D01      L6 CT
    ## F2S4_161207_004_E01      L6 CT
    ## F2S4_161207_004_F01      L6 CT
    ## F2S4_161207_004_G01      L6 CT
    ## F2S4_161207_005_A01      L6 CT
    ## F2S4_161207_005_B01      L6 CT
    ## F2S4_161207_005_D01      L6 CT
    ## F2S4_161207_005_E01      L6 CT
    ## F2S4_161207_005_F01      L6 CT
    ## F2S4_161207_005_G01      L6 CT
    ## F2S4_161207_005_H01      L6 CT
    ## F2S4_161207_006_A01      L6 CT
    ## F2S4_161207_006_B01      L6 CT
    ## F2S4_161207_006_C01      L6 CT
    ## F2S4_161207_006_D01      L6 CT
    ## F2S4_161207_006_E01      L6 CT
    ## F2S4_161207_006_F01      L6 CT
    ## F2S4_161207_006_G01      L6 CT
    ## F2S4_161207_006_H01      L6 CT
    ## F2S4_161207_007_A01      L6 CT
    ## F2S4_161207_007_B01      L6 CT
    ## F2S4_161207_007_C01      L6 CT
    ## F2S4_161207_007_D01      L6 CT
    ## F2S4_161207_007_E01      L6 CT
    ## F2S4_161207_007_F01      L6 CT
    ## F2S4_161207_007_G01      L6 CT
    ## F2S4_161207_008_A01      L6 CT
    ## F2S4_161207_008_B01      L6 CT
    ## F2S4_161207_008_C01      L6 CT
    ## F2S4_161207_008_D01      L6 CT
    ## F2S4_161207_008_E01      L6 CT
    ## F2S4_161207_008_F01      L6 CT
    ## F2S4_161207_008_G01      L6 CT
    ## F2S4_161207_008_H01      L6 CT
    ## F2S4_161207_009_A01      L6 CT
    ## F2S4_161207_009_B01      L6 CT
    ## F2S4_161207_009_C01      L6 CT
    ## F2S4_161207_009_D01      L6 CT
    ## F2S4_161207_009_E01      L6 CT
    ## F2S4_161207_009_F01      L6 CT
    ## F2S4_161207_009_G01      L6 CT
    ## F2S4_161207_009_H01      L6 CT
    ## F2S4_161207_010_A01      L6 CT
    ## F2S4_161207_010_B01      L6 CT
    ## F2S4_161207_010_C01      L6 CT
    ## F2S4_161207_010_D01      L6 CT
    ## F2S4_161207_010_E01      L6 CT
    ## F2S4_161207_010_F01      L5 PT
    ## F2S4_161207_010_G01      L6 CT
    ## F2S4_161207_011_A01      L6 CT
    ## F2S4_161207_011_B01      L6 CT
    ## F2S4_161207_011_C01      L6 CT
    ## F2S4_161207_011_D01      L6 CT
    ## F2S4_161207_011_E01      L6 CT
    ## F2S4_161207_011_F01      L6 CT
    ## F2S4_161207_011_G01      L6 CT
    ## F2S4_161207_011_H01      L6 CT
    ## F2S4_161207_012_A01      L6 CT
    ## F2S4_161207_012_B01      L6 CT
    ## F2S4_161207_012_C01      L6 CT
    ## F2S4_161207_012_D01      L6 CT
    ## F2S4_161207_012_E01      L6 CT
    ## F2S4_161207_012_F01      L6 CT
    ## F2S4_161207_012_G01      L6 CT
    ## F2S4_161207_012_H01      L6 CT
    ## F2S4_161207_013_B01      L6 CT
    ## F2S4_161207_013_C01      L6 CT
    ## F2S4_161207_013_D01      L6 CT
    ## F2S4_161207_013_E01      L6 CT
    ## F2S4_161207_013_F01      L6 CT
    ## F2S4_161207_013_G01      L6 CT
    ## F2S4_161207_013_H01      L6 CT
    ## F2S4_161207_014_A01      L6 CT
    ## F2S4_161207_014_B01      L6 CT
    ## F2S4_161207_014_C01      L6 CT
    ## F2S4_161207_014_D01      L6 CT
    ## F2S4_161207_014_E01      L6 CT
    ## F2S4_161207_014_F01      L6 CT
    ## F2S4_161207_014_G01      L6 CT
    ## F2S4_161207_014_H01      L6 CT
    ## F2S4_161207_015_A01      L6 CT
    ## F2S4_161207_015_B01      L6 CT
    ## F2S4_161207_015_C01      L6 CT
    ## F2S4_161207_015_E01      L6 CT
    ## F2S4_161207_015_F01      L6 CT
    ## F2S4_161207_015_G01      L6 CT
    ## F2S4_161207_015_H01      L6 CT
    ## F2S4_161207_016_A01      L6 CT
    ## F2S4_161207_016_B01      L6 CT
    ## F2S4_161207_016_C01      L6 CT
    ## F2S4_161207_016_D01      L6 CT
    ## F2S4_161207_016_E01      L6 CT
    ## F2S4_161207_016_G01      L6 CT
    ## F2S4_161207_016_H01      L6 CT
    ## F2S4_161207_017_A01      L6 CT
    ## F2S4_161207_017_B01      L6 CT
    ## F2S4_161207_017_C01      L6 CT
    ## F2S4_161207_017_D01      L6 CT
    ## F2S4_161207_017_E01      L6 CT
    ## F2S4_161207_017_F01      L6 CT
    ## F2S4_161207_017_G01      L6 CT
    ## F2S4_161207_017_H01      L6 CT
    ## F2S4_161207_018_A01      L6 CT
    ## F2S4_161207_018_B01      L6 CT
    ## F2S4_161207_018_C01      L6 CT
    ## F2S4_161207_018_D01      L6 CT
    ## F2S4_161207_018_E01      L6 CT
    ## F2S4_161207_018_F01      L6 CT
    ## F2S4_161207_018_G01      L6 CT
    ## F2S4_161207_018_H01      L6 CT
    ## F2S4_161207_019_A01      L6 CT
    ## F2S4_161207_019_B01      L6 CT
    ## F2S4_161207_019_C01      L6 CT
    ## F2S4_161207_019_D01      L6 CT
    ## F2S4_161207_019_E01      L6 CT
    ## F2S4_161207_019_F01      L6 CT
    ## F2S4_161207_019_G01      L6 CT
    ## F2S4_161207_019_H01      L6 CT
    ## F2S4_161207_020_A01      L6 CT
    ## F2S4_161207_020_B01      L6 CT
    ## F2S4_161207_020_C01      L6 CT
    ## F2S4_161207_020_D01      L6 CT
    ## F2S4_161207_020_E01      L6 CT
    ## F2S4_161207_020_F01      L6 CT
    ## F2S4_161207_020_G01      L6 CT
    ## F2S4_161207_020_H01      L6 CT
    ## F2S4_161207_021_A01      L6 CT
    ## F2S4_161207_021_B01      L6 CT
    ## F2S4_161207_021_C01      L6 CT
    ## F2S4_161207_021_D01      L6 CT
    ## F2S4_161207_021_F01      L6 CT
    ## F2S4_161207_021_G01      L6 CT
    ## F2S4_161207_021_H01      L6 CT
    ## F2S4_161207_022_A01      L6 CT
    ## F2S4_161207_022_B01      L6 CT
    ## F2S4_161207_022_C01      L6 CT
    ## F2S4_161207_022_D01      L6 CT
    ## F2S4_161207_022_E01      L6 CT
    ## F2S4_161207_022_F01      L6 CT
    ## F2S4_161207_022_G01      L6 CT
    ## F2S4_161207_022_H01      L6 CT
    ## F2S4_161207_034_A01      L6 CT
    ## F2S4_161207_034_B01      L6 CT
    ## F2S4_161207_034_C01      L6 CT
    ## F2S4_161207_034_D01      L6 CT
    ## F2S4_161207_034_E01      L6 CT
    ## F2S4_161207_034_F01      L6 CT
    ## F2S4_161207_034_G01      L6 CT
    ## F2S4_161207_034_H01      L6 CT
    ## F2S4_161207_035_A01      L6 CT
    ## F2S4_161207_035_B01      L6 CT
    ## F2S4_161207_035_C01      L6 CT
    ## F2S4_161207_035_D01      L6 CT
    ## F2S4_161207_035_E01      L6 CT
    ## F2S4_170104_001_A01      L6 IT
    ## F2S4_170104_001_B01         L4
    ## F2S4_170104_001_C01      L6 IT
    ## F2S4_170104_001_D01    L2/3 IT
    ## F2S4_170104_001_E01    L2/3 IT
    ## F2S4_170104_001_F01      L6 IT
    ## F2S4_170104_001_G01      L6 IT
    ## F2S4_170104_001_H01      L6 IT
    ## F2S4_170104_002_A01    L2/3 IT
    ## F2S4_170104_002_B01    L2/3 IT
    ## F2S4_170104_002_C01      L6 IT
    ## F2S4_170104_002_D01      L5 IT
    ## F2S4_170104_002_E01      L6 IT
    ## F2S4_170104_002_F01      L6 IT
    ## F2S4_170104_002_G01         L4
    ## F2S4_170104_002_H01    L2/3 IT
    ## F2S4_170104_003_D01      L6 IT
    ## F2S4_170104_003_E01         L4
    ## F2S4_170104_003_F01      L6 IT
    ## F2S4_170104_003_G01    L2/3 IT
    ## F2S4_170104_003_H01      L6 IT
    ## F2S4_170105_001_C01      L5 PT
    ## F2S4_170105_001_D01      L6 IT
    ## F2S4_170105_001_E01      L5 PT
    ## F2S4_170105_001_F01      L6 IT
    ## F2S4_170105_001_G01      L6 IT
    ## F2S4_170105_002_A01      L6 IT
    ## F2S4_170105_002_B01      L5 PT
    ## F2S4_170105_002_C01    L2/3 IT
    ## F2S4_170105_002_D01      L5 PT
    ## F2S4_170105_002_E01      L5 PT
    ## F2S4_170105_002_F01      L6 IT
    ## F2S4_170105_002_G01      L5 PT
    ## F2S4_170105_002_H01    L2/3 IT
    ## F2S4_170105_003_A01      L5 PT
    ## F2S4_170105_003_C01      L6 IT
    ## F2S4_170105_003_D01      L5 IT
    ## F2S4_170105_003_E01      L5 PT
    ## F2S4_170105_003_F01      L5 PT
    ## F2S4_170105_003_G01      L5 PT
    ## F2S4_170105_003_H01      L5 PT
    ## F2S4_170105_004_B01      L6 IT
    ## F2S4_170105_004_C01      L5 PT
    ## F2S4_170105_004_D01      L6 IT
    ## F2S4_170105_004_E01      L5 PT
    ## F2S4_170105_004_F01      L5 PT
    ## F2S4_170105_004_G01      L5 PT
    ## F2S4_170105_005_B01    L2/3 IT
    ## F2S4_170105_005_C01      L5 PT
    ## F2S4_170105_005_D01      L5 PT
    ## F2S4_170105_005_E01      L5 PT
    ## F2S4_170105_005_F01      L5 PT
    ## F2S4_170105_005_H01      L6 IT
    ## F2S4_170105_006_A01      L5 PT
    ## F2S4_170105_006_C01      L6 IT
    ## F2S4_170105_006_F01      L5 PT
    ## F2S4_170105_006_G01      L5 PT
    ## F2S4_170105_006_H01      L5 PT
    ## F2S4_170105_007_A01      L5 PT
    ## F2S4_170105_007_D01      L6 IT
    ## F2S4_170105_007_G01      L5 PT
    ## F2S4_170105_007_H01      L5 PT
    ## F2S4_170105_008_A01      L5 PT
    ## F2S4_170105_008_B01      L5 PT
    ## F2S4_170105_008_C01      L5 PT
    ## F2S4_170105_008_E01      L6 IT
    ## F2S4_170105_008_F01    L2/3 IT
    ## F2S4_170105_009_A01      L5 PT
    ## F2S4_170105_009_B01      L5 PT
    ## F2S4_170105_009_C01      L5 PT
    ## F2S4_170105_009_D01    L2/3 IT
    ## F2S4_170105_009_E01      L6 IT
    ## F2S4_170105_009_F01      L5 PT
    ## F2S4_170105_009_H01      L6 IT
    ## F2S4_170105_010_C01      L5 PT
    ## F2S4_170105_010_D01      L5 PT
    ## F2S4_170105_010_F01      L6 IT
    ## F2S4_170105_010_H01      L6 IT
    ## F2S4_170105_011_B01      L5 PT
    ## F2S4_170105_011_C01      L5 PT
    ## F2S4_170105_011_D01      L6 IT
    ## F2S4_170105_011_F01      L6 IT
    ## F2S4_170105_011_G01      L5 PT
    ## F2S4_170105_011_H01      L5 PT
    ## F2S4_170105_012_A01      L6 IT
    ## F2S4_170105_012_B01      L5 PT
    ## F2S4_170105_012_C01      L5 PT
    ## F2S4_170105_012_E01    L2/3 IT
    ## F2S4_170105_012_F01    L2/3 IT
    ## F2S4_170105_012_G01      L5 PT
    ## F2S4_170105_013_A01      L6 IT
    ## F2S4_170105_013_B01      L5 PT
    ## F2S4_170105_013_D01      L5 PT
    ## F2S4_170116_001_A01         L4
    ## F2S4_170116_001_B01         L4
    ## F2S4_170116_001_C01         L4
    ## F2S4_170116_001_D01         L4
    ## F2S4_170116_001_E01      L5 IT
    ## F2S4_170116_001_F01      L5 IT
    ## F2S4_170116_001_G01         L4
    ## F2S4_170116_001_H01         L4
    ## F2S4_170116_002_A01         L4
    ## F2S4_170116_002_B01      L5 IT
    ## F2S4_170116_002_C01         L4
    ## F2S4_170116_002_D01         L4
    ## F2S4_170116_002_E01         L4
    ## F2S4_170116_002_F01         L4
    ## F2S4_170116_002_G01         L4
    ## F2S4_170116_002_H01         L4
    ## F2S4_170116_003_A01         L4
    ## F2S4_170116_003_B01         L4
    ## F2S4_170116_003_D01         L4
    ## F2S4_170116_003_E01         L4
    ## F2S4_170116_003_F01         L4
    ## F2S4_170116_003_G01         L4
    ## F2S4_170116_003_H01         L4
    ## F2S4_170116_004_A01         L4
    ## F2S4_170116_004_B01         L4
    ## F2S4_170116_004_D01         L4
    ## F2S4_170116_004_E01      L5 IT
    ## F2S4_170116_004_F01         L4
    ## F2S4_170116_004_G01         L4
    ## F2S4_170116_005_A01         L4
    ## F2S4_170116_005_B01         L4
    ## F2S4_170116_005_C01      L5 IT
    ## F2S4_170116_005_D01         L4
    ## F2S4_170116_005_E01         L4
    ## F2S4_170116_005_F01         L4
    ## F2S4_170116_005_G01         L4
    ## F2S4_170116_005_H01         L4
    ## F2S4_170116_006_A01         L4
    ## F2S4_170116_006_B01         L4
    ## F2S4_170116_006_C01         L4
    ## F2S4_170116_006_D01         L4
    ## F2S4_170116_006_E01         L4
    ## F2S4_170116_006_F01         L4
    ## F2S4_170116_006_G01         L4
    ## F2S4_170116_006_H01         L4
    ## F2S4_170116_007_A01         L4
    ## F2S4_170116_007_B01         L4
    ## F2S4_170116_007_C01         L4
    ## F2S4_170116_007_D01         L4
    ## F2S4_170116_007_E01         L4
    ## F2S4_170116_007_G01         L4
    ## F2S4_170116_007_H01         L4
    ## F2S4_170116_008_A01      L5 IT
    ## F2S4_170116_008_B01      L5 IT
    ## F2S4_170116_008_C01      L5 IT
    ## F2S4_170116_008_D01      L5 IT
    ## F2S4_170116_008_E01      L5 IT
    ## F2S4_170116_008_H01         L4
    ## F2S4_170116_009_A01         L4
    ## F2S4_170116_009_B01         L4
    ## F2S4_170116_009_C01         L4
    ## F2S4_170116_009_D01         L4
    ## F2S4_170116_009_E01         L4
    ## F2S4_170116_009_F01         L4
    ## F2S4_170116_009_G01         L4
    ## F2S4_170116_009_H01         L4
    ## F2S4_170116_010_A01         L4
    ## F2S4_170116_010_B01         L4
    ## F2S4_170116_010_C01         L4
    ## F2S4_170116_010_D01         L4
    ## F2S4_170116_010_E01         L4
    ## F2S4_170116_010_F01         L4
    ## F2S4_170116_010_G01         L4
    ## F2S4_170116_010_H01         L4
    ## F2S4_170116_011_A01         L4
    ## F2S4_170116_011_B01         L4
    ## F2S4_170116_011_C01         L4
    ## F2S4_170116_011_D01         L4
    ## F2S4_170116_011_E01         L4
    ## F2S4_170116_011_F01         L4
    ## F2S4_170116_011_G01         L4
    ## F2S4_170116_011_H01         L4
    ## F2S4_170116_012_A01      L5 IT
    ## F2S4_170116_012_B01         L4
    ## F2S4_170116_012_C01         L4
    ## F2S4_170116_012_D01         L4
    ## F2S4_170116_012_E01         L4
    ## F2S4_170116_012_F01         L4
    ## F2S4_170116_012_G01         L4
    ## F2S4_170116_012_H01         L4
    ## F2S4_170116_013_A01         L4
    ## F2S4_170116_013_B01         L4
    ## F2S4_170116_013_C01         L4
    ## F2S4_170116_013_D01         L4
    ## F2S4_170116_013_E01         L4
    ## F2S4_170116_013_F01         L4
    ## F2S4_170116_013_G01         L4
    ## F2S4_170116_013_H01      L5 IT
    ## F2S4_170116_014_A01         L4
    ## F2S4_170116_014_B01         L4
    ## F2S4_170116_014_C01         L4
    ## F2S4_170116_014_D01         L4
    ## F2S4_170116_014_E01         L4
    ## F2S4_170116_014_F01         L4
    ## F2S4_170116_014_G01         L4
    ## F2S4_170116_014_H01         L4
    ## F2S4_170116_015_A01      Astro
    ## F2S4_170116_015_B01      L5 IT
    ## F2S4_170116_015_C01      L5 IT
    ## F2S4_170116_015_D01      L5 IT
    ## F2S4_170116_015_E01      L5 IT
    ## F2S4_170116_015_F01      L5 IT
    ## F2S4_170116_015_G01      L5 IT
    ## F2S4_170116_015_H01      L5 IT
    ## F2S4_170116_016_H01      L5 IT
    ## F2S4_170123_001_A01      L5 IT
    ## F2S4_170123_001_B01      L5 IT
    ## F2S4_170123_001_C01    L2/3 IT
    ## F2S4_170123_001_D01    L2/3 IT
    ## F2S4_170123_001_E01      L5 IT
    ## F2S4_170123_001_G01      L5 IT
    ## F2S4_170123_001_H01    L2/3 IT
    ## F2S4_170123_002_A01      L6 IT
    ## F2S4_170123_002_B01      L6 IT
    ## F2S4_170123_002_C01    L2/3 IT
    ## F2S4_170124_001_A01    L2/3 IT
    ## F2S4_170124_001_B01      L6 IT
    ## F2S4_170124_001_C01      L5 IT
    ## F2S4_170124_001_D01    L2/3 IT
    ## F2S4_170124_001_E01      L5 IT
    ## F2S4_170124_001_F01      L5 IT
    ## F2S4_170124_001_G01      L6 IT
    ## F2S4_170124_001_H01    L2/3 IT
    ## F2S4_170124_002_A01      L6 IT
    ## F2S4_170124_002_B01      L6 IT
    ## F2S4_170124_002_C01      L6 IT
    ## F2S4_170124_002_D01      L6 IT
    ## F2S4_170124_002_F01      L6 IT
    ## F2S4_170124_002_G01    L2/3 IT
    ## F2S4_170124_002_H01      L6 IT
    ## F2S4_170125_001_A01      L5 PT
    ## F2S4_170125_001_B01      L5 PT
    ## F2S4_170125_001_C01      L6 CT
    ## F2S4_170125_001_D01      Oligo
    ## F2S4_170125_001_E01      L6 CT
    ## F2S4_170125_001_F01        L6b
    ## F2S4_170125_001_G01      L6 CT
    ## F2S4_170125_001_H01         L4
    ## F2S4_170125_002_A01    L2/3 IT
    ## F2S4_170125_002_B01      L6 CT
    ## F2S4_170125_002_C01      L6 IT
    ## F2S4_170125_002_E01        L6b
    ## F2S4_170125_002_F01    L2/3 IT
    ## F2S4_170125_002_G01      L6 CT
    ## F2S4_170125_002_H01      L6 CT
    ## F2S4_170126_001_B01         L4
    ## F2S4_170126_001_C01      L5 IT
    ## F2S4_170126_001_D01      L6 IT
    ## F2S4_170126_001_G01         L4
    ## F2S4_170126_001_H01      L5 IT
    ## F2S4_170126_002_A01      L6 IT
    ## F2S4_170126_002_B01         L4
    ## F2S4_170126_002_E01      L6 IT
    ## F2S4_170126_002_G01         L4
    ## F2S4_170126_002_H01      L5 IT
    ## F2S4_170126_003_A01      L5 IT
    ## F2S4_170126_003_B01      L6 IT
    ## F2S4_170126_003_C01      L6 IT
    ## F2S4_170126_003_D01      L6 IT
    ## F2S4_170126_003_E01      L6 IT
    ## F2S4_170126_003_F01         L4
    ## F2S4_170126_004_A01         L4
    ## F2S4_170126_004_C01      L5 IT
    ## F2S4_170126_004_D01         L4
    ## F2S4_170126_004_E01      L6 IT
    ## F2S4_170126_004_F01    L2/3 IT
    ## F2S4_170126_004_G01      L6 IT
    ## F2S4_170126_004_H01         L4
    ## F2S4_170126_005_A01      L5 IT
    ## F2S4_170126_005_B01      L6 IT
    ## F2S4_170126_005_C01      L6 IT
    ## F2S4_170126_005_D01      L5 IT
    ## F2S4_170126_005_E01      L6 IT
    ## F2S4_170126_005_G01    L2/3 IT
    ## F2S4_170126_005_H01      L5 IT
    ## F2S4_170126_006_A01      L5 IT
    ## F2S4_170126_006_D01         L4
    ## F2S4_170126_006_E01      L6 IT
    ## F2S4_170126_006_F01      Lamp5
    ## F2S4_170126_006_G01      L6 IT
    ## F2S4_170126_006_H01         L4
    ## F2S4_170126_007_A01      L6 CT
    ## F2S4_170126_007_B01      L5 IT
    ## F2S4_170126_007_C01         L4
    ## F2S4_170126_007_D01         L4
    ## F2S4_170126_007_E01         L4
    ## F2S4_170126_007_F01         L4
    ## F2S4_170126_007_G01         L4
    ## F2S4_170126_007_H01         L4
    ## F2S4_170126_008_A01         L4
    ## F2S4_170126_008_B01         L4
    ## F2S4_170126_008_C01      L5 IT
    ## F2S4_170126_008_D01         L4
    ## F2S4_170126_008_G01         L4
    ## F2S4_170126_008_H01         L4
    ## F2S4_170126_009_A01         L4
    ## F2S4_170126_009_B01         L4
    ## F2S4_170126_009_C01         L4
    ## F2S4_170126_009_E01         L4
    ## F2S4_170126_009_F01         L4
    ## F2S4_170126_009_G01         L4
    ## F2S4_170126_010_C01         L4
    ## F2S4_170126_010_D01         L4
    ## F2S4_170126_010_F01         L4
    ## F2S4_170126_010_G01         L4
    ## F2S4_170126_010_H01         L4
    ## F2S4_170126_011_B01      L5 IT
    ## F2S4_170126_011_C01         L4
    ## F2S4_170126_011_D01         L4
    ## F2S4_170126_011_E01         L4
    ## F2S4_170126_011_F01         L4
    ## F2S4_170126_011_G01         L4
    ## F2S4_170126_011_H01      L5 IT
    ## F2S4_170126_012_A01         L4
    ## F2S4_170126_012_D01         L4
    ## F2S4_170126_012_E01         L4
    ## F2S4_170126_012_F01         L4
    ## F2S4_170126_012_G01         L4
    ## F2S4_170126_012_H01         L4
    ## F2S4_170126_013_B01         L4
    ## F2S4_170126_013_C01         L4
    ## F2S4_170126_013_D01      L5 IT
    ## F2S4_170126_013_E01         L4
    ## F2S4_170126_013_F01         L4
    ## F2S4_170126_013_G01         L4
    ## F2S4_170126_013_H01      L5 IT
    ## F2S4_170126_014_A01         L4
    ## F2S4_170126_014_D01         L4
    ## F2S4_170126_014_E01         L4
    ## F2S4_170126_014_F01         L4
    ## F2S4_170126_014_H01         L4
    ## F2S4_170126_015_B01      L5 IT
    ## F2S4_170126_015_C01         L4
    ## F2S4_170126_015_D01         L4
    ## F2S4_170126_015_E01         L4
    ## F2S4_170126_015_F01         L4
    ## F2S4_170126_015_G01         L4
    ## F2S4_170126_015_H01         L4
    ## F2S4_170126_016_A01         L4
    ## F2S4_170126_016_B01         L4
    ## F2S4_170126_016_C01         L4
    ## F2S4_170126_016_D01         L4
    ## F2S4_170126_016_E01         L4
    ## F2S4_170126_016_F01         L4
    ## F2S4_170126_016_G01         L4
    ## F2S4_170126_016_H01         L4
    ## F2S4_170126_017_A01         L4
    ## F2S4_170126_017_B01         L4
    ## F2S4_170126_017_C01         L4
    ## F2S4_170126_017_D01         L4
    ## F2S4_170126_017_F01         L4
    ## F2S4_170126_017_G01         L4
    ## F2S4_170126_017_H01         L4
    ## F2S4_170126_018_A01         L4
    ## F2S4_170126_018_B01         L4
    ## F2S4_170126_018_C01      L5 IT
    ## F2S4_170126_018_D01         L4
    ## F2S4_170126_018_F01         L4
    ## F2S4_170126_019_A01      L5 IT
    ## F2S4_170126_019_B01         L4
    ## F2S4_170126_019_C01      L5 IT
    ## F2S4_170126_019_D01         L4
    ## F2S4_170126_019_E01         L4
    ## F2S4_170126_019_F01         L4
    ## F2S4_170126_019_G01         L4
    ## F2S4_170126_019_H01         L4
    ## F2S4_170126_020_A01         L4
    ## F2S4_170126_020_B01         L4
    ## F2S4_170126_020_D01         L4
    ## F2S4_170126_020_E01         L4
    ## F2S4_170126_020_F01         L4
    ## F2S4_170126_020_G01         L4
    ## F2S4_170126_020_H01         L4
    ## F2S4_170126_021_A01         L4
    ## F2S4_170126_021_B01         L4
    ## F2S4_170126_021_C01         L4
    ## F2S4_170126_021_D01         L4
    ## F2S4_170126_021_E01         L4
    ## F2S4_170126_021_F01         L4
    ## F2S4_170126_021_G01         L4
    ## F2S4_170126_021_H01         L4
    ## F2S4_170126_022_A01         L4
    ## F2S4_170126_022_B01         L4
    ## F2S4_170126_022_C01         L4
    ## F2S4_170126_022_D01         L4
    ## F2S4_170126_022_E01         L4
    ## F2S4_170126_022_F01         L4
    ## F2S4_170126_022_G01         L4
    ## F2S4_170126_022_H01         L4
    ## F2S4_170126_023_A01      Pvalb
    ## F2S4_170126_023_B01        Sst
    ## F2S4_170126_023_C01      Pvalb
    ## F2S4_170126_023_D01      Pvalb
    ## F2S4_170126_023_E01        Sst
    ## F2S4_170126_023_F01      Pvalb
    ## F2S4_170126_023_G01      L5 IT
    ## F2S4_170126_024_G01        Sst
    ## F2S4_170126_024_H01        Sst
    ## F2S4_170130_001_B01      L5 IT
    ## F2S4_170130_001_C01      L6 IT
    ## F2S4_170130_001_F01      L6 IT
    ## F2S4_170130_001_G01    L2/3 IT
    ## F2S4_170130_001_H01      L6 IT
    ## F2S4_170130_002_C01      L6 IT
    ## F2S4_170130_002_D01      L6 IT
    ## F2S4_170130_002_E01      L5 PT
    ## F2S4_170130_002_G01    L2/3 IT
    ## F2S4_170130_002_H01      L6 IT
    ## F2S4_170130_003_A01    L2/3 IT
    ## F2S4_170130_003_D01      L6 IT
    ## F2S4_170130_003_E01    L2/3 IT
    ## F2S4_170130_003_F01      L6 IT
    ## F2S4_170130_003_H01      L6 IT
    ## F2S4_170130_004_A01      L6 IT
    ## F2S4_170130_004_B01      L6 IT
    ## F2S4_170130_004_C01      L6 IT
    ## F2S4_170130_004_D01      L5 PT
    ## F2S4_170130_004_E01         L4
    ## F2S4_170130_004_F01      L6 IT
    ## F2S4_170130_004_G01    L2/3 IT
    ## F2S4_170130_004_H01      L6 IT
    ## F2S4_170130_005_A01      L6 IT
    ## F2S4_170130_005_B01      L6 IT
    ## F2S4_170130_005_C01      L6 IT
    ## F2S4_170130_005_D01      L6 IT
    ## F2S4_170130_005_E01      L6 IT
    ## F2S4_170130_005_F01    L2/3 IT
    ## F2S4_170130_005_G01      L6 IT
    ## F2S4_170130_005_H01      L6 IT
    ## F2S4_170130_006_E01    L2/3 IT
    ## F2S4_170130_006_F01      L5 IT
    ## F2S4_170130_006_G01    L2/3 IT
    ## F2S4_170130_006_H01      L6 IT
    ## F2S4_170130_007_A01      L6 IT
    ## F2S4_170130_007_B01      L6 IT
    ## F2S4_170130_007_C01      L5 PT
    ## F2S4_170130_007_D01      Oligo
    ## F2S4_170130_008_C01      L6 IT
    ## F2S4_170130_008_D01    L2/3 IT
    ## F2S4_170130_008_E01    L2/3 IT
    ## F2S4_170130_009_A01      L6 IT
    ## F2S4_170130_010_B01      L6 IT
    ## F2S4_170130_010_D01      L5 PT
    ## F2S4_170130_010_H01    L2/3 IT
    ## F2S4_170131_001_A01      L6 IT
    ## F2S4_170131_001_B01      L6 IT
    ## F2S4_170131_001_C01      L6 IT
    ## F2S4_170201_001_A01      L5 PT
    ## F2S4_170201_001_C01      L5 PT
    ## F2S4_170201_001_D01      L5 PT
    ## F2S4_170201_001_H01      L5 PT
    ## F2S4_170201_002_A01      L6 CT
    ## F2S4_170201_002_E01      L6 CT
    ## F2S4_170201_002_H01      L5 PT
    ## F2S4_170201_003_A01      L5 PT
    ## F2S4_170201_003_C01      L5 PT
    ## F2S4_170201_003_E01        L6b
    ## F2S4_170201_003_G01      L5 PT
    ## F2S4_170201_003_H01      L5 PT
    ## F2S4_170201_004_A01      L6 CT
    ## F2S4_170201_004_B01      L5 PT
    ## F2S4_170201_004_C01      L6 IT
    ## F2S4_170201_004_D01      L5 PT
    ## F2S4_170201_004_E01        L6b
    ## F2S4_170201_004_F01      L6 CT
    ## F2S4_170201_004_G01      L6 CT
    ## F2S4_170201_004_H01      L6 CT
    ## F2S4_170201_005_A01      L5 PT
    ## F2S4_170201_005_E01      L6 CT
    ## F2S4_170201_005_F01      L5 PT
    ## F2S4_170201_006_B01      L6 CT
    ## F2S4_170201_006_C01      L6 CT
    ## F2S4_170201_006_D01      L5 PT
    ## F2S4_170201_006_E01      L6 CT
    ## F2S4_170201_006_F01        L6b
    ## F2S4_170201_007_C01      L5 PT
    ## F2S4_170201_007_F01      L6 CT
    ## F2S4_170201_007_G01      L5 PT
    ## F2S4_170201_008_E01      L5 PT
    ## F2S4_170201_008_F01      L5 PT
    ## F2S4_170201_009_A01      L6 CT
    ## F2S4_170201_009_C01        L6b
    ## F2S4_170201_009_D01      L5 PT
    ## F2S4_170201_009_E01      L6 CT
    ## F2S4_170201_009_F01      L5 PT
    ## F2S4_170201_010_D01      L5 PT
    ## F2S4_170201_010_F01      L5 PT
    ## F2S4_170201_011_A01      L5 PT
    ## F2S4_170201_011_D01        L6b
    ## F2S4_170201_011_E01    L2/3 IT
    ## F2S4_170201_012_B01      L5 PT
    ## F2S4_170201_012_C01      L5 PT
    ## F2S4_170201_012_E01      L5 PT
    ## F2S4_170201_012_F01      L5 PT
    ## F2S4_170201_012_G01      L6 CT
    ## F2S4_170201_013_A01      L6 CT
    ## F2S4_170201_013_B01      L5 PT
    ## F2S4_170201_013_C01      L5 PT
    ## F2S4_170201_013_H01      L5 PT
    ## F2S4_170201_014_A01      L6 CT
    ## F2S4_170201_014_B01      L5 PT
    ## F2S4_170201_014_C01      L6 CT
    ## F2S4_170201_014_D01      L5 PT
    ## F2S4_170201_014_F01      L6 CT
    ## F2S4_170201_014_G01      Oligo
    ## F2S4_170201_014_H01    L2/3 IT
    ## F2S4_170201_015_A01      L6 CT
    ## F2S4_170201_015_B01        L6b
    ## F2S4_170201_015_C01      L5 PT
    ## F2S4_170201_015_D01      L5 PT
    ## F2S4_170201_015_E01      L6 CT
    ## F2S4_170201_016_E01      L6 CT
    ## F2S4_170201_016_F01      L5 PT
    ## F2S4_170201_016_G01      L5 PT
    ## F2S4_170201_017_C01      L5 PT
    ## F2S4_170201_017_D01      L6 CT
    ## F2S4_170201_017_F01      L5 PT
    ## F2S4_170201_017_G01      L6 CT
    ## F2S4_170201_018_A01      L6 CT
    ## F2S4_170201_018_B01      L6 CT
    ## F2S4_170201_018_D01      L6 CT
    ## F2S4_170201_018_E01      L5 PT
    ## F2S4_170201_018_F01      L6 CT
    ## F2S4_170201_018_H01    L2/3 IT
    ## F2S4_170201_019_A01      Oligo
    ## F2S4_170201_019_C01      L6 CT
    ## F2S4_170202_001_A01      L5 PT
    ## F2S4_170202_001_G01      L5 PT
    ## F2S4_170202_001_H01    L2/3 IT
    ## F2S4_170202_002_F01      L5 PT
    ## F2S4_170202_002_G01      L5 PT
    ## F2S4_170202_003_B01      L5 PT
    ## F2S4_170202_003_E01      L5 PT
    ## F2S4_170202_004_A01      L5 PT
    ## F2S4_170202_004_D01      L5 PT
    ## F2S4_170202_004_E01      L5 PT
    ## F2S4_170202_005_A01      L5 PT
    ## F2S4_170202_005_B01      L5 PT
    ## F2S4_170202_005_G01      L5 PT
    ## F2S4_170202_006_A01      L6 IT
    ## F2S4_170202_006_B01      L6 IT
    ## F2S4_170202_006_C01      L6 IT
    ## F2S4_170202_006_D01      L6 IT
    ## F2S4_170202_006_E01      L5 IT
    ## F2S4_170202_006_F01    L2/3 IT
    ## F2S4_170202_007_A01      L6 IT
    ## F2S4_170202_007_B01         L4
    ## F2S4_170202_007_C01      L6 IT
    ## F2S4_170202_007_D01      L6 IT
    ## F2S4_170202_007_E01      L6 IT
    ## F2S4_170202_007_F01      L6 IT
    ## F2S4_170202_007_G01      L5 IT
    ## F2S4_170202_007_H01      L6 IT
    ## F2S4_170202_008_A01      L6 IT
    ## F2S4_170202_008_B01      L6 IT
    ## F2S4_170202_008_C01      L6 IT
    ## F2S4_170202_008_D01      L5 IT
    ## F2S4_170202_008_E01      L6 IT
    ## F2S4_170202_008_F01      L6 IT
    ## F2S4_170202_008_G01      L6 IT
    ## F2S4_170202_008_H01      L6 IT
    ## F2S4_170202_009_A01      L6 IT
    ## F2S4_170202_009_B01      L5 IT
    ## F2S4_170202_009_C01      L6 IT
    ## F2S4_170202_009_D01      L6 IT
    ## F2S4_170202_009_E01      L6 IT
    ## F2S4_170202_009_F01      L6 IT
    ## F2S4_170202_009_G01      L6 IT
    ## F2S4_170202_009_H01      L6 IT
    ## F2S4_170202_010_A01      L6 IT
    ## F2S4_170202_010_B01      L5 IT
    ## F2S4_170202_010_C01      L6 IT
    ## F2S4_170202_010_D01      L6 IT
    ## F2S4_170202_010_E01      L6 IT
    ## F2S4_170202_010_F01      L6 IT
    ## F2S4_170202_010_G01      L6 IT
    ## F2S4_170202_010_H01      L6 IT
    ## F2S4_170202_011_A01      L6 IT
    ## F2S4_170202_011_B01      L6 IT
    ## F2S4_170202_011_C01      L6 IT
    ## F2S4_170202_011_D01      L5 IT
    ## F2S4_170202_011_E01      L6 IT
    ## F2S4_170202_011_F01      L5 IT
    ## F2S4_170202_011_G01      L6 IT
    ## F2S4_170202_011_H01      L6 IT
    ## F2S4_170202_012_A01      L6 IT
    ## F2S4_170202_012_B01      L6 IT
    ## F2S4_170202_012_C01      L6 IT
    ## F2S4_170202_012_D01      L6 IT
    ## F2S4_170202_012_E01      L6 IT
    ## F2S4_170202_012_F01      L6 IT
    ## F2S4_170202_012_G01      L6 IT
    ## F2S4_170202_012_H01      L6 IT
    ## F2S4_170202_013_C01      L5 IT
    ## F2S4_170202_013_D01      L6 IT
    ## F2S4_170202_013_E01      L6 IT
    ## F2S4_170202_013_F01      L6 IT
    ## F2S4_170202_013_G01      L6 IT
    ## F2S4_170202_013_H01      L6 IT
    ## F2S4_170203_001_A01      L5 PT
    ## F2S4_170203_001_B01      L5 PT
    ## F2S4_170203_001_E01      L5 PT
    ## F2S4_170203_001_G01      L5 PT
    ## F2S4_170203_002_E01      L5 PT
    ## F2S4_170203_003_A01      L5 PT
    ## F2S4_170203_004_A01      L6 IT
    ## F2S4_170203_004_B01      L6 IT
    ## F2S4_170203_004_C01      L6 IT
    ## F2S4_170203_004_D01      L6 IT
    ## F2S4_170203_004_E01      L6 IT
    ## F2S4_170203_004_F01      L6 IT
    ## F2S4_170203_004_G01      L6 IT
    ## F2S4_170203_004_H01      L6 IT
    ## F2S4_170203_005_A01      L5 IT
    ## F2S4_170203_005_C01      L6 IT
    ## F2S4_170203_005_D01      L6 IT
    ## F2S4_170203_005_E01      L5 IT
    ## F2S4_170203_005_F01      L6 IT
    ## F2S4_170203_005_G01      L6 IT
    ## F2S4_170203_005_H01      L6 IT
    ## F2S4_170203_006_A01      L6 IT
    ## F2S4_170203_006_B01      L6 IT
    ## F2S4_170203_006_C01      L6 IT
    ## F2S4_170203_006_E01      L6 IT
    ## F2S4_170203_006_F01         L4
    ## F2S4_170203_006_G01      L6 IT
    ## F2S4_170203_006_H01      L6 IT
    ## F2S4_170203_007_A01      L6 IT
    ## F2S4_170203_007_B01      L6 IT
    ## F2S4_170203_007_C01         L4
    ## F2S4_170203_007_D01    L2/3 IT
    ## F2S4_170203_007_E01         L4
    ## F2S4_170203_007_F01    L2/3 IT
    ## F2S4_170203_007_G01      L6 IT
    ## F2S4_170203_007_H01      L6 IT
    ## F2S4_170203_008_A01      L6 IT
    ## F2S4_170203_008_B01      L6 IT
    ## F2S4_170203_008_C01      L6 IT
    ## F2S4_170203_008_D01      L6 IT
    ## F2S4_170203_008_E01      L6 IT
    ## F2S4_170203_008_F01      L6 IT
    ## F2S4_170203_008_G01      L6 IT
    ## F2S4_170203_008_H01         L4
    ## F2S4_170203_009_A01      L5 IT
    ## F2S4_170203_009_B01    L2/3 IT
    ## F2S4_170203_009_C01      L6 IT
    ## F2S4_170203_009_D01         L4
    ## F2S4_170203_009_E01      L6 IT
    ## F2S4_170203_009_F01      L5 IT
    ## F2S4_170203_009_G01      L6 IT
    ## F2S4_170203_009_H01      L6 IT
    ## F2S4_170203_010_A01         L4
    ## F2S4_170203_010_B01      L6 IT
    ## F2S4_170203_010_C01      L6 IT
    ## F2S4_170203_010_D01         L4
    ## F2S4_170203_010_E01      L6 IT
    ## F2S4_170203_010_F01      L6 IT
    ## F2S4_170203_010_G01    L2/3 IT
    ## F2S4_170203_010_H01      L6 IT
    ## F2S4_170203_011_A01         L4
    ## F2S4_170203_011_B01         L4
    ## F2S4_170203_011_C01    L2/3 IT
    ## F2S4_170203_011_D01      L6 IT
    ## F2S4_170203_011_E01    L2/3 IT
    ## F2S4_170203_011_F01      L6 IT
    ## F2S4_170203_011_G01      L5 IT
    ## F2S4_170203_011_H01         L4
    ## F2S4_170203_012_A01      L6 IT
    ## F2S4_170203_012_B01      L6 IT
    ## F2S4_170203_012_C01         L4
    ## F2S4_170203_012_D01         L4
    ## F2S4_170203_012_E01      L6 IT
    ## F2S4_170203_012_F01         L4
    ## F2S4_170203_012_G01      L6 IT
    ## F2S4_170203_012_H01      L6 IT
    ## F2S4_170203_013_A01    L2/3 IT
    ## F2S4_170203_013_B01      L5 IT
    ## F2S4_170203_013_C01      L6 IT
    ## F2S4_170203_013_D01      L6 IT
    ## F2S4_170203_013_E01      L6 IT
    ## F2S4_170203_013_F01         L4
    ## F2S4_170203_013_G01      L6 IT
    ## F2S4_170203_013_H01      L6 IT
    ## F2S4_170203_014_A01      L6 IT
    ## F2S4_170203_014_B01      L6 IT
    ## F2S4_170203_014_C01      L6 IT
    ## F2S4_170203_014_D01    L2/3 IT
    ## F2S4_170203_014_E01      L6 IT
    ## F2S4_170203_014_F01         L4
    ## F2S4_170203_014_G01      L6 IT
    ## F2S4_170203_014_H01      L6 IT
    ## F2S4_170203_015_A01      L6 IT
    ## F2S4_170203_015_B01         L4
    ## F2S4_170203_015_C01      L6 IT
    ## F2S4_170203_015_D01      L6 IT
    ## F2S4_170203_015_F01      L6 IT
    ## F2S4_170203_015_G01      L6 IT
    ## F2S4_170203_015_H01      L6 IT
    ## F2S4_170203_016_A01         L4
    ## F2S4_170203_016_B01      L5 IT
    ## F2S4_170203_016_C01      L6 IT
    ## F2S4_170203_016_D01      L6 IT
    ## F2S4_170203_016_E01      L6 IT
    ## F2S4_170203_016_F01      L6 IT
    ## F2S4_170203_016_G01      L6 IT
    ## F2S4_170203_016_H01    L2/3 IT
    ## F2S4_170203_017_A01      L6 IT
    ## F2S4_170203_017_B01      L6 IT
    ## F2S4_170203_017_C01         L4
    ## F2S4_170203_017_D01      L6 IT
    ## F2S4_170203_017_E01      L6 IT
    ## F2S4_170203_017_F01      L6 IT
    ## F2S4_170203_017_G01      L6 IT
    ## F2S4_170203_017_H01      L6 CT
    ## F2S4_170203_018_A01      L6 IT
    ## F2S4_170203_018_B01      L6 IT
    ## F2S4_170203_018_C01      L6 IT
    ## F2S4_170203_018_D01      L5 IT
    ## F2S4_170203_018_E01      L6 IT
    ## F2S4_170203_018_F01      L6 IT
    ## F2S4_170203_018_G01      L6 IT
    ## F2S4_170203_018_H01      L6 IT
    ## F2S4_170203_019_A01    L2/3 IT
    ## F2S4_170203_019_B01      L6 IT
    ## F2S4_170203_019_C01      L6 IT
    ## F2S4_170203_019_D01         L4
    ## F2S4_170203_019_E01      L6 IT
    ## F2S4_170203_019_F01         L4
    ## F2S4_170203_019_G01      L6 IT
    ## F2S4_170203_019_H01         L4
    ## F2S4_170203_020_A01         L4
    ## F2S4_170203_020_B01      L5 IT
    ## F2S4_170203_020_C01      L5 IT
    ## F2S4_170203_020_D01      L6 IT
    ## F2S4_170203_020_E01      L6 IT
    ## F2S4_170203_020_F01    L2/3 IT
    ## F2S4_170203_020_G01      L6 IT
    ## F2S4_170203_020_H01      L6 IT
    ## F2S4_170203_021_A01      L5 IT
    ## F2S4_170203_021_B01      L6 IT
    ## F2S4_170203_021_C01      L6 IT
    ## F2S4_170203_021_D01      L6 IT
    ## F2S4_170203_021_E01      L5 IT
    ## F2S4_170203_021_F01      L6 IT
    ## F2S4_170203_021_G01      L6 IT
    ## F2S4_170203_021_H01    L2/3 IT
    ## F2S4_170203_022_A01      L6 IT
    ## F2S4_170203_022_B01      L6 IT
    ## F2S4_170203_022_C01      L6 IT
    ## F2S4_170203_022_D01      L6 IT
    ## F2S4_170203_022_E01      L6 IT
    ## F2S4_170203_022_F01      L6 IT
    ## F2S4_170203_022_G01      L6 IT
    ## F2S4_170203_022_H01      L6 IT
    ## F2S4_170203_023_A01      L6 IT
    ## F2S4_170203_023_B01      L6 IT
    ## F2S4_170203_023_C01      L5 IT
    ## F2S4_170203_023_E01      L6 IT
    ## F2S4_170203_023_F01    L2/3 IT
    ## F2S4_170203_023_G01      L6 IT
    ## F2S4_170203_023_H01      L6 IT
    ## F2S4_170203_024_A01      L6 IT
    ## F2S4_170203_024_B01      L6 IT
    ## F2S4_170203_024_C01      L6 IT
    ## F2S4_170203_024_D01      L6 IT
    ## F2S4_170203_024_E01      L6 IT
    ## F2S4_170203_024_F01      L6 IT
    ## F2S4_170203_024_G01      L6 IT
    ## F2S4_170203_024_H01      L6 IT
    ## F2S4_170203_025_A01      L6 IT
    ## F2S4_170203_025_B01      L6 IT
    ## F2S4_170203_025_C01      L6 IT
    ## F2S4_170203_025_D01         L4
    ## F2S4_170203_025_E01      L6 IT
    ## F2S4_170203_025_F01      L6 IT
    ## F2S4_170203_025_G01      L6 IT
    ## F2S4_170203_025_H01      L6 IT
    ## F2S4_170203_026_A01      L6 IT
    ## F2S4_170203_026_B01      L6 IT
    ## F2S4_170203_026_C01      L6 IT
    ## F2S4_170203_026_D01      L5 IT
    ## F2S4_170203_026_E01      L6 IT
    ## F2S4_170203_026_F01         L4
    ## F2S4_170203_026_G01      L6 IT
    ## F2S4_170203_026_H01      L6 IT
    ## F2S4_170320_001_A01      Astro
    ## F2S4_170320_001_B01         L4
    ## F2S4_170320_001_C01      Astro
    ## F2S4_170320_001_D01      Astro
    ## F2S4_170320_001_E01        SMC
    ## F2S4_170320_001_F01      Astro
    ## F2S4_170320_001_H01      Astro
    ## F2S4_170320_002_A01    L2/3 IT
    ## F2S4_170320_002_B01      Astro
    ## F2S4_170320_002_C01         NP
    ## F2S4_170320_002_D01        SMC
    ## F2S4_170320_002_E01         NP
    ## F2S4_170320_002_F01        Vip
    ## F2S4_170320_002_G01      Astro
    ## F2S4_170320_002_H01        Vip
    ## F2S4_170320_003_A01      Astro
    ## F2S4_170320_003_B01        SMC
    ## F2S4_170320_003_C01      Astro
    ## F2S4_170320_003_D01      Astro
    ## F2S4_170320_003_E01        SMC
    ## F2S4_170320_003_G01      Astro
    ## F2S4_170320_003_H01      Lamp5
    ## F2S4_170320_004_A01      Astro
    ## F2S4_170320_004_B01        SMC
    ## F2S4_170320_004_C01      Astro
    ## F2S4_170320_004_D01        SMC
    ## F2S4_170320_004_E01        Vip
    ## F2S4_170320_004_F01        Sst
    ## F2S4_170320_004_H01      Astro
    ## F2S4_170320_005_A01      Meis2
    ## F2S4_170320_005_B01      Astro
    ## F2S4_170320_005_C01      Astro
    ## F2S4_170320_005_D01        Vip
    ## F2S4_170320_005_E01        SMC
    ## F2S4_170320_005_F01      Lamp5
    ## F2S4_170320_005_G01        SMC
    ## F2S4_170320_005_H01        SMC
    ## F2S4_170320_006_A01        Vip
    ## F2S4_170320_006_B01        SMC
    ## F2S4_170320_006_C01      Astro
    ## F2S4_170320_006_D01      Pvalb
    ## F2S4_170320_006_E01        Vip
    ## F2S4_170320_006_F01         L4
    ## F2S4_170320_006_G01      Astro
    ## F2S4_170320_006_H01        SMC
    ## F2S4_170320_007_B01      Astro
    ## F2S4_170320_007_C01      Astro
    ## F2S4_170320_007_D01       VLMC
    ## F2S4_170320_007_E01      Astro
    ## F2S4_170320_007_F01      Astro
    ## F2S4_170320_007_G01        Vip
    ## F2S4_170320_007_H01        SMC
    ## F2S4_170320_008_A01      Meis2
    ## F2S4_170320_008_B01    L2/3 IT
    ## F2S4_170320_008_C01       VLMC
    ## F2S4_170320_008_D01      Meis2
    ## F2S4_170320_008_F01       VLMC
    ## F2S4_170320_008_G01        Vip
    ## F2S4_170320_008_H01        SMC
    ## F2S4_170323_001_A01 Macrophage
    ## F2S4_170323_001_B01      Astro
    ## F2S4_170323_001_C01      Astro
    ## F2S4_170323_001_D01      Astro
    ## F2S4_170323_001_E01      Astro
    ## F2S4_170323_001_F01      Astro
    ## F2S4_170323_001_H01      Oligo
    ## F2S4_170323_002_B01      Astro
    ## F2S4_170323_002_C01      Astro
    ## F2S4_170323_002_D01      Astro
    ## F2S4_170323_002_E01      Astro
    ## F2S4_170323_002_F01      Oligo
    ## F2S4_170323_002_G01      Oligo
    ## F2S4_170323_003_A01       Endo
    ## F2S4_170323_003_B01 Macrophage
    ## F2S4_170323_003_C01      Oligo
    ## F2S4_170323_003_D01      Astro
    ## F2S4_170323_003_E01      Oligo
    ## F2S4_170323_003_F01      Oligo
    ## F2S4_170323_003_G01 Macrophage
    ## F2S4_170323_003_H01      Astro
    ## F2S4_170323_010_A01      Oligo
    ## F2S4_170323_010_C01      Astro
    ## F2S4_170323_010_D01      Oligo
    ## F2S4_170323_010_E01 Macrophage
    ## F2S4_170323_010_F01      Astro
    ## F2S4_170323_010_G01      Oligo
    ## F2S4_170323_011_A01      Oligo
    ## F2S4_170323_011_B01      Astro
    ## F2S4_170323_011_C01      Astro
    ## F2S4_170323_011_D01 Macrophage
    ## F2S4_170323_011_E01      Astro
    ## F2S4_170323_011_F01       Endo
    ## F2S4_170323_011_G01      Oligo
    ## F2S4_170323_012_A01      Astro
    ## F2S4_170323_012_B01        Vip
    ## F2S4_170323_012_C01      Astro
    ## F2S4_170323_012_D01      Oligo
    ## F2S4_170323_012_E01      Astro
    ## F2S4_170323_012_F01 Macrophage
    ## F2S4_170323_012_G01      Astro
    ## F2S4_170323_012_H01      Astro
    ## F2S4_170323_013_B01      Oligo
    ## F2S4_170323_013_C01 Macrophage
    ## F2S4_170323_013_D01 Macrophage
    ## F2S4_170323_013_E01 Macrophage
    ## F2S4_170323_013_F01 Macrophage
    ## F2S4_170323_013_H01      Astro
    ## F2S4_170323_014_A01      Astro
    ## F2S4_170323_014_B01      Oligo
    ## F2S4_170323_014_C01      Oligo
    ## F2S4_170323_014_D01 Macrophage
    ## F2S4_170323_014_E01      Oligo
    ## F2S4_170323_014_F01      Astro
    ## F2S4_170323_014_H01 Macrophage
    ## F2S4_170323_015_B01 Macrophage
    ## F2S4_170323_015_C01      Astro
    ## F2S4_170323_015_D01 Macrophage
    ## F2S4_170323_015_E01      Astro
    ## F2S4_170323_015_F01      Astro
    ## F2S4_170323_015_G01       Peri
    ## F2S4_170323_015_H01 Macrophage
    ## F2S4_170323_016_A01      Astro
    ## F2S4_170323_016_B01      Astro
    ## F2S4_170323_016_D01      Oligo
    ## F2S4_170323_016_F01 Macrophage
    ## F2S4_170323_016_G01 Macrophage
    ## F2S4_170323_016_H01 Macrophage
    ## F2S4_170323_017_B01      Astro
    ## F2S4_170323_017_C01      Oligo
    ## F2S4_170323_017_D01 Macrophage
    ## F2S4_170323_017_E01       Peri
    ## F2S4_170323_017_F01      Oligo
    ## F2S4_170323_017_G01 Macrophage
    ## F2S4_170323_017_H01 Macrophage
    ## F2S4_170323_018_A01 Macrophage
    ## F2S4_170323_018_B01      Astro
    ## F2S4_170323_018_D01      Astro
    ## F2S4_170323_018_F01 Macrophage
    ## F2S4_170323_018_G01       Peri
    ## F2S4_170323_018_H01 Macrophage
    ## F2S4_170323_025_A01      Oligo
    ## F2S4_170323_025_B01      Oligo
    ## F2S4_170323_025_C01      Astro
    ## F2S4_170323_025_D01      Astro
    ## F2S4_170323_025_E01      Astro
    ## F2S4_170323_025_F01      Astro
    ## F2S4_170323_025_G01 Macrophage
    ## F2S4_170323_025_H01 Macrophage
    ## F2S4_170323_026_B01      Oligo
    ## F2S4_170323_026_C01      Astro
    ## F2S4_170323_026_D01      Oligo
    ## F2S4_170323_026_E01      Oligo
    ## F2S4_170323_026_F01 Macrophage
    ## F2S4_170323_026_G01      Astro
    ## F2S4_170323_026_H01      Astro
    ## F2S4_170323_027_A01       Endo
    ## F2S4_170323_027_B01      Astro
    ## F2S4_170323_027_D01      Oligo
    ## F2S4_170323_027_E01      Oligo
    ## F2S4_170323_027_F01      Oligo
    ## F2S4_170323_027_G01      Oligo
    ## F2S4_170323_027_H01       Endo
    ## F2S4_170323_028_B01      Oligo
    ## F2S4_170323_028_C01      Astro
    ## F2S4_170323_028_D01      Astro
    ## F2S4_170323_028_E01      Oligo
    ## F2S4_170323_028_F01      Astro
    ## F2S4_170323_028_H01      Astro
    ## F2S4_170323_035_A01      Astro
    ## F2S4_170323_035_B01      Oligo
    ## F2S4_170323_035_D01      Oligo
    ## F2S4_170323_035_E01 Macrophage
    ## F2S4_170323_035_F01 Macrophage
    ## F2S4_170323_035_G01      Astro
    ## F2S4_170323_036_B01      Oligo
    ## F2S4_170323_036_C01      Astro
    ## F2S4_170323_036_E01      Astro
    ## F2S4_170323_036_F01      Astro
    ## F2S4_170323_036_H01      Oligo
    ## F2S4_170323_037_A01      Oligo
    ## F2S4_170323_037_B01      Astro
    ## F2S4_170323_037_C01      Oligo
    ## F2S4_170323_037_D01 Macrophage
    ## F2S4_170323_037_E01      Oligo
    ## F2S4_170323_037_G01 Macrophage
    ## F2S4_170323_037_H01      Oligo
    ## F2S4_170323_038_A01      Oligo
    ## F2S4_170323_038_B01      Astro
    ## F2S4_170323_038_C01      Oligo
    ## F2S4_170323_038_D01 Macrophage
    ## F2S4_170323_038_E01      Oligo
    ## F2S4_170323_038_F01      Oligo
    ## F2S4_170323_038_G01      Astro
    ## F2S4_170323_038_H01      Astro
    ## F2S4_170323_039_A01      Astro
    ## F2S4_170323_039_B01      Oligo
    ## F2S4_170323_039_C01      Oligo
    ## F2S4_170323_039_D01      Oligo
    ## F2S4_170323_039_E01      Astro
    ## F2S4_170323_039_F01      Oligo
    ## F2S4_170323_039_G01      Oligo
    ## F2S4_170323_039_H01      Astro
    ## F2S4_170323_040_A01      Oligo
    ## F2S4_170323_040_B01 Macrophage
    ## F2S4_170323_040_C01      Oligo
    ## F2S4_170323_040_E01      Oligo
    ## F2S4_170323_040_F01      Oligo
    ## F2S4_170323_040_G01      Astro
    ## F2S4_170323_040_H01      Astro
    ## F2S4_170323_041_A01      Astro
    ## F2S4_170323_041_B01      Oligo
    ## F2S4_170323_041_D01      Astro
    ## F2S4_170323_041_E01      Oligo
    ## F2S4_170323_041_F01      Oligo
    ## F2S4_170323_041_G01      Astro
    ## F2S4_170323_041_H01      Oligo
    ## F2S4_170323_048_A01      Oligo
    ## F2S4_170323_048_B01      Oligo
    ## F2S4_170323_048_C01      Astro
    ## F2S4_170323_048_D01      Oligo
    ## F2S4_170323_048_E01       Peri
    ## F2S4_170323_048_F01      Oligo
    ## F2S4_170323_048_G01      Oligo
    ## F2S4_170323_048_H01      Oligo
    ## F2S4_170323_049_A01 Macrophage
    ## F2S4_170323_049_B01      Astro
    ## F2S4_170323_049_C01      Astro
    ## F2S4_170323_049_D01      Astro
    ## F2S4_170323_049_E01      Oligo
    ## F2S4_170323_049_F01 Macrophage
    ## F2S4_170323_049_G01      Astro
    ## F2S4_170323_049_H01      Astro
    ## F2S4_170323_050_A01      Astro
    ## F2S4_170323_050_B01      Astro
    ## F2S4_170323_050_C01      Oligo
    ## F2S4_170323_050_D01      Oligo
    ## F2S4_170323_050_E01 Macrophage
    ## F2S4_170323_050_F01      Oligo
    ## F2S4_170323_050_G01      Oligo
    ## F2S4_170323_050_H01      Astro
    ## F2S4_170323_051_A01      Oligo
    ## F2S4_170323_051_B01      Oligo
    ## F2S4_170323_051_C01      Astro
    ## F2S4_170323_051_D01      Oligo
    ## F2S4_170323_051_E01 Macrophage
    ## F2S4_170323_051_F01      Astro
    ## F2S4_170323_051_G01       Peri
    ## F2S4_170323_051_H01 Macrophage
    ## F2S4_170323_052_A01      Astro
    ## F2S4_170323_052_B01      Oligo
    ## F2S4_170323_052_C01      Oligo
    ## F2S4_170323_052_D01      Astro
    ## F2S4_170323_052_E01 Macrophage
    ## F2S4_170323_052_F01 Macrophage
    ## F2S4_170323_052_G01      Astro
    ## F2S4_170323_052_H01      Astro
    ## F2S4_170323_059_A01      Astro
    ## F2S4_170323_059_B01      Oligo
    ## F2S4_170323_059_D01      Astro
    ## F2S4_170323_059_E01 Macrophage
    ## F2S4_170323_059_F01      Astro
    ## F2S4_170323_059_G01      Oligo
    ## F2S4_170323_059_H01       Peri
    ## F2S4_170323_060_A01      Oligo
    ## F2S4_170323_060_B01      Astro
    ## F2S4_170323_060_C01      Astro
    ## F2S4_170323_060_D01      Oligo
    ## F2S4_170323_060_E01 Macrophage
    ## F2S4_170323_060_F01      Astro
    ## F2S4_170323_060_G01      Oligo
    ## F2S4_170406_001_A01      Pvalb
    ## F2S4_170406_001_B01      Pvalb
    ## F2S4_170406_001_C01      L6 IT
    ## F2S4_170406_001_E01      L5 PT
    ## F2S4_170406_001_F01      L6 IT
    ## F2S4_170406_001_G01      Pvalb
    ## F2S4_170406_001_H01      Pvalb
    ## F2S4_170406_002_A01      L6 IT
    ## F2S4_170406_002_B01      Pvalb
    ## F2S4_170406_002_C01      Pvalb
    ## F2S4_170406_002_D01      L6 IT
    ## F2S4_170406_002_E01      Pvalb
    ## F2S4_170406_002_F01        Sst
    ## F2S4_170406_002_G01      Pvalb
    ## F2S4_170406_002_H01      L6 IT
    ## F2S4_170406_003_A01      Pvalb
    ## F2S4_170406_003_B01      L6 IT
    ## F2S4_170406_003_C01      Pvalb
    ## F2S4_170406_003_D01      Pvalb
    ## F2S4_170406_003_E01      Pvalb
    ## F2S4_170406_003_F01      Pvalb
    ## F2S4_170406_003_G01      Pvalb
    ## F2S4_170406_003_H01      Pvalb
    ## F2S4_170406_004_A01      Pvalb
    ## F2S4_170406_004_B01      Pvalb
    ## F2S4_170406_004_C01      Astro
    ## F2S4_170406_004_D01      Pvalb
    ## F2S4_170406_004_E01      Pvalb
    ## F2S4_170406_004_F01      Pvalb
    ## F2S4_170406_004_G01      Pvalb
    ## F2S4_170406_004_H01      Pvalb
    ## F2S4_170406_005_A01      Pvalb
    ## F2S4_170406_005_B01      Pvalb
    ## F2S4_170406_005_C01      Pvalb
    ## F2S4_170406_005_D01        L6b
    ## F2S4_170406_005_E01      Pvalb
    ## F2S4_170406_005_F01      Pvalb
    ## F2S4_170406_005_G01      Pvalb
    ## F2S4_170406_005_H01      Pvalb
    ## F2S4_170406_006_A01      L6 IT
    ## F2S4_170406_006_B01      Pvalb
    ## F2S4_170406_006_C01      Pvalb
    ## F2S4_170406_006_D01      Pvalb
    ## F2S4_170406_006_E01        L6b
    ## F2S4_170406_006_F01      Pvalb
    ## F2S4_170406_006_G01      Pvalb
    ## F2S4_170406_007_A01      Pvalb
    ## F2S4_170406_007_B01      L5 PT
    ## F2S4_170406_007_C01      Pvalb
    ## F2S4_170406_007_D01      Lamp5
    ## F2S4_170406_007_E01      L6 IT
    ## F2S4_170406_007_F01      Pvalb
    ## F2S4_170406_007_G01      Pvalb
    ## F2S4_170406_007_H01      Pvalb
    ## F2S4_170406_008_B01        Sst
    ## F2S4_170406_008_C01      Pvalb
    ## F2S4_170406_008_D01      Lamp5
    ## F2S4_170406_008_E01      L5 PT
    ## F2S4_170406_008_F01        Vip
    ## F2S4_170406_008_G01      Pvalb
    ## F2S4_170406_008_H01      Pvalb
    ## F2S4_170406_009_A01      Lamp5
    ## F2S4_170406_009_B01      Pvalb
    ## F2S4_170406_009_C01      Pvalb
    ## F2S4_170406_009_D01      Pvalb
    ## F2S4_170406_009_E01      Pvalb
    ## F2S4_170406_009_F01      Pvalb
    ## F2S4_170406_009_G01      Pvalb
    ## F2S4_170406_009_H01      Pvalb
    ## F2S4_170406_010_A01      Pvalb
    ## F2S4_170406_010_B01        Sst
    ## F2S4_170406_010_C01        Sst
    ## F2S4_170406_010_D01      Pvalb
    ## F2S4_170406_010_E01      Pvalb
    ## F2S4_170406_010_F01      Pvalb
    ## F2S4_170406_010_G01      L6 IT
    ## F2S4_170406_010_H01   Serpinf1
    ## F2S4_170406_011_A01        Sst
    ## F2S4_170406_011_B01        Sst
    ## F2S4_170406_011_C01      Pvalb
    ## F2S4_170406_011_D01      Pvalb
    ## F2S4_170406_011_E01      Pvalb
    ## F2S4_170406_011_F01      Pvalb
    ## F2S4_170406_011_G01      Pvalb
    ## F2S4_170406_011_H01      Pvalb
    ## F2S4_170406_012_A01      Pvalb
    ## F2S4_170406_012_B01      Pvalb
    ## F2S4_170406_012_C01      Pvalb
    ## F2S4_170406_012_D01      Pvalb
    ## F2S4_170406_012_E01      Pvalb
    ## F2S4_170406_012_F01      Pvalb
    ## F2S4_170406_012_G01      Pvalb
    ## F2S4_170406_012_H01      Pvalb
    ## F2S4_170406_013_A01        Vip
    ## F2S4_170406_013_B01        Sst
    ## F2S4_170406_013_C01      Pvalb
    ## F2S4_170406_013_D01        Sst
    ## F2S4_170406_013_E01      Pvalb
    ## F2S4_170406_013_F01      Pvalb
    ## F2S4_170406_013_H01      Pvalb
    ## F2S4_170406_014_A01      Lamp5
    ## F2S4_170406_014_B01      Pvalb
    ## F2S4_170406_014_C01      Lamp5
    ## F2S4_170406_014_D01      Pvalb
    ## F2S4_170406_014_E01        Sst
    ## F2S4_170406_014_F01        Sst
    ## F2S4_170406_014_G01      Astro
    ## F2S4_170406_014_H01        Sst
    ## F2S4_170406_015_A01      Pvalb
    ## F2S4_170406_015_B01        Sst
    ## F2S4_170406_015_C01        Sst
    ## F2S4_170406_015_D01      Pvalb
    ## F2S4_170406_015_E01        Vip
    ## F2S4_170406_015_F01      Pvalb
    ## F2S4_170406_015_G01      Pvalb
    ## F2S4_170406_015_H01      Pvalb
    ## F2S4_170406_016_A01        Vip
    ## F2S4_170406_016_B01      Pvalb
    ## F2S4_170406_016_C01        Sst
    ## F2S4_170406_016_D01      Pvalb
    ## F2S4_170406_016_E01        Sst
    ## F2S4_170406_016_F01        Vip
    ## F2S4_170406_016_G01      Pvalb
    ## F2S4_170406_016_H01      Pvalb
    ## F2S4_170406_017_A01        Sst
    ## F2S4_170406_017_B01      Pvalb
    ## F2S4_170406_017_C01        Vip
    ## F2S4_170406_017_D01      Pvalb
    ## F2S4_170406_017_E01        Vip
    ## F2S4_170406_017_F01      Pvalb
    ## F2S4_170406_017_G01        Vip
    ## F2S4_170406_017_H01      Pvalb
    ## F2S4_170406_018_A01      Pvalb
    ## F2S4_170406_018_B01      Pvalb
    ## F2S4_170406_018_C01        Sst
    ## F2S4_170406_018_D01      Pvalb
    ## F2S4_170406_018_E01        Vip
    ## F2S4_170406_018_F01      Pvalb
    ## F2S4_170406_018_G01      Lamp5
    ## F2S4_170410_001_A01      L6 IT
    ## F2S4_170410_001_B01      L6 IT
    ## F2S4_170410_001_C01    L2/3 IT
    ## F2S4_170410_001_D01      L6 IT
    ## F2S4_170410_001_E01      L6 IT
    ## F2S4_170410_001_F01      L6 IT
    ## F2S4_170410_001_G01      L5 IT
    ## F2S4_170410_001_H01    L2/3 IT
    ## F2S4_170410_002_A01      L6 IT
    ## F2S4_170410_002_C01    L2/3 IT
    ## F2S4_170410_002_D01    L2/3 IT
    ## F2S4_170410_002_E01    L2/3 IT
    ## F2S4_170410_002_F01    L2/3 IT
    ## F2S4_170410_002_G01    L2/3 IT
    ## F2S4_170410_002_H01      L6 IT
    ## F2S4_170410_003_A01      L6 IT
    ## F2S4_170410_003_B01      L6 IT
    ## F2S4_170410_003_C01      L6 IT
    ## F2S4_170410_003_D01      L6 IT
    ## F2S4_170410_003_E01    L2/3 IT
    ## F2S4_170410_003_F01      L6 IT
    ## F2S4_170410_003_G01      L6 IT
    ## F2S4_170410_003_H01      L5 IT
    ## F2S4_170410_004_A01      L6 IT
    ## F2S4_170410_004_B01    L2/3 IT
    ## F2S4_170410_004_C01      L5 IT
    ## F2S4_170410_004_D01    L2/3 IT
    ## F2S4_170410_004_E01      L6 IT
    ## F2S4_170410_004_F01      L6 IT
    ## F2S4_170410_004_G01      L5 IT
    ## F2S4_170410_004_H01      L6 IT
    ## F2S4_170410_005_A01      L5 PT
    ## F2S4_170410_005_B01      L5 PT
    ## F2S4_170410_005_C01      L5 PT
    ## F2S4_170410_005_E01      L5 PT
    ## F2S4_170410_005_G01      L5 PT
    ## F2S4_170410_006_A01      L5 IT
    ## F2S4_170410_006_B01      L5 PT
    ## F2S4_170411_001_C01      L6 IT
    ## F2S4_170411_001_D01      L6 IT
    ## F2S4_170411_001_F01      L5 PT
    ## F2S4_170411_001_G01      L6 IT
    ## F2S4_170411_002_A01      L5 PT
    ## F2S4_170411_002_C01      L6 IT
    ## F2S4_170411_002_D01      L5 PT
    ## F2S4_170411_002_F01      L5 PT
    ## F2S4_170411_002_G01      L6 IT
    ## F2S4_170411_003_A01    L2/3 IT
    ## F2S4_170411_003_C01      L6 IT
    ## F2S4_170411_003_D01      L6 IT
    ## F2S4_170411_003_F01      L5 PT
    ## F2S4_170411_004_C01    L2/3 IT
    ## F2S4_170411_004_D01      L6 IT
    ## F2S4_170411_004_F01    L2/3 IT
    ## F2S4_170411_005_A01    L2/3 IT
    ## F2S4_170411_005_B01    L2/3 IT
    ## F2S4_170411_005_C01      L6 IT
    ## F2S4_170411_005_D01      L6 IT
    ## F2S4_170411_005_F01      L5 PT
    ## F2S4_170411_006_B01    L2/3 IT
    ## F2S4_170411_006_C01      L6 IT
    ## F2S4_170411_006_D01      L6 IT
    ## F2S4_170411_006_F01      L5 PT
    ## F2S4_170411_006_G01      L5 PT
    ## F2S4_170411_007_A01      L6 IT
    ## F2S4_170411_007_B01    L2/3 IT
    ## F2S4_170411_007_E01    L2/3 IT
    ## F2S4_170411_007_G01      L5 IT
    ## F2S4_170411_008_A01      L6 IT
    ## F2S4_170411_008_C01      L5 PT
    ## F2S4_170411_008_D01      L5 IT
    ## F2S4_170411_008_E01      L6 IT
    ## F2S4_170411_008_F01      L5 PT
    ## F2S4_170411_008_H01      L5 PT
    ## F2S4_170411_009_C01      L6 IT
    ## F2S4_170411_009_E01      L6 IT
    ## F2S4_170411_009_F01    L2/3 IT
    ## F2S4_170411_010_B01      L5 PT
    ## F2S4_170411_010_C01      L6 IT
    ## F2S4_170411_010_D01    L2/3 IT
    ## F2S4_170411_010_E01      L6 IT
    ## F2S4_170411_010_F01    L2/3 IT
    ## F2S4_170411_010_G01    L2/3 IT
    ## F2S4_170412_001_A01      Pvalb
    ## F2S4_170412_001_B01        Vip
    ## F2S4_170412_001_C01        Sst
    ## F2S4_170412_001_D01      Lamp5
    ## F2S4_170412_001_E01      Meis2
    ## F2S4_170412_001_F01      Meis2
    ## F2S4_170412_001_G01      Meis2
    ## F2S4_170412_001_H01      Meis2
    ## F2S4_170412_002_B01      Meis2
    ## F2S4_170412_002_C01      Meis2
    ## F2S4_170412_002_D01      Meis2
    ## F2S4_170412_002_E01      Meis2
    ## F2S4_170412_002_F01      Meis2
    ## F2S4_170412_002_G01      Pvalb
    ## F2S4_170412_002_H01      Meis2
    ## F2S4_170412_003_A01      Meis2
    ## F2S4_170412_003_B01      Pvalb
    ## F2S4_170412_003_C01      Meis2
    ## F2S4_170412_003_D01      Meis2
    ## F2S4_170412_003_E01      Meis2
    ## F2S4_170412_003_F01      Meis2
    ## F2S4_170412_003_G01      Meis2
    ## F2S4_170412_003_H01        Vip
    ## F2S4_170412_004_A01      Pvalb
    ## F2S4_170412_004_B01      Pvalb
    ## F2S4_170412_004_C01        Vip
    ## F2S4_170412_004_D01      Pvalb
    ## F2S4_170412_004_E01      Pvalb
    ## F2S4_170420_007_A01      Lamp5
    ## F2S4_170420_007_B01      Pvalb
    ## F2S4_170420_007_C01      Lamp5
    ## F2S4_170420_007_D01      Pvalb
    ## F2S4_170420_007_E01      Pvalb
    ## F2S4_170420_007_F01      Pvalb
    ## F2S4_170420_007_G01      Pvalb
    ## F2S4_170420_007_H01      Pvalb
    ## F2S4_170420_008_A01      Pvalb
    ## F2S4_170420_008_B01      Pvalb
    ## F2S4_170420_008_C01      Pvalb
    ## F2S4_170420_008_D01      Pvalb
    ## F2S4_170420_008_E01      Pvalb
    ## F2S4_170420_008_F01      Pvalb
    ## F2S4_170420_008_G01      Pvalb
    ## F2S4_170420_008_H01      Pvalb
    ## F2S4_170420_009_B01        Sst
    ## F2S4_170420_009_C01      Pvalb
    ## F2S4_170420_009_D01        Sst
    ## F2S4_170420_009_E01      Pvalb
    ## F2S4_170420_009_F01      Pvalb
    ## F2S4_170420_009_G01      Pvalb
    ## F2S4_170420_009_H01      Pvalb
    ## F2S4_170420_010_A01      Pvalb
    ## F2S4_170420_010_B01      Pvalb
    ## F2S4_170420_010_C01      Pvalb
    ## F2S4_170420_010_D01        Sst
    ## F2S4_170420_010_E01      Pvalb
    ## F2S4_170420_010_H01      Pvalb
    ## F2S4_170420_011_A01      Pvalb
    ## F2S4_170420_011_B01      Pvalb
    ## F2S4_170420_011_D01      Lamp5
    ## F2S4_170420_011_E01      Pvalb
    ## F2S4_170420_011_F01      Lamp5
    ## F2S4_170420_011_G01      Lamp5
    ## F2S4_170420_011_H01      Lamp5
    ## F2S4_170420_012_A01      Pvalb
    ## F2S4_170420_012_B01      Lamp5
    ## F2S4_170420_012_C01      Lamp5
    ## F2S4_170511_001_C01        Vip
    ## F2S4_170511_001_D01      Pvalb
    ## F2S4_170511_001_E01        Vip
    ## F2S4_170511_001_F01        Vip
    ## F2S4_170511_001_G01      Pvalb
    ## F2S4_170511_001_H01      Pvalb
    ## F2S4_170511_002_B01      Meis2
    ## F2S4_170511_002_E01      Pvalb
    ## F2S4_170511_002_F01      Meis2
    ## F2S4_170516_001_A01         L4
    ## F2S4_170516_001_B01         L4
    ## F2S4_170516_001_C01         L4
    ## F2S4_170516_001_D01         L4
    ## F2S4_170516_001_E01    L2/3 IT
    ## F2S4_170516_001_F01         L4
    ## F2S4_170516_001_G01         L4
    ## F2S4_170516_001_H01         L4
    ## F2S4_170516_002_A01         L4
    ## F2S4_170516_002_B01         L4
    ## F2S4_170516_002_C01         L4
    ## F2S4_170516_002_D01         L4
    ## F2S4_170516_002_E01         L4
    ## F2S4_170516_002_F01         L4
    ## F2S4_170516_002_G01         L4
    ## F2S4_170516_002_H01         L4
    ## F2S4_170516_003_A01         L4
    ## F2S4_170516_003_B01         L4
    ## F2S4_170516_003_C01         L4
    ## F2S4_170516_003_D01         L4
    ## F2S4_170516_003_E01         L4
    ## F2S4_170516_003_F01         L4
    ## F2S4_170516_003_G01         L4
    ## F2S4_170516_003_H01         L4
    ## F2S4_170516_004_A01         L4
    ## F2S4_170516_004_B01         L4
    ## F2S4_170516_004_C01         L4
    ## F2S4_170516_004_D01         L4
    ## F2S4_170516_004_E01         L4
    ## F2S4_170516_004_F01         L4
    ## F2S4_170516_004_G01         L4
    ## F2S4_170516_004_H01         L4
    ## F2S4_170516_005_A01         L4
    ## F2S4_170516_005_B01         L4
    ## F2S4_170516_005_C01         L4
    ## F2S4_170516_005_E01         L4
    ## F2S4_170516_005_F01         L4
    ## F2S4_170516_005_G01         L4
    ## F2S4_170516_005_H01         L4
    ## F2S4_170516_006_A01         L4
    ## F2S4_170516_006_B01         L4
    ## F2S4_170516_006_C01         L4
    ## F2S4_170516_006_D01    L2/3 IT
    ## F2S4_170516_006_E01         L4
    ## F2S4_170516_006_F01         L4
    ## F2S4_170516_006_G01      L5 IT
    ## F2S4_170516_006_H01         L4
    ## F2S4_170516_007_A01         L4
    ## F2S4_170516_007_B01         L4
    ## F2S4_170516_007_C01      L5 IT
    ## F2S4_170516_007_D01         L4
    ## F2S4_170516_007_E01         L4
    ## F2S4_170516_007_F01         L4
    ## F2S4_170516_007_G01      L5 IT
    ## F2S4_170516_007_H01         L4
    ## F2S4_170516_008_A01         L4
    ## F2S4_170516_008_B01         L4
    ## F2S4_170516_008_C01         L4
    ## F2S4_170516_008_D01         L4
    ## F2S4_170516_008_E01         L4
    ## F2S4_170516_008_F01         L4
    ## F2S4_170516_008_G01         L4
    ## F2S4_170516_008_H01         L4
    ## F2S4_170516_009_A01         L4
    ## F2S4_170516_009_B01         L4
    ## F2S4_170516_009_C01      L5 IT
    ## F2S4_170516_009_D01         L4
    ## F2S4_170516_009_E01         L4
    ## F2S4_170516_009_F01         L4
    ## F2S4_170516_009_G01         L4
    ## F2S4_170516_009_H01         L4
    ## F2S4_170516_010_A01         L4
    ## F2S4_170516_010_B01         L4
    ## F2S4_170516_010_C01         L4
    ## F2S4_170516_010_D01      L5 IT
    ## F2S4_170516_010_E01         L4
    ## F2S4_170516_010_G01         L4
    ## F2S4_170516_010_H01         L4
    ## F2S4_170516_011_A01         L4
    ## F2S4_170516_011_B01         L4
    ## F2S4_170516_011_C01         L4
    ## F2S4_170516_011_D01         L4
    ## F2S4_170516_011_E01         L4
    ## F2S4_170516_011_F01      L5 IT
    ## F2S4_170516_011_G01         L4
    ## F2S4_170516_011_H01         L4
    ## F2S4_170516_012_A01         L4
    ## F2S4_170516_012_B01         L4
    ## F2S4_170516_012_C01         L4
    ## F2S4_170516_012_D01         L4
    ## F2S4_170516_012_E01         L4
    ## F2S4_170516_012_F01         L4
    ## F2S4_170516_012_G01         L4
    ## F2S4_170516_013_A01         L4
    ## F2S4_170516_013_B01         L4
    ## F2S4_170516_013_C01         L4
    ## F2S4_170516_013_D01         L4
    ## F2S4_170516_013_E01         L4
    ## F2S4_170516_013_F01         L4
    ## F2S4_170516_013_H01         L4
    ## F2S4_170516_014_G01         L4
    ## F2S4_170516_014_H01         L4
    ## F2S4_170518_001_A01    L2/3 IT
    ## F2S4_170518_001_B01        Vip
    ## F2S4_170518_001_C01        Vip
    ## F2S4_170518_001_D01        Vip
    ## F2S4_170518_001_E01        Vip
    ## F2S4_170518_001_F01         CR
    ## F2S4_170518_001_G01        Vip
    ## F2S4_170518_001_H01        Vip
    ## F2S4_170518_002_A01        Vip
    ## F2S4_170518_002_B01        Vip
    ## F2S4_170518_002_C01        Vip
    ## F2S4_170518_002_D01        Vip
    ## F2S4_170518_002_E01        Vip
    ## F2S4_170518_002_F01        Vip
    ## F2S4_170518_002_G01        Vip
    ## F2S4_170518_002_H01    L2/3 IT
    ## F2S4_170518_006_A01      Meis2
    ## F2S4_170518_006_B01      L6 IT
    ## F2S4_170518_006_C01      L6 IT
    ## F2S4_170518_006_D01      L6 IT
    ## F2S4_170518_006_E01      L6 IT
    ## F2S4_170518_006_F01      L6 IT
    ## F2S4_170518_006_G01        Vip
    ## F2S4_170518_006_H01        Sst
    ## F2S4_170518_007_A01      L6 IT
    ## F2S4_170518_007_B01      L6 IT
    ## F2S4_170518_007_C01      L6 IT
    ## F2S4_170518_007_D01      L6 IT
    ## F2S4_170518_007_E01      L6 IT
    ## F2S4_170518_007_F01      L6 IT
    ## F2S4_170518_007_G01        Vip
    ## F2S4_170518_007_H01      L6 IT
    ## F2S4_170518_011_A01      L6 IT
    ## F2S4_170518_011_B01      Meis2
    ## F2S4_170518_011_C01      L6 IT
    ## F2S4_170518_011_D01      L6 IT
    ## F2S4_170518_011_E01      L6 IT
    ## F2S4_170518_011_F01      L6 IT
    ## F2S4_170518_011_G01      Pvalb
    ## F2S4_170518_011_H01      Meis2
    ## F2S4_170518_012_A01      Meis2
    ## F2S4_170518_012_B01      L6 IT
    ## F2S4_170518_012_C01        Sst
    ## F2S4_170518_012_D01      L6 IT
    ## F2S4_170518_012_E01      L6 IT
    ## F2S4_170518_012_F01      L6 IT
    ## F2S4_170518_012_G01        Sst
    ## F2S4_170518_012_H01      L6 IT
    ## F2S4_170518_016_A01        Vip
    ## F2S4_170518_016_B01        Vip
    ## F2S4_170518_016_C01        Vip
    ## F2S4_170518_016_D01        Vip
    ## F2S4_170518_016_E01        Vip
    ## F2S4_170518_016_F01        Vip
    ## F2S4_170518_016_G01        Vip
    ## F2S4_170518_016_H01        Vip
    ## F2S4_170518_017_A01        Vip
    ## F2S4_170518_017_B01        Vip
    ## F2S4_170518_017_C01        Vip
    ## F2S4_170518_017_D01        Vip
    ## F2S4_170518_017_E01        Vip
    ## F2S4_170518_017_F01        Vip
    ## F2S4_170518_017_G01        Vip
    ## F2S4_170518_017_H01        Vip
    ## F2S4_170519_001_A01        Vip
    ## F2S4_170519_001_B01        Vip
    ## F2S4_170519_001_C01        Vip
    ## F2S4_170519_001_D01        Vip
    ## F2S4_170519_001_E01    L2/3 IT
    ## F2S4_170519_001_F01        Vip
    ## F2S4_170519_001_G01        Vip
    ## F2S4_170519_001_H01    L2/3 IT
    ## F2S4_170519_002_A01        Vip
    ## F2S4_170519_002_B01        Vip
    ## F2S4_170519_002_C01        Vip
    ## F2S4_170519_002_D01        Vip
    ## F2S4_170519_002_E01        Vip
    ## F2S4_170519_002_F01        Vip
    ## F2S4_170519_002_G01        Vip
    ## F2S4_170519_002_H01        Vip
    ## F2S4_170519_003_A01        Vip
    ## F2S4_170519_003_B01        Vip
    ## F2S4_170519_003_C01        Vip
    ## F2S4_170519_003_D01    L2/3 IT
    ## F2S4_170519_003_E01    L2/3 IT
    ## F2S4_170519_003_F01    L2/3 IT
    ## F2S4_170519_003_G01    L2/3 IT
    ## F2S4_170519_003_H01        Vip
    ## F2S4_170519_004_A01        Vip
    ## F2S4_170519_004_B01        Vip
    ## F2S4_170519_004_C01        Vip
    ## F2S4_170519_004_D01        Vip
    ## F2S4_170519_004_E01        Vip
    ## F2S4_170519_004_F01        Vip
    ## F2S4_170519_004_G01    L2/3 IT
    ## F2S4_170519_004_H01      Lamp5
    ## F2S4_170519_005_A01      L6 IT
    ## F2S4_170519_005_B01      L6 IT
    ## F2S4_170519_005_C01      L6 IT
    ## F2S4_170519_005_D01      L6 IT
    ## F2S4_170519_005_E01      L6 IT
    ## F2S4_170519_005_F01      L6 IT
    ## F2S4_170519_005_G01      L6 IT
    ## F2S4_170519_005_H01        Vip
    ## F2S4_170519_006_A01      L6 IT
    ## F2S4_170519_006_B01      L6 IT
    ## F2S4_170519_006_C01      L6 IT
    ## F2S4_170519_006_D01      L6 IT
    ## F2S4_170519_006_E01      Meis2
    ## F2S4_170519_006_F01      L6 IT
    ## F2S4_170519_006_G01      L6 IT
    ## F2S4_170519_006_H01        Vip
    ## F2S4_170519_007_A01      L6 IT
    ## F2S4_170519_007_B01      L6 IT
    ## F2S4_170519_007_C01        Vip
    ## F2S4_170519_007_D01      Meis2
    ## F2S4_170519_007_E01        Vip
    ## F2S4_170519_007_F01      Meis2
    ## F2S4_170519_007_G01        Sst
    ## F2S4_170519_007_H01        Vip
    ## F2S4_170519_008_A01      L6 IT
    ## F2S4_170519_008_B01      L6 IT
    ## F2S4_170519_008_C01      L6 IT
    ## F2S4_170519_008_D01      L6 IT
    ## F2S4_170519_008_E01      L6 IT
    ## F2S4_170519_008_F01      L6 IT
    ## F2S4_170519_008_G01      L6 IT
    ## F2S4_170519_008_H01        Vip
    ## F2S4_170526_001_A01         NP
    ## F2S4_170526_001_B01         NP
    ## F2S4_170526_001_C01      L6 IT
    ## F2S4_170526_001_D01        Sst
    ## F2S4_170526_001_E01       Peri
    ## F2S4_170526_001_F01         NP
    ## F2S4_170526_001_G01         NP
    ## F2S4_170526_001_H01         NP
    ## F2S4_170526_002_A01      Lamp5
    ## F2S4_170526_002_B01         NP
    ## F2S4_170526_002_C01         NP
    ## F2S4_170526_002_D01         NP
    ## F2S4_170526_002_E01         NP
    ## F2S4_170526_002_F01      Pvalb
    ## F2S4_170526_002_G01         NP
    ## F2S4_170526_002_H01      Lamp5
    ## F2S4_170526_003_A01         NP
    ## F2S4_170526_003_B01         NP
    ## F2S4_170526_003_C01         NP
    ## F2S4_170526_003_D01      Pvalb
    ## F2S4_170526_003_F01        Sst
    ## F2S4_170526_003_G01        Sst
    ## F2S4_170526_003_H01      Pvalb
    ## F2S4_170526_004_A01      L6 IT
    ## F2S4_170526_004_B01         NP
    ## F2S4_170526_004_C01      Pvalb
    ## F2S4_170526_004_D01         NP
    ## F2S4_170526_004_E01      Pvalb
    ## F2S4_170526_004_G01      L6 IT
    ## F2S4_170526_004_H01      L6 IT
    ## F2S4_170526_005_A01      L6 IT
    ## F2S4_170526_005_B01         NP
    ## F2S4_170526_005_C01         NP
    ## F2S4_170526_005_D01      Pvalb
    ## F2S4_170526_005_E01      Pvalb
    ## F2S4_170526_005_G01         NP
    ## F2S4_170526_006_A01         NP
    ## F2S4_170526_006_B01      L6 IT
    ## F2S4_170526_006_C01        Sst
    ## F2S4_170526_006_D01      L6 IT
    ## F2S4_170526_006_E01         NP
    ## F2S4_170526_006_F01      L6 CT
    ## F2S4_170526_006_H01         NP
    ## F2S4_170526_007_A01         NP
    ## F2S4_170526_007_B01        Sst
    ## F2S4_170526_007_C01      Lamp5
    ## F2S4_170526_007_D01         NP
    ## F2S4_170526_007_E01        Sst
    ## F2S4_170526_007_F01      Pvalb
    ## F2S4_170526_007_G01      Pvalb
    ## F2S4_170526_007_H01         NP
    ## F2S4_170526_008_A01        Sst
    ## F2S4_170526_008_B01        Sst
    ## F2S4_170526_008_C01      Pvalb
    ## F2S4_170526_008_D01      Lamp5
    ## F2S4_170526_008_E01        Sst
    ## F2S4_170526_008_F01         NP
    ## F2S4_170526_008_G01        Sst
    ## F2S4_170526_008_H01        Sst
    ## F2S4_170526_009_A01       Sncg
    ## F2S4_170526_009_B01      Pvalb
    ## F2S4_170526_009_C01      L6 IT
    ## F2S4_170526_009_D01        Sst
    ## F2S4_170526_009_E01      L6 CT
    ## F2S4_170526_009_F01      Pvalb
    ## F2S4_170526_009_G01         NP
    ## F2S4_170526_009_H01        Sst
    ## F2S4_170526_010_A01        Sst
    ## F2S4_170526_010_B01      Pvalb
    ## F2S4_170526_010_C01      Pvalb
    ## F2S4_170526_010_D01         NP
    ## F2S4_170526_010_E01      Pvalb
    ## F2S4_170526_010_F01         NP
    ## F2S4_170526_010_G01      L6 IT
    ## F2S4_170526_010_H01      Pvalb
    ## F2S4_170526_011_A01        Sst
    ## F2S4_170526_011_B01      L5 PT
    ## F2S4_170526_011_C01         NP
    ## F2S4_170526_011_D01        Sst
    ## F2S4_170526_011_E01         NP
    ## F2S4_170526_011_F01         NP
    ## F2S4_170526_011_G01         NP
    ## F2S4_170526_011_H01      L6 IT
    ## F2S4_170526_012_B01       Sncg
    ## F2S4_170526_012_D01      Lamp5
    ## F2S4_170526_012_E01       Sncg
    ## F2S4_170526_012_F01         NP
    ## F2S4_170526_012_H01         NP
    ## F2S4_170526_013_A01         NP
    ## F2S4_170526_013_C01         NP
    ## F2S4_170526_013_D01        Sst
    ## F2S4_170526_013_E01        Sst
    ## F2S4_170526_013_F01      Pvalb
    ## F2S4_170526_013_G01        Sst
    ## F2S4_170526_013_H01         NP
    ## F2S4_170526_014_A01      Lamp5
    ## F2S4_170526_014_C01         NP
    ## F2S4_170526_014_D01        Sst
    ## F2S4_170526_014_E01      Lamp5
    ## F2S4_170526_014_F01         NP
    ## F2S4_170526_014_G01         NP
    ## F2S4_170526_014_H01      L5 PT
    ## F2S4_170526_015_A01      L6 IT
    ## F2S4_170526_015_B01         NP
    ## F2S4_170526_015_C01        L6b
    ## F2S4_170526_015_D01      Pvalb
    ## F2S4_170526_015_H01         NP
    ## F2S4_170526_016_A01      Pvalb
    ## F2S4_170526_016_B01         NP
    ## F2S4_170526_016_C01         NP
    ## F2S4_170526_016_E01         NP
    ## F2S4_170526_016_F01      L6 IT
    ## F2S4_170526_016_G01         NP
    ## F2S4_170526_016_H01        Sst
    ## F2S4_170526_017_B01         NP
    ## F2S4_170526_017_C01         NP
    ## F2S4_170526_017_D01         NP
    ## F2S4_170526_017_E01        L6b
    ## F2S4_170526_017_F01      Pvalb
    ## F2S4_170526_017_G01      L6 IT
    ## F2S4_170526_017_H01         NP
    ## F2S4_170526_018_A01      Pvalb
    ## F2S4_170526_018_B01        Sst
    ## F2S4_170526_018_C01      L6 IT
    ## F2S4_170526_018_D01      Pvalb
    ## F2S4_170526_018_E01         NP
    ## F2S4_170526_018_F01      L6 IT
    ## F2S4_170526_018_G01         NP
    ## F2S4_170526_018_H01      Pvalb
    ## F2S4_170526_019_A01         NP
    ## F2S4_170526_019_B01      L6 IT
    ## F2S4_170526_019_E01      Pvalb
    ## F2S4_170526_019_F01         NP
    ## F2S4_170526_019_G01      L6 IT
    ## F2S4_170526_019_H01         NP
    ## F2S4_170526_020_A01         NP
    ## F2S4_170526_020_B01      L6 CT
    ## F2S4_170526_020_C01      L6 IT
    ## F2S4_170526_020_D01         NP
    ## F2S4_170526_020_E01         NP
    ## F2S4_170526_020_F01      Pvalb
    ## F2S4_170526_020_G01         NP
    ## F2S4_170526_020_H01         NP
    ## F2S4_170605_001_A01        L6b
    ## F2S4_170605_001_B01        L6b
    ## F2S4_170605_001_C01        L6b
    ## F2S4_170605_001_D01        L6b
    ## F2S4_170605_001_E01        L6b
    ## F2S4_170605_001_F01        L6b
    ## F2S4_170605_001_G01        L6b
    ## F2S4_170605_001_H01        L6b
    ## F2S4_170605_002_B01        L6b
    ## F2S4_170605_002_C01        L6b
    ## F2S4_170605_002_D01        L6b
    ## F2S4_170605_002_E01        L6b
    ## F2S4_170605_002_F01        L6b
    ## F2S4_170605_002_G01        L6b
    ## F2S4_170605_002_H01        L6b
    ## F2S4_170605_003_A01        L6b
    ## F2S4_170605_003_B01        L6b
    ## F2S4_170605_003_C01        L6b
    ## F2S4_170605_003_D01        L6b
    ## F2S4_170605_003_E01        L6b
    ## F2S4_170605_003_F01        L6b
    ## F2S4_170605_003_G01        L6b
    ## F2S4_170605_003_H01        L6b
    ## F2S4_170605_004_A01        L6b
    ## F2S4_170605_004_B01        L6b
    ## F2S4_170605_004_C01        L6b
    ## F2S4_170605_004_D01        L6b
    ## F2S4_170605_004_E01        L6b
    ## F2S4_170605_004_F01        L6b
    ## F2S4_170605_004_G01        L6b
    ## F2S4_170605_004_H01        L6b
    ## F2S4_170605_005_A01        L6b
    ## F2S4_170605_005_B01        L6b
    ## F2S4_170605_005_C01        L6b
    ## F2S4_170605_005_D01        L6b
    ## F2S4_170605_005_E01        L6b
    ## F2S4_170605_005_F01        L6b
    ## F2S4_170605_005_G01        L6b
    ## F2S4_170605_005_H01        L6b
    ## F2S4_170605_006_A01        L6b
    ## F2S4_170605_006_B01        L6b
    ## F2S4_170605_006_C01        L6b
    ## F2S4_170605_006_D01        L6b
    ## F2S4_170605_006_E01        L6b
    ## F2S4_170605_006_F01        L6b
    ## F2S4_170605_006_G01        L6b
    ## F2S4_170605_006_H01        L6b
    ## F2S4_170605_007_A01        L6b
    ## F2S4_170605_007_B01        L6b
    ## F2S4_170605_007_C01        L6b
    ## F2S4_170605_007_D01        L6b
    ## F2S4_170605_007_E01        L6b
    ## F2S4_170605_007_F01        L6b
    ## F2S4_170605_007_G01        L6b
    ## F2S4_170605_007_H01        L6b
    ## F2S4_170605_008_A01        L6b
    ## F2S4_170605_008_B01        L6b
    ## F2S4_170605_008_C01        L6b
    ## F2S4_170605_008_D01        L6b
    ## F2S4_170605_008_E01        L6b
    ## F2S4_170605_008_F01        L6b
    ## F2S4_170605_008_G01        L6b
    ## F2S4_170605_008_H01        L6b
    ## F2S4_170605_009_A01        L6b
    ## F2S4_170605_009_B01        L6b
    ## F2S4_170605_009_C01        L6b
    ## F2S4_170605_009_D01        L6b
    ## F2S4_170605_009_E01        L6b
    ## F2S4_170605_009_F01        L6b
    ## F2S4_170605_009_G01        L6b
    ## F2S4_170605_009_H01        L6b
    ## F2S4_170605_010_A01        L6b
    ## F2S4_170605_010_B01        L6b
    ## F2S4_170605_010_C01        L6b
    ## F2S4_170605_010_D01        L6b
    ## F2S4_170605_010_E01        L6b
    ## F2S4_170605_010_F01        L6b
    ## F2S4_170605_010_G01        L6b
    ## F2S4_170605_010_H01        L6b
    ## F2S4_170605_011_B01        L6b
    ## F2S4_170605_011_C01        L6b
    ## F2S4_170605_011_D01        L6b
    ## F2S4_170605_011_E01        L6b
    ## F2S4_170605_011_F01        L6b
    ## F2S4_170605_011_G01        L6b
    ## F2S4_170605_011_H01        L6b
    ## F2S4_170605_012_A01        L6b
    ## F2S4_170605_012_B01        L6b
    ## F2S4_170605_012_C01        L6b
    ## F2S4_170605_012_D01        L6b
    ## F2S4_170605_012_E01        L6b
    ## F2S4_170605_012_F01        L6b
    ## F2S4_170605_012_G01        L6b
    ## F2S4_170605_012_H01        L6b
    ## F2S4_170605_013_A01        L6b
    ## F2S4_170605_013_B01        L6b
    ## F2S4_170605_013_C01        L6b
    ## F2S4_170605_013_D01        L6b
    ## F2S4_170605_013_E01        L6b
    ## F2S4_170605_013_F01        L6b
    ## F2S4_170605_013_G01        L6b
    ## F2S4_170605_013_H01        L6b
    ## F2S4_170605_014_A01        L6b
    ## F2S4_170605_014_B01        L6b
    ## F2S4_170605_014_C01        L6b
    ## F2S4_170605_014_D01        L6b
    ## F2S4_170605_014_E01        L6b
    ## F2S4_170605_014_F01        L6b
    ## F2S4_170605_014_G01        L6b
    ## F2S4_170605_014_H01        L6b
    ## F2S4_170605_015_A01        L6b
    ## F2S4_170605_015_B01        L6b
    ## F2S4_170605_015_C01        L6b
    ## F2S4_170605_015_D01      L6 CT
    ## F2S4_170605_015_E01        L6b
    ## F2S4_170605_015_G01        L6b
    ## F2S4_170605_015_H01        L6b
    ## F2S4_170605_016_A01        L6b
    ## F2S4_170605_016_B01        L6b
    ## F2S4_170605_016_C01        L6b
    ## F2S4_170605_016_D01        L6b
    ## F2S4_170605_016_E01        L6b
    ## F2S4_170605_016_F01        L6b
    ## F2S4_170605_016_G01        L6b
    ## F2S4_170605_016_H01        L6b
    ## F2S4_170605_017_B01        L6b
    ## F2S4_170605_017_C01        L6b
    ## F2S4_170605_017_D01        L6b
    ## F2S4_170605_017_E01        L6b
    ## F2S4_170605_017_F01        L6b
    ## F2S4_170605_017_G01        L6b
    ## F2S4_170605_017_H01        L6b
    ## F2S4_170605_018_C01        L6b
    ## F2S4_170605_018_D01        L6b
    ## F2S4_170605_018_E01        L6b
    ## F2S4_170605_018_F01        L6b
    ## F2S4_170605_018_G01        L6b
    ## F2S4_170605_018_H01        L6b
    ## F2S4_170609_001_A01      L5 PT
    ## F2S4_170609_001_C01      L5 PT
    ## F2S4_170609_001_D01      L5 PT
    ## F2S4_170609_001_E01      L5 PT
    ## F2S4_170609_001_F01      L5 PT
    ## F2S4_170609_001_H01        L6b
    ## F2S4_170609_002_A01        L6b
    ## F2S4_170609_002_B01        L6b
    ## F2S4_170609_002_C01        L6b
    ## F2S4_170609_002_D01        L6b
    ## F2S4_170609_002_E01        L6b
    ## F2S4_170609_002_F01        L6b
    ## F2S4_170609_002_G01        L6b
    ## F2S4_170609_002_H01        L6b
    ## F2S4_170609_003_A01        L6b
    ## F2S4_170609_003_B01        L6b
    ## F2S4_170609_003_C01        L6b
    ## F2S4_170609_003_D01        L6b
    ## F2S4_170609_003_E01        L6b
    ## F2S4_170609_003_F01        L6b
    ## F2S4_170609_003_G01        L6b
    ## F2S4_170609_003_H01        L6b
    ## F2S4_170609_004_A01    L2/3 IT
    ## F2S4_170609_004_B01      L5 IT
    ## F2S4_170609_004_C01      L6 IT
    ## F2S4_170609_004_D01    L2/3 IT
    ## F2S4_170609_004_E01      L6 IT
    ## F2S4_170609_004_F01    L2/3 IT
    ## F2S4_170609_004_G01    L2/3 IT
    ## F2S4_170609_004_H01    L2/3 IT
    ## F2S4_170609_005_A01    L2/3 IT
    ## F2S4_170609_005_B01    L2/3 IT
    ## F2S4_170609_005_C01      L6 IT
    ## F2S4_170609_005_D01      L6 IT
    ## F2S4_170609_005_E01    L2/3 IT
    ## F2S4_170609_005_F01      L6 IT
    ## F2S4_170609_005_G01      L6 IT
    ## F2S4_170609_005_H01    L2/3 IT
    ## F2S4_170609_006_A01      L6 IT
    ## F2S4_170609_006_B01    L2/3 IT
    ## F2S4_170609_006_C01    L2/3 IT
    ## F2S4_170609_006_D01      L6 IT
    ## F2S4_170609_006_E01    L2/3 IT
    ## F2S4_170609_006_F01    L2/3 IT
    ## F2S4_170609_006_G01      L6 IT
    ## F2S4_170609_006_H01    L2/3 IT
    ## F2S4_170609_007_A01        L6b
    ## F2S4_170609_007_B01        L6b
    ## F2S4_170609_007_C01        L6b
    ## F2S4_170609_007_D01    L2/3 IT
    ## F2S4_170609_007_E01    L2/3 IT
    ## F2S4_170609_007_F01    L2/3 IT
    ## F2S4_170609_007_G01    L2/3 IT
    ## F2S4_170609_007_H01    L2/3 IT
    ## F2S4_170612_001_A01      L6 IT
    ## F2S4_170612_001_B01      L6 IT
    ## F2S4_170612_001_C01      L6 IT
    ## F2S4_170612_001_E01      L5 PT
    ## F2S4_170612_001_F01      L5 PT
    ## F2S4_170612_001_G01      L5 PT
    ## F2S4_170612_001_H01      L5 PT
    ## F2S4_170612_002_A01    L2/3 IT
    ## F2S4_170706_001_A01         L4
    ## F2S4_170706_001_B01      Astro
    ## F2S4_170706_001_E01      Astro
    ## F2S4_170707_003_A01      Lamp5
    ## F2S4_170707_003_B01      Lamp5
    ## F2S4_170707_003_C01        Vip
    ## F2S4_170707_003_D01         CR
    ## F2S4_170707_003_E01      Lamp5
    ## F2S4_170707_003_F01      Lamp5
    ## F2S4_170707_003_H01      Lamp5
    ## F2S4_170707_004_A01      Lamp5
    ## F2S4_170707_004_B01         CR
    ## F2S4_170707_004_C01      Lamp5
    ## F2S4_170707_004_D01      Lamp5
    ## F2S4_170707_004_E01      Lamp5
    ## F2S4_170707_004_G01      Lamp5
    ## F2S4_170707_004_H01      Lamp5
    ## F2S4_170707_005_A01        Vip
    ## F2S4_170707_005_B01        Vip
    ## F2S4_170707_005_C01        Vip
    ## F2S4_170707_005_D01        Vip
    ## F2S4_170707_005_E01      Lamp5
    ## F2S4_170707_005_F01        Vip
    ## F2S4_170707_005_G01       Sncg
    ## F2S4_170707_005_H01        Vip
    ## F2S4_170707_006_A01      Lamp5
    ## F2S4_170707_006_B01       Sncg
    ## F2S4_170707_006_C01        Vip
    ## F2S4_170707_006_D01      Lamp5
    ## F2S4_170707_006_E01        Vip
    ## F2S4_170707_006_F01       Sncg
    ## F2S4_170707_006_G01      Lamp5
    ## F2S4_170707_006_H01        Vip
    ## F2S4_170707_011_A01        Vip
    ## F2S4_170707_011_B01        Vip
    ## F2S4_170707_011_C01        Vip
    ## F2S4_170707_011_D01        L6b
    ## F2S4_170707_011_E01        L6b
    ## F2S4_170707_011_F01        Vip
    ## F2S4_170707_011_G01        Vip
    ## F2S4_170707_012_A01        Sst
    ## F2S4_170707_012_B01        Vip
    ## F2S4_170707_012_C01   Serpinf1
    ## F2S4_170707_012_D01        Vip
    ## F2S4_170707_012_E01      L6 CT
    ## F2S4_170707_012_F01        L6b
    ## F2S4_170707_012_G01   Serpinf1
    ## F2S4_170707_012_H01        Vip
    ## F2S4_170707_013_A01        Vip
    ## F2S4_170707_013_B01        Vip
    ## F2S4_170707_013_C01        Vip
    ## F2S4_170707_013_D01        Vip
    ## F2S4_170707_013_E01        Vip
    ## F2S4_170707_013_F01        Vip
    ## F2S4_170707_013_G01        Vip
    ## F2S4_170707_013_H01        Vip
    ## F2S4_170707_014_A01        Vip
    ## F2S4_170707_014_B01        Vip
    ## F2S4_170707_014_C01        Vip
    ## F2S4_170707_014_D01        Vip
    ## F2S4_170707_014_E01        Vip
    ## F2S4_170707_014_F01        Vip
    ## F2S4_170707_014_H01        Vip
    ## F2S4_170707_019_A01      Pvalb
    ## F2S4_170707_019_B01        Sst
    ## F2S4_170707_019_C01      L6 IT
    ## F2S4_170707_019_D01        Vip
    ## F2S4_170707_019_E01        Sst
    ## F2S4_170707_019_F01       Sncg
    ## F2S4_170707_019_G01        Sst
    ## F2S4_170707_019_H01        Sst
    ## F2S4_170707_020_A01        Vip
    ## F2S4_170707_020_B01        Vip
    ## F2S4_170707_020_C01        Vip
    ## F2S4_170707_020_D01       Sncg
    ## F2S4_170707_020_E01        Sst
    ## F2S4_170707_020_F01       Sncg
    ## F2S4_170707_020_G01        Vip
    ## F2S4_170707_020_H01      Pvalb
    ## F2S4_170707_021_C01      Lamp5
    ## F2S4_170707_021_D01      Lamp5
    ## F2S4_170707_021_E01        Vip
    ## F2S4_170707_021_F01        Vip
    ## F2S4_170707_021_G01      Lamp5
    ## F2S4_170707_021_H01      Lamp5
    ## F2S4_170707_022_B01      Lamp5
    ## F2S4_170707_022_C01      Lamp5
    ## F2S4_170707_022_D01      Meis2
    ## F2S4_170707_022_E01      Lamp5
    ## F2S4_170707_022_F01      Lamp5
    ## F2S4_170707_022_G01      Lamp5
    ## F2S4_170707_022_H01      Lamp5
    ## F2S4_170707_027_A01        Sst
    ## F2S4_170707_027_B01       Sncg
    ## F2S4_170707_027_C01        Sst
    ## F2S4_170707_027_D01   Serpinf1
    ## F2S4_170707_027_E01       Sncg
    ## F2S4_170707_027_F01        Sst
    ## F2S4_170707_027_G01        Sst
    ## F2S4_170707_027_H01        Vip
    ## F2S4_170707_028_A01        Sst
    ## F2S4_170707_028_B01        Sst
    ## F2S4_170707_028_C01        Vip
    ## F2S4_170707_028_D01      Pvalb
    ## F2S4_170707_028_E01      Lamp5
    ## F2S4_170707_028_F01       Sncg
    ## F2S4_170707_028_G01        Sst
    ## F2S4_170707_029_A01   Serpinf1
    ## F2S4_170707_029_B01      L6 CT
    ## F2S4_170707_029_C01        L6b
    ## F2S4_170707_029_D01        Vip
    ## F2S4_170707_029_E01        Vip
    ## F2S4_170707_029_F01        Vip
    ## F2S4_170707_029_G01      L6 CT
    ## F2S4_170707_029_H01      Pvalb
    ## F2S4_170707_030_A01      Meis2
    ## F2S4_170707_030_B01        Vip
    ## F2S4_170707_030_C01      Lamp5
    ## F2S4_170707_030_D01        Sst
    ## F2S4_170707_030_E01        Vip
    ## F2S4_170707_030_F01        Sst
    ## F2S4_170707_030_G01        Vip
    ## F2S4_170707_030_H01        L6b
    ## F2S4_170707_035_A01        Vip
    ## F2S4_170707_035_B01        Vip
    ## F2S4_170707_035_C01        Vip
    ## F2S4_170707_035_D01        Vip
    ## F2S4_170707_035_E01        Vip
    ## F2S4_170707_035_F01        Vip
    ## F2S4_170707_035_G01      Lamp5
    ## F2S4_170707_035_H01        Vip
    ## F2S4_170707_036_A01        Vip
    ## F2S4_170707_036_B01      Lamp5
    ## F2S4_170707_036_C01        Vip
    ## F2S4_170707_036_D01        Vip
    ## F2S4_170707_036_E01        Vip
    ## F2S4_170707_036_F01        Vip
    ## F2S4_170707_036_G01      Lamp5
    ## F2S4_170707_036_H01        Vip
    ## F2S4_170707_037_A01        Vip
    ## F2S4_170707_037_B01      Lamp5
    ## F2S4_170707_037_C01        Vip
    ## F2S4_170707_037_D01        Vip
    ## F2S4_170707_037_E01        Vip
    ## F2S4_170707_037_F01        Vip
    ## F2S4_170707_037_G01        Vip
    ## F2S4_170707_037_H01        Vip
    ## F2S4_170707_038_A01        Vip
    ## F2S4_170707_038_B01        Vip
    ## F2S4_170707_038_C01        Vip
    ## F2S4_170707_038_D01        Vip
    ## F2S4_170707_038_E01        Vip
    ## F2S4_170707_038_F01        Vip
    ## F2S4_170707_038_G01        Vip
    ## F2S4_170707_038_H01        Vip
    ## F2S4_170803_001_A01      Pvalb
    ## F2S4_170803_001_B01       VLMC
    ## F2S4_170803_001_C01       Peri
    ## F2S4_170803_001_D01      L6 IT
    ## F2S4_170803_001_E01        Sst
    ## F2S4_170803_001_F01        Sst
    ## F2S4_170803_001_G01       Peri
    ## F2S4_170803_001_H01        Sst
    ## F2S4_170803_002_A01        Sst
    ## F2S4_170803_002_B01        Sst
    ## F2S4_170803_002_C01        Sst
    ## F2S4_170803_002_D01        Sst
    ## F2S4_170803_002_E01        Sst
    ## F2S4_170803_002_F01       VLMC
    ## F2S4_170803_002_G01        Sst
    ## F2S4_170803_002_H01      L6 IT
    ## F2S4_170803_003_A01        Sst
    ## F2S4_170803_003_B01      L5 PT
    ## F2S4_170803_003_C01       Peri
    ## F2S4_170803_003_D01        Sst
    ## F2S4_170803_003_E01        Sst
    ## F2S4_170803_003_F01        Sst
    ## F2S4_170803_003_G01        Sst
    ## F2S4_170803_003_H01      L5 PT
    ## F2S4_170803_004_A01      Pvalb
    ## F2S4_170803_004_C01        Sst
    ## F2S4_170803_004_D01        Sst
    ## F2S4_170803_004_E01        Sst
    ## F2S4_170803_004_F01      Pvalb
    ## F2S4_170803_004_G01      L6 IT
    ## F2S4_170803_004_H01        Sst
    ## F2S4_170803_005_A01      Pvalb
    ## F2S4_170803_005_B01      Pvalb
    ## F2S4_170803_005_C01       Peri
    ## F2S4_170803_005_D01      L6 IT
    ## F2S4_170803_005_E01      L6 IT
    ## F2S4_170803_005_F01      L6 IT
    ## F2S4_170803_005_G01        Sst
    ## F2S4_170803_005_H01      L5 PT
    ## F2S4_170803_006_A01      Pvalb
    ## F2S4_170803_006_B01        Sst
    ## F2S4_170803_006_C01      L6 IT
    ## F2S4_170803_006_D01        Sst
    ## F2S4_170803_006_E01       Peri
    ## F2S4_170803_006_F01        Sst
    ## F2S4_170803_006_G01        Sst
    ## F2S4_170803_007_A01      L6 IT
    ## F2S4_170803_007_B01        Sst
    ## F2S4_170803_007_D01       Peri
    ## F2S4_170803_007_E01      L6 IT
    ## F2S4_170803_007_F01        Sst
    ## F2S4_170803_007_G01      L6 IT
    ## F2S4_170803_007_H01      Pvalb
    ## F2S4_170803_008_A01      L6 IT
    ## F2S4_170803_008_B01      L6 IT
    ## F2S4_170803_008_C01      L6 IT
    ## F2S4_170803_008_D01      L6 IT
    ## F2S4_170803_008_E01      L6 IT
    ## F2S4_170803_008_F01      L6 IT
    ## F2S4_170803_008_G01        Sst
    ## F2S4_170803_008_H01      L6 IT
    ## F2S4_170803_009_A01      L6 IT
    ## F2S4_170803_009_B01      L6 IT
    ## F2S4_170803_009_C01      L6 IT
    ## F2S4_170803_009_D01      L6 IT
    ## F2S4_170803_009_E01        Sst
    ## F2S4_170803_009_F01      L6 IT
    ## F2S4_170803_009_G01      L6 IT
    ## F2S4_170803_009_H01      L6 IT
    ## F2S4_170803_010_A01       Peri
    ## F2S4_170803_010_B01       Peri
    ## F2S4_170803_010_C01        Sst
    ## F2S4_170803_010_D01      L6 IT
    ## F2S4_170803_010_E01      L6 IT
    ## F2S4_170803_010_F01        Sst
    ## F2S4_170803_010_G01       Peri
    ## F2S4_170803_010_H01      L6 IT
    ## F2S4_170803_011_A01       Peri
    ## F2S4_170803_011_B01      L6 IT
    ## F2S4_170803_011_C01       Peri
    ## F2S4_170803_011_D01       Peri
    ## F2S4_170803_011_E01      L6 IT
    ## F2S4_170803_011_F01      L6 IT
    ## F2S4_170803_011_G01       Peri
    ## F2S4_170803_012_A01        Sst
    ## F2S4_170803_012_B01      L6 IT
    ## F2S4_170803_012_C01        Sst
    ## F2S4_170803_012_D01      L6 IT
    ## F2S4_170803_012_E01      L6 IT
    ## F2S4_170803_012_F01      L6 IT
    ## F2S4_170803_012_G01        Sst
    ## F2S4_170803_012_H01      L6 IT
    ## F2S4_170803_013_A01      L6 IT
    ## F2S4_170803_013_B01      L5 PT
    ## F2S4_170803_013_C01        Sst
    ## F2S4_170803_013_D01      L6 IT
    ## F2S4_170803_013_E01        SMC
    ## F2S4_170803_013_F01        Sst
    ## F2S4_170803_013_G01        Sst
    ## F2S4_170803_013_H01      Pvalb
    ## F2S4_170803_014_A01       Peri
    ## F2S4_170803_014_B01        Sst
    ## F2S4_170803_014_C01        Sst
    ## F2S4_170803_014_D01       Peri
    ## F2S4_170803_014_E01        Sst
    ## F2S4_170803_014_F01      L6 IT
    ## F2S4_170803_014_G01      L5 PT
    ## F2S4_170803_014_H01      Pvalb
    ## F2S4_170803_015_A01       Peri
    ## F2S4_170803_015_B01       Peri
    ## F2S4_170803_015_C01      L6 IT
    ## F2S4_170803_015_D01      L6 IT
    ## F2S4_170803_015_E01      L6 IT
    ## F2S4_170803_015_F01      Pvalb
    ## F2S4_170803_015_G01        Sst
    ## F2S4_170803_015_H01      Pvalb
    ## F2S4_170803_016_B01      L6 IT
    ## F2S4_170803_016_C01        Sst
    ## F2S4_170803_016_D01        Sst
    ## F2S4_170803_016_E01      L6 IT
    ## F2S4_170803_016_F01       Sncg
    ## F2S4_170803_016_H01      Pvalb
    ## F2S4_170803_017_A01      Pvalb
    ## F2S4_170803_017_B01        Sst
    ## F2S4_170803_017_C01        Sst
    ## F2S4_170803_017_D01        Sst
    ## F2S4_170803_017_E01       Peri
    ## F2S4_170803_017_F01      Pvalb
    ## F2S4_170803_017_G01       Peri
    ## F2S4_170803_017_H01       Peri
    ## F2S4_170803_018_A01        Sst
    ## F2S4_170803_018_B01        Sst
    ## F2S4_170803_018_C01        Sst
    ## F2S4_170803_018_D01        Sst
    ## F2S4_170803_018_E01      L6 IT
    ## F2S4_170803_018_F01       Peri
    ## F2S4_170803_018_G01      L6 IT
    ## F2S4_170803_018_H01        Sst
    ## F2S4_170803_019_B01      L6 IT
    ## F2S4_170803_019_C01        Sst
    ## F2S4_170803_019_D01      L6 IT
    ## F2S4_170803_019_F01      Pvalb
    ## F2S4_170803_019_G01      L6 IT
    ## F2S4_170803_019_H01      L6 IT
    ## F2S4_170803_020_A01      L6 IT
    ## F2S4_170803_020_B01       Peri
    ## F2S4_170803_020_C01        Sst
    ## F2S4_170803_020_D01      L6 IT
    ## F2S4_170803_020_E01        Sst
    ## F2S4_170803_020_G01        Sst
    ## F2S4_170803_020_H01      L6 IT
    ## F2S4_170803_021_B01       Peri
    ## F2S4_170803_021_D01        Sst
    ## F2S4_170803_021_E01      L6 IT
    ## F2S4_170803_021_F01      L6 IT
    ## F2S4_170803_021_G01      L6 IT
    ## F2S4_170803_021_H01      L6 IT
    ## F2S4_170807_001_G01       Sncg
    ## F2S4_170807_001_H01        Vip
    ## F2S4_170807_005_A01      Lamp5
    ## F2S4_170807_005_B01        Sst
    ## F2S4_170807_005_C01      Lamp5
    ## F2S4_170807_005_D01      Pvalb
    ## F2S4_170807_005_E01        Vip
    ## F2S4_170807_005_F01      Lamp5
    ## F2S4_170807_005_G01        Vip
    ## F2S4_170807_005_H01        Sst
    ## F2S4_170807_006_A01        Vip
    ## F2S4_170807_006_B01    L2/3 IT
    ## F2S4_170807_006_C01        Vip
    ## F2S4_170807_006_D01        Vip
    ## F2S4_170807_006_E01      Pvalb
    ## F2S4_170807_006_F01      Lamp5
    ## F2S4_170807_006_G01        Vip
    ## F2S4_170807_006_H01      Lamp5
    ## F2S4_170807_007_A01        Vip
    ## F2S4_170807_007_B01        Vip
    ## F2S4_170807_007_C01        Vip
    ## F2S4_170807_007_D01      Lamp5
    ## F2S4_170807_007_E01        Vip
    ## F2S4_170807_007_F01      Lamp5
    ## F2S4_170807_007_G01        Vip
    ## F2S4_170807_007_H01        Vip
    ## F2S4_170807_008_A01    L2/3 IT
    ## F2S4_170807_008_B01      Lamp5
    ## F2S4_170807_008_C01        Vip
    ## F2S4_170807_008_D01        Sst
    ## F2S4_170807_008_E01        Vip
    ## F2S4_170807_008_F01        Vip
    ## F2S4_170807_008_G01        Sst
    ## F2S4_170807_008_H01        Vip
    ## F2S4_170807_009_A01        Vip
    ## F2S4_170807_009_B01    L2/3 IT
    ## F2S4_170807_009_C01      Lamp5
    ## F2S4_170807_009_D01        Vip
    ## F2S4_170807_009_E01        Vip
    ## F2S4_170807_009_F01      Lamp5
    ## F2S4_170807_009_G01   Serpinf1
    ## F2S4_170807_009_H01        Vip
    ## F2S4_170807_010_A01      Pvalb
    ## F2S4_170807_010_B01      Pvalb
    ## F2S4_170807_010_C01      Lamp5
    ## F2S4_170807_010_D01        Vip
    ## F2S4_170807_010_E01       Peri
    ## F2S4_170807_010_F01        Vip
    ## F2S4_170807_010_G01        Vip
    ## F2S4_170807_010_H01      Lamp5
    ## F2S4_170807_011_A01        Vip
    ## F2S4_170807_011_B01    L2/3 IT
    ## F2S4_170807_011_C01      Lamp5
    ## F2S4_170807_011_D01      Lamp5
    ## F2S4_170807_011_E01        Vip
    ## F2S4_170807_011_F01      Lamp5
    ## F2S4_170807_011_G01        Vip
    ## F2S4_170807_011_H01        Sst
    ## F2S4_170807_012_A01      Lamp5
    ## F2S4_170807_012_B01        Vip
    ## F2S4_170807_012_C01        Vip
    ## F2S4_170807_012_D01        Vip
    ## F2S4_170807_012_E01    L2/3 IT
    ## F2S4_170807_012_F01      Lamp5
    ## F2S4_170807_012_G01        Sst
    ## F2S4_170807_012_H01        Vip
    ## F2S4_170807_013_A01        Sst
    ## F2S4_170807_013_B01      L6 IT
    ## F2S4_170807_013_C01        Sst
    ## F2S4_170807_013_D01      Pvalb
    ## F2S4_170807_013_E01        Sst
    ## F2S4_170807_013_F01        Sst
    ## F2S4_170807_013_G01        Sst
    ## F2S4_170807_013_H01      L6 IT
    ## F2S4_170807_014_A01        Sst
    ## F2S4_170807_014_B01      L6 IT
    ## F2S4_170807_014_C01      Pvalb
    ## F2S4_170807_014_D01        Sst
    ## F2S4_170807_014_E01        Vip
    ## F2S4_170807_014_F01        Sst
    ## F2S4_170807_014_G01        Sst
    ## F2S4_170807_014_H01        Sst
    ## F2S4_170807_015_A01        Sst
    ## F2S4_170807_015_B01        Sst
    ## F2S4_170807_015_C01        Sst
    ## F2S4_170807_015_D01      L6 IT
    ## F2S4_170807_015_E01      Pvalb
    ## F2S4_170807_015_F01      L6 IT
    ## F2S4_170807_015_G01        Sst
    ## F2S4_170807_015_H01        Sst
    ## F2S4_170807_016_A01      L6 IT
    ## F2S4_170807_016_B01        Sst
    ## F2S4_170807_016_C01        Sst
    ## F2S4_170807_016_D01        Sst
    ## F2S4_170807_016_E01        Sst
    ## F2S4_170807_016_F01        Sst
    ## F2S4_170807_016_G01      L6 IT
    ## F2S4_170807_016_H01      L6 IT
    ## F2S4_170807_020_B01        Sst
    ## F2S4_170807_020_C01      Lamp5
    ## F2S4_170807_020_D01        Sst
    ## F2S4_170807_020_E01        Sst
    ## F2S4_170807_021_A01      Lamp5
    ## F2S4_170807_021_B01        Sst
    ## F2S4_170807_021_C01      L6 IT
    ## F2S4_170807_021_D01      L6 IT
    ## F2S4_170807_021_E01      L6 IT
    ## F2S4_170807_021_F01      L6 IT
    ## F2S4_170807_021_H01        Sst
    ## F2S4_170807_022_A01      L6 IT
    ## F2S4_170807_022_B01        Sst
    ## F2S4_170807_022_C01        Sst
    ## F2S4_170807_022_D01        Sst
    ## F2S4_170807_022_E01        Sst
    ## F2S4_170807_022_F01      L6 IT
    ## F2S4_170807_022_G01        Sst
    ## F2S4_170807_022_H01        Sst
    ## F2S4_170807_023_A01      L6 IT
    ## F2S4_170807_023_B01        Sst
    ## F2S4_170807_023_C01        Sst
    ## F2S4_170807_023_D01      L6 IT
    ## F2S4_170807_023_E01        Sst
    ## F2S4_170807_023_F01      L6 IT
    ## F2S4_170807_023_G01        Sst
    ## F2S4_170807_023_H01        Sst
    ## F2S4_170807_024_A01      L6 IT
    ## F2S4_170807_024_B01      L6 IT
    ## F2S4_170807_024_C01        Sst
    ## F2S4_170807_024_D01      Meis2
    ## F2S4_170807_024_E01      L6 IT
    ## F2S4_170807_024_F01      L6 IT
    ## F2S4_170807_024_G01        Sst
    ## F2S4_170807_024_H01      Pvalb
    ## F2S4_170828_001_A01        Vip
    ## F2S4_170828_001_B01      L6 IT
    ## F2S4_170828_001_C01        Sst
    ## F2S4_170828_001_D01      Pvalb
    ## F2S4_170828_001_E01        Sst
    ## F2S4_170828_001_F01      Pvalb
    ## F2S4_170828_001_G01      L5 PT
    ## F2S4_170828_001_H01      L5 PT
    ## F2S4_170828_002_A01       Sncg
    ## F2S4_170828_002_B01      Pvalb
    ## F2S4_170828_002_C01      Pvalb
    ## F2S4_170828_002_D01        Vip
    ## F2S4_170828_002_E01      L5 PT
    ## F2S4_170828_002_G01         NP
    ## F2S4_170828_002_H01       VLMC
    ## F2S4_170828_003_A01        Sst
    ## F2S4_170828_003_B01      Lamp5
    ## F2S4_170828_003_C01      L5 PT
    ## F2S4_170828_003_F01      L5 PT
    ## F2S4_170828_003_G01      L6 IT
    ## F2S4_170828_003_H01      L6 IT
    ## F2S4_170828_004_A01        Sst
    ## F2S4_170828_004_B01      L5 PT
    ## F2S4_170828_004_C01      Pvalb
    ## F2S4_170828_004_D01      L5 PT
    ## F2S4_170828_004_E01      Pvalb
    ## F2S4_170828_004_F01        Sst
    ## F2S4_170828_004_G01        Sst
    ## F2S4_170828_004_H01      Pvalb
    ## F2S4_170828_005_A01       Sncg
    ## F2S4_170828_005_B01      L5 PT
    ## F2S4_170828_005_D01        Sst
    ## F2S4_170828_005_E01        Sst
    ## F2S4_170828_005_F01      Pvalb
    ## F2S4_170828_005_G01        Sst
    ## F2S4_170828_005_H01       Sncg
    ## F2S4_170828_006_A01      L5 PT
    ## F2S4_170828_006_B01      Lamp5
    ## F2S4_170828_006_C01        Sst
    ## F2S4_170828_006_D01         NP
    ## F2S4_170828_006_E01         NP
    ## F2S4_170828_006_F01      L5 PT
    ## F2S4_170828_006_G01      Lamp5
    ## F2S4_170828_006_H01      Pvalb
    ## F2S4_170828_007_A01   Serpinf1
    ## F2S4_170828_007_B01      L5 PT
    ## F2S4_170828_007_C01      L5 PT
    ## F2S4_170828_007_D01      Lamp5
    ## F2S4_170828_007_E01      L5 PT
    ## F2S4_170828_007_F01        Sst
    ## F2S4_170828_007_G01        Vip
    ## F2S4_170828_007_H01        Vip
    ## F2S4_170828_008_A01         L4
    ## F2S4_170828_008_B01      L5 PT
    ## F2S4_170828_008_C01      Pvalb
    ## F2S4_170828_008_D01         NP
    ## F2S4_170828_008_E01      L5 PT
    ## F2S4_170828_008_F01      L5 PT
    ## F2S4_170828_008_G01        Sst
    ## F2S4_170828_008_H01         NP
    ## F2S4_170830_001_D01        Sst
    ## F2S4_170830_001_E01        Sst
    ## F2S4_170830_001_F01        Sst
    ## F2S4_170830_001_G01        Sst
    ## F2S4_170830_001_H01        Sst
    ## F2S4_170830_002_A01        Sst
    ## F2S4_170830_002_B01        Sst
    ## F2S4_170830_002_D01        Sst
    ## F2S4_170830_002_E01        Sst
    ## F2S4_170830_002_F01        Sst
    ## F2S4_170830_002_G01        Sst
    ## F2S4_170830_003_A01        Sst
    ## F2S4_170830_003_B01        Sst
    ## F2S4_170830_003_C01        Sst
    ## F2S4_170830_003_D01        Sst
    ## F2S4_170830_003_E01      Pvalb
    ## F2S4_170830_003_F01        Sst
    ## F2S4_170830_003_G01        Sst
    ## F2S4_170830_003_H01        Sst
    ## F2S4_170830_004_G01        Sst
    ## F2S4_170830_004_H01        Sst
    ## F2S4_170907_001_A01        Sst
    ## F2S4_170907_001_B01        Sst
    ## F2S4_170907_001_C01        Sst
    ## F2S4_170907_001_D01      Pvalb
    ## F2S4_170907_001_E01      Pvalb
    ## F2S4_170907_001_F01        Sst
    ## F2S4_170907_001_G01        Sst
    ## F2S4_170907_001_H01        Sst
    ## F2S4_170907_002_A01        Sst
    ## F2S4_170907_002_B01      Pvalb
    ## F2S4_170907_002_C01      Pvalb
    ## F2S4_170907_002_D01      Pvalb
    ## F2S4_170907_002_E01        Sst
    ## F2S4_170907_002_F01        Sst
    ## F2S4_170907_002_G01        Sst
    ## F2S4_170907_002_H01        Sst
    ## F2S4_170907_003_A01        Sst
    ## F2S4_170907_003_B01        Sst
    ## F2S4_170907_003_C01        Sst
    ## F2S4_170907_003_D01        Sst
    ## F2S4_170907_003_E01        Sst
    ## F2S4_170907_003_F01        Sst
    ## F2S4_170907_003_G01        Sst
    ## F2S4_170907_003_H01        Sst
    ## F2S4_170907_004_A01      Pvalb
    ## F2S4_170907_004_B01      Pvalb
    ## F2S4_170907_004_C01      Pvalb
    ## F2S4_170907_004_D01        Sst
    ## F2S4_170907_004_E01      Pvalb
    ## F2S4_170907_004_F01      Pvalb
    ## F2S4_170907_004_G01        Sst
    ## F2S4_170907_004_H01        Sst
    ## F2S4_170907_005_A01      Pvalb
    ## F2S4_170908_001_D01        Sst
    ## F2S4_170908_001_E01        Sst
    ## F2S4_170908_001_F01        Sst
    ## F2S4_170908_001_G01        Sst
    ## F2S4_170908_001_H01        Sst
    ## F2S4_170908_002_A01        Sst
    ## F2S4_170908_002_B01        Sst
    ## F2S4_170908_002_D01        Sst
    ## F2S4_170908_002_E01        Sst
    ## F2S4_170908_002_F01        Sst
    ## F2S4_170908_002_G01        Sst
    ## F2S4_170908_002_H01        Sst
    ## F2S4_170908_003_A01        Sst
    ## F2S4_170908_003_B01        Sst
    ## F2S4_170908_003_C01        Sst
    ## F2S4_170908_003_D01        Sst
    ## F2S4_170908_003_E01        Sst
    ## F2S4_170908_003_G01        Sst
    ## F2S4_170908_003_H01        Sst
    ## F2S4_170908_004_A01        Sst
    ## F2S4_170908_004_B01        Sst
    ## F2S4_170908_004_C01        Sst
    ## F2S4_170908_004_D01        Sst
    ## F2S4_170908_004_E01        Sst
    ## F2S4_170908_004_F01        Sst
    ## F2S4_170908_004_G01        Sst
    ## F2S4_170908_004_H01        Sst
    ## F2S4_170908_005_A01        Sst
    ## F2S4_170908_005_B01        Sst
    ## F2S4_170908_005_C01        Sst
    ## F2S4_170908_005_D01        Sst
    ## F2S4_170908_005_E01        Sst
    ## F2S4_170908_005_F01        Sst
    ## F2S4_170908_005_G01        Sst
    ## F2S4_170908_005_H01        Sst
    ## F2S4_170908_006_A01        Sst
    ## F2S4_170908_006_B01        Sst
    ## F2S4_170908_006_C01        Sst
    ## F2S4_170908_006_D01        Sst
    ## F2S4_170908_006_E01        Sst
    ## F2S4_170908_006_F01        Sst
    ## F2S4_170908_006_G01        Sst
    ## F2S4_170908_006_H01        Sst
    ## F2S4_170908_007_A01        Sst
    ## F2S4_170908_007_B01        Sst
    ## F2S4_170908_007_C01        Sst
    ## F2S4_170908_007_D01        Sst
    ## F2S4_170908_007_F01        Sst
    ## F2S4_170908_007_G01        Sst
    ## F2S4_170908_007_H01        Sst
    ## F2S4_170908_008_A01        Sst
    ## F2S4_170908_008_B01        Sst
    ## F2S4_170908_008_C01        Sst
    ## F2S4_170908_008_D01        Sst
    ## F2S4_170908_008_E01        Sst
    ## F2S4_170908_008_F01        Sst
    ## F2S4_170908_008_G01        Sst
    ## F2S4_170908_008_H01        Sst
    ## F2S4_170908_009_A01        Sst
    ## F2S4_170908_009_B01        Sst
    ## F2S4_170908_009_C01        Sst
    ## F2S4_170908_009_D01        Sst
    ## F2S4_170908_009_E01        Sst
    ## F2S4_170908_009_F01        Sst
    ## F2S4_170908_009_G01        Sst
    ## F2S4_170908_009_H01        Sst
    ## F2S4_170908_010_A01        Sst
    ## F2S4_170908_010_B01        Sst
    ## F2S4_170908_010_C01        Sst
    ## F2S4_170908_010_D01        Sst
    ## F2S4_170908_010_E01        Sst
    ## F2S4_170908_010_F01        Sst
    ## F2S4_170908_010_G01        Sst
    ## F2S4_170908_010_H01        Sst
    ## F2S4_170908_011_A01        Sst
    ## F2S4_170908_011_B01        Sst
    ## F2S4_170908_011_C01        Sst
    ## F2S4_170908_011_D01        Sst
    ## F2S4_170908_011_E01        Sst
    ## F2S4_170908_011_F01        Sst
    ## F2S4_170908_011_G01        Sst
    ## F2S4_170908_011_H01        Sst
    ## F2S4_170908_012_A01        Sst
    ## F2S4_170908_012_B01        Sst
    ## F2S4_170908_012_C01        Sst
    ## F2S4_170908_012_D01        Sst
    ## F2S4_170908_012_E01        Sst
    ## F2S4_170908_012_F01        Sst
    ## F2S4_170908_012_G01        Sst
    ## F2S4_170908_012_H01        Sst
    ## F2S4_170908_013_A01        Sst
    ## F2S4_170908_013_B01        Sst
    ## F2S4_170908_013_C01        Sst
    ## F2S4_170908_013_D01        Sst
    ## F2S4_170908_013_E01        Sst
    ## F2S4_170908_013_F01        Sst
    ## F2S4_170908_013_G01        Sst
    ## F2S4_170908_013_H01        Sst
    ## F2S4_170908_014_D01        Sst
    ## F2S4_170908_014_E01        Sst
    ## F2S4_170908_014_F01        Sst
    ## F2S4_170908_014_G01        Sst
    ## F2S4_170908_014_H01        Sst
    ## F2S4_170914_001_A01      L6 IT
    ## F2S4_170914_001_B01      L6 IT
    ## F2S4_170914_001_C01      L6 IT
    ## F2S4_170914_001_D01      L6 IT
    ## F2S4_170914_001_E01      L6 IT
    ## F2S4_170914_001_F01      L6 IT
    ## F2S4_170914_001_G01      L6 IT
    ## F2S4_170914_001_H01      L6 IT
    ## F2S4_170914_002_A01      L6 IT
    ## F2S4_170914_002_B01      L6 IT
    ## F2S4_170914_002_C01      L6 IT
    ## F2S4_170914_002_D01      L6 IT
    ## F2S4_170914_002_E01      L6 IT
    ## F2S4_170914_002_F01      L6 IT
    ## F2S4_170914_002_G01      L6 IT
    ## F2S4_170914_003_A01      L6 IT
    ## F2S4_170914_003_B01      L6 IT
    ## F2S4_170914_003_C01      L6 IT
    ## F2S4_170914_003_D01      L6 IT
    ## F2S4_170914_003_E01      L6 IT
    ## F2S4_170918_001_A01      L6 IT
    ## F2S4_170918_001_B01      L6 IT
    ## F2S4_170918_001_C01      L6 IT
    ## F2S4_170918_001_D01      L6 IT
    ## F2S4_170918_001_E01      L6 IT
    ## F2S4_170918_001_F01      L6 IT
    ## F2S4_170918_001_G01      L6 IT
    ## F2S4_170918_001_H01      L6 IT
    ## F2S4_170918_002_A01      L6 IT
    ## F2S4_170918_002_B01      L6 IT
    ## F2S4_170918_002_C01      L6 IT
    ## F2S4_170918_002_D01      L6 IT
    ## F2S4_170918_002_E01      L6 IT
    ## F2S4_170918_002_F01      L6 IT
    ## F2S4_170918_002_G01      L6 IT
    ## F2S4_170918_002_H01      L6 IT
    ## F2S4_170918_003_A01      L6 IT
    ## F2S4_171003_007_A01    L2/3 IT
    ## F2S4_171003_007_B01      Lamp5
    ## F2S4_171003_007_C01        Vip
    ## F2S4_171003_007_D01      Lamp5
    ## F2S4_171003_007_E01      Pvalb
    ## F2S4_171003_007_F01        Sst
    ## F2S4_171003_007_G01    L2/3 IT
    ## F2S4_171003_007_H01        Vip
    ## F2S4_171003_016_A01      L6 IT
    ## F2S4_171003_016_B01      Lamp5
    ## F2S4_171003_016_C01      L6 IT
    ## F2S4_171003_016_D01        Sst
    ## F2S4_171003_016_E01        L6b
    ## F2S4_171003_016_F01      L6 IT
    ## F2S4_171003_016_G01        Sst
    ## F2S4_171003_016_H01        Sst
    ## F2S4_171003_017_A01      L6 IT
    ## F2S4_171003_017_B01      L6 IT
    ## F2S4_171003_017_C01      L6 IT
    ## F2S4_171003_017_D01        Sst
    ## F2S4_171003_017_E01        Sst
    ## F2S4_171003_017_F01        Sst
    ## F2S4_171003_017_G01        Sst
    ## F2S4_171003_017_H01      L6 IT
    ## F2S4_171004_003_A01        Sst
    ## F2S4_171004_003_B01      L5 PT
    ## F2S4_171004_003_C01      L5 PT
    ## F2S4_171004_003_D01      L5 PT
    ## F2S4_171004_003_E01   Serpinf1
    ## F2S4_171004_003_F01      L5 PT
    ## F2S4_171004_003_G01      L5 PT
    ## F2S4_171004_003_H01      L5 PT
    ## F2S4_171004_004_A01      L5 PT
    ## F2S4_171004_004_B01      L5 PT
    ## F2S4_171004_004_C01         NP
    ## F2S4_171004_004_D01      L5 PT
    ## F2S4_171004_004_E01         NP
    ## F2S4_171004_004_F01      L5 PT
    ## F2S4_171004_004_G01      L5 PT
    ## F2S4_171004_004_H01        Sst
    ## F2S4_171018_001_A01      Pvalb
    ## F2S4_171018_001_B01      Pvalb
    ## F2S4_171018_001_C01      Pvalb
    ## F2S4_171018_001_D01      Pvalb
    ## F2S4_171018_001_E01      Pvalb
    ## F2S4_171018_001_F01      Pvalb
    ## F2S4_171018_001_G01      Pvalb
    ## F2S4_171018_001_H01      Pvalb
    ## F2S4_171018_002_A01      Pvalb
    ## F2S4_171018_002_B01      Pvalb
    ## F2S4_171018_002_C01      Pvalb
    ## F2S4_171018_002_D01        Sst
    ## F2S4_171018_002_E01      Pvalb
    ## F2S4_171018_002_F01      Pvalb
    ## F2S4_171018_002_G01      Pvalb
    ## F2S4_171018_002_H01      Pvalb
    ## F2S4_171018_003_A01      Pvalb
    ## F2S4_171018_003_B01      Pvalb
    ## F2S4_171018_003_C01      Pvalb
    ## F2S4_171018_003_D01      Pvalb
    ## F2S4_171018_003_E01      Pvalb
    ## F2S4_171018_003_F01      Pvalb
    ## F2S4_171027_001_A01      L5 PT
    ## F2S4_171027_001_B01      L5 IT
    ## F2S4_171027_001_C01         NP
    ## F2S4_171027_001_E01        Sst
    ## F2S4_171027_001_F01      L5 PT
    ## F2S4_171027_001_G01         NP
    ## F2S4_171027_002_A01      L6 CT
    ## F2S4_171027_002_B01      L6 CT
    ## F2S4_171027_002_D01      L5 PT
    ## F2S4_171027_002_E01      L5 IT
    ## F2S4_171027_002_F01         NP
    ## F2S4_171027_002_G01        Sst
    ## F2S4_171027_002_H01        Sst
    ## F2S4_171027_003_A01      L6 CT
    ## F2S4_171027_003_B01        Sst
    ## F2S4_171027_003_D01      L5 IT
    ## F2S4_171027_003_E01         NP
    ## F2S4_171027_003_G01      L5 PT
    ## F2S4_171027_003_H01         NP
    ## F2S4_171027_004_A01         NP
    ## F2S4_171027_004_B01         NP
    ## F2S4_171027_004_C01         NP
    ## F2S4_171027_004_D01      L6 CT
    ## F2S4_171027_004_E01      L6 IT
    ## F2S4_171027_004_F01         NP
    ## F2S4_171027_004_G01         NP
    ## F2S4_171027_004_H01      L6 IT
    ## F2S4_171027_005_A01         NP
    ## F2S4_171027_005_B01         NP
    ## F2S4_171027_005_C01         NP
    ## F2S4_171027_005_D01         NP
    ## F2S4_171027_005_E01      L6 CT
    ## F2S4_171027_005_F01         NP
    ## F2S4_171027_005_G01         NP
    ## F2S4_171027_005_H01      L6 CT
    ## F2S4_171027_006_G01         NP
    ## F2S4_171108_001_B01        Sst
    ## F2S4_171108_001_C01        Sst
    ## F2S4_171108_001_D01        Sst
    ## F2S4_171108_001_E01        Sst
    ## F2S4_171108_001_F01        Sst
    ## F2S4_171108_001_G01        Sst
    ## F2S4_171108_002_A01        Sst
    ## F2S4_171108_002_B01        Sst
    ## F2S4_171108_002_C01        Sst
    ## F2S4_171108_002_D01        Sst
    ## F2S4_171108_002_E01        Sst
    ## F2S4_171108_002_F01        Sst
    ## F2S4_171108_002_G01        Sst
    ## F2S4_171108_002_H01        Sst
    ## F2S4_171108_003_A01        Sst
    ## F2S4_171108_003_B01        Sst
    ## F2S4_171108_003_C01        Sst
    ## F2S4_171108_003_D01        Sst
    ## F2S4_171108_003_E01        Sst
    ## F2S4_171108_003_F01        Sst
    ## F2S4_171108_003_G01        Sst
    ## F2S4_171108_003_H01        Sst
    ## F2S4_171108_004_B01        Sst
    ## F2S4_171108_004_C01        Sst
    ## F2S4_171108_004_D01        Sst
    ## F2S4_171108_004_E01        Sst
    ## F2S4_171108_004_F01        Sst
    ## F2S4_171108_004_G01        Sst
    ## F2S4_171108_004_H01        Sst
    ## F2S4_171109_007_A01        Sst
    ## F2S4_171109_007_B01        Sst
    ## F2S4_171109_007_C01        Sst
    ## F2S4_171109_007_D01      Pvalb
    ## F2S4_171109_007_E01        Sst
    ## F2S4_171109_007_G01        Sst
    ## F2S4_171109_007_H01        Sst
    ## F2S4_171109_008_B01        Sst
    ## F2S4_171109_008_E01        Sst
    ## F2S4_171109_008_F01        Sst
    ## F2S4_171109_008_G01        Sst
    ## F2S4_171109_008_H01        Sst
    ## F2S4_171109_009_A01        Sst
    ## F2S4_171109_009_B01        Sst
    ## F2S4_171109_009_C01        Sst
    ## F2S4_171109_009_D01        Sst
    ## F2S4_171109_009_E01        Sst
    ## F2S4_171110_007_B01      Pvalb
    ## F2S4_171114_001_A01      Lamp5
    ## F2S4_171114_001_B01      Lamp5
    ## F2S4_171114_001_C01      Lamp5
    ## F2S4_171114_001_D01        Vip
    ## F2S4_171114_001_E01      Lamp5
    ## F2S4_171114_001_F01      Lamp5
    ## F2S4_171120_001_A01      Astro
    ## F2S4_171120_001_B01      Astro
    ## F2S4_171120_001_C01      Astro
    ## F2S4_171120_001_D01      Astro
    ## F2S4_171120_001_E01      Astro
    ## F2S4_171120_001_F01      Astro
    ## F2S4_171120_001_G01      Astro
    ## F2S4_171120_001_H01      Astro
    ## F2S4_171120_002_A01      Astro
    ## F2S4_171120_002_B01      Astro
    ## F2S4_171120_002_C01      Astro
    ## F2S4_171120_002_D01         L4
    ## F2S4_171120_002_E01      Astro
    ## F2S4_171120_002_F01      Astro
    ## F2S4_171120_002_H01      Astro
    ## F2S4_171120_003_A01      Astro
    ## F2S4_171120_003_B01      Astro
    ## F2S4_171120_003_C01      Astro
    ## F2S4_171120_003_D01      Astro
    ## F2S4_171120_003_E01      Astro
    ## F2S4_171120_003_F01      Astro
    ## F2S4_171120_003_G01      Astro
    ## F2S4_171120_003_H01      Astro
    ## F2S4_171120_004_A01      Astro
    ## F2S4_171120_004_B01      Astro
    ## F2S4_171120_004_C01      Astro
    ## F2S4_171120_004_D01      Astro
    ## F2S4_171120_004_E01      Astro
    ## F2S4_171120_004_F01      Astro
    ## F2S4_171120_004_G01      Astro
    ## F2S4_171120_004_H01      Astro
    ## F2S4_171120_005_A01      Astro
    ## F2S4_171120_005_B01      Astro
    ## F2S4_171120_005_C01      Astro
    ## F2S4_171120_005_D01      Astro
    ## F2S4_171120_005_E01      Astro
    ## F2S4_171120_005_G01      Astro
    ## F2S4_171120_005_H01      Astro
    ## F2S4_171120_006_A01      Astro
    ## F2S4_171120_006_C01      Astro
    ## F2S4_171120_006_D01      Astro
    ## F2S4_171120_006_E01      Astro
    ## F2S4_171120_006_G01      Astro
    ## F2S4_171120_006_H01      Astro
    ## F2S4_171120_009_B01      Astro
    ## F2S4_171120_009_C01      Astro
    ## F2S4_171120_009_D01      L5 IT
    ## F2S4_171120_009_E01      Astro
    ## F2S4_171120_009_F01      Astro
    ## F2S4_171120_009_G01      Astro
    ## F2S4_171120_009_H01      Astro
    ## F2S4_171120_010_A01      Astro
    ## F2S4_171120_010_B01      Astro
    ## F2S4_171120_010_C01      Astro
    ## F2S4_171120_010_D01         L4
    ## F2S4_171120_010_E01         L4
    ## F2S4_171120_010_F01      Astro
    ## F2S4_171120_010_G01         L4
    ## F2S4_171120_010_H01      Astro
    ## F2S4_171120_011_B01         L4
    ## F2S4_171120_011_C01         L4
    ## F2S4_171120_011_D01         L4
    ## F2S4_171120_011_E01      Astro
    ## F2S4_171120_011_F01      Astro
    ## F2S4_171120_011_G01      Astro
    ## F2S4_171120_011_H01      Astro
    ## F2S4_171120_012_A01         L4
    ## F2S4_171120_012_B01         L4
    ## F2S4_171120_012_C01      Astro
    ## F2S4_171120_012_D01      Astro
    ## F2S4_171120_012_E01      Astro
    ## F2S4_171120_012_F01         L4
    ## F2S4_171120_012_G01      Astro
    ## F2S4_171120_012_H01         L4
    ## F2S4_171120_013_A01         L4
    ## F2S4_171120_013_B01         L4
    ## F2S4_171120_013_C01         L4
    ## F2S4_171120_013_D01         L4
    ## F2S4_171120_013_F01      Astro
    ## F2S4_171120_013_G01      Astro
    ## F2S4_171120_013_H01      Astro
    ## F2S4_171120_014_B01         L4
    ## F2S4_171120_014_C01      Astro
    ## F2S4_171120_014_D01      Astro
    ## F2S4_171120_014_E01      Astro
    ## F2S4_171120_014_F01      Astro
    ## F2S4_171120_014_H01         L4
    ## F2S4_171120_015_A01      Astro
    ## F2S4_171120_015_B01      L5 IT
    ## F2S4_171120_015_C01      Astro
    ## F2S4_171120_015_D01      Astro
    ## F2S4_171120_015_E01      L5 IT
    ## F2S4_171120_015_F01      Astro
    ## F2S4_171120_015_G01      Astro
    ## F2S4_171120_015_H01      Astro
    ## F2S4_171120_016_A01      Astro
    ## F2S4_171120_016_B01      Astro
    ## F2S4_171120_016_C01      Astro
    ## F2S4_171120_016_E01      Astro
    ## F2S4_171120_016_F01      Astro
    ## F2S4_171120_016_G01      Astro
    ## F2S4_171120_016_H01      Astro
    ## F2S4_171120_017_A01      Astro
    ## F2S4_171120_017_B01         NP
    ## F2S4_171120_017_C01      Astro
    ## F2S4_171120_017_D01      Astro
    ## F2S4_171120_017_F01      L5 IT
    ## F2S4_171120_017_G01      Astro
    ## F2S4_171120_017_H01      Astro
    ## F2S4_171120_018_C01      Astro
    ## F2S4_171120_018_D01      L5 IT
    ## F2S4_171120_018_F01      Astro
    ## F2S4_171120_018_G01      L5 IT
    ## F2S4_171120_018_H01      Astro
    ## F2S4_171120_019_A01      Astro
    ## F2S4_171120_019_B01      Astro
    ## F2S4_171120_019_C01      Astro
    ## F2S4_171120_019_D01      Astro
    ## F2S4_171120_019_F01      Astro
    ## F2S4_171120_019_G01      Astro
    ## F2S4_171120_019_H01      Astro
    ## F2S4_171120_020_A01      Astro
    ## F2S4_171120_020_B01      Astro
    ## F2S4_171120_020_C01      L5 IT
    ## F2S4_171120_020_D01      Astro
    ## F2S4_171120_020_E01         NP
    ## F2S4_171120_020_F01      Astro
    ## F2S4_171120_020_G01      L5 IT
    ## F2S4_171120_020_H01      Astro
    ## F2S4_171120_023_A01      Astro
    ## F2S4_171120_023_B01      Astro
    ## F2S4_171120_023_C01      Astro
    ## F2S4_171120_023_D01      Astro
    ## F2S4_171120_023_E01      Astro
    ## F2S4_171120_023_F01      Astro
    ## F2S4_171120_023_G01      Astro
    ## F2S4_171120_023_H01      Astro
    ## F2S4_171120_024_A01      Astro
    ## F2S4_171120_024_B01      Astro
    ## F2S4_171120_024_C01      Astro
    ## F2S4_171120_024_D01      Astro
    ## F2S4_171120_024_E01      Astro
    ## F2S4_171120_024_F01      Astro
    ## F2S4_171120_024_G01      Astro
    ## F2S4_171120_024_H01      Astro
    ## F2S4_171120_025_A01      Astro
    ## F2S4_171120_025_B01      Astro
    ## F2S4_171120_025_C01      Astro
    ## F2S4_171120_025_D01      Astro
    ## F2S4_171120_025_E01      Astro
    ## F2S4_171120_025_F01      Astro
    ## F2S4_171120_025_G01      Astro
    ## F2S4_171120_026_A01      Lamp5
    ## F2S4_171120_026_B01      Astro
    ## F2S4_171120_026_C01      Astro
    ## F2S4_171120_026_D01      Astro
    ## F2S4_171120_026_E01      Astro
    ## F2S4_171120_026_F01      Astro
    ## F2S4_171120_026_G01      Astro
    ## F2S4_171120_026_H01      Astro
    ## F2S4_171120_027_A01      Astro
    ## F2S4_171120_027_B01      Astro
    ## F2S4_171120_027_C01      Astro
    ## F2S4_171120_027_D01      Astro
    ## F2S4_171120_027_E01      L6 CT
    ## F2S4_171120_027_F01      Astro
    ## F2S4_171120_027_G01      Astro
    ## F2S4_171120_027_H01      Astro
    ## F2S4_171120_028_A01      Astro
    ## F2S4_171120_028_B01      Astro
    ## F2S4_171120_028_C01      Astro
    ## F2S4_171120_028_D01      Astro
    ## F2S4_171120_028_E01      Astro
    ## F2S4_171120_028_F01      Astro
    ## F2S4_171120_028_G01      Astro
    ## F2S4_171120_028_H01      Astro
    ## F2S4_171122_007_A01      Lamp5
    ## F2S4_171122_007_B01      Lamp5
    ## F2S4_171122_007_C01      Lamp5
    ## F2S4_171122_007_D01      Lamp5
    ## F2S4_171122_007_E01      Lamp5
    ## F2S4_171122_007_F01      Lamp5
    ## F2S4_171122_007_G01      Lamp5
    ## F2S4_171122_007_H01      Lamp5
    ## F2S4_171122_008_A01      Lamp5
    ## F2S4_171122_008_B01      Lamp5
    ## F2S4_171122_008_C01      Lamp5
    ## F2S4_171122_008_D01      Lamp5
    ## F2S4_171122_008_E01      Lamp5
    ## F2S4_171122_008_F01      Lamp5
    ## F2S4_171122_008_G01      Lamp5
    ## F2S4_171122_008_H01      Lamp5
    ## F2S4_171122_009_A01      Lamp5
    ## F2S4_171122_009_B01      Lamp5
    ## F2S4_171122_009_C01      Lamp5
    ## F2S4_171122_009_D01      Lamp5
    ## F2S4_171122_009_E01      Lamp5
    ## F2S4_171122_009_F01      Lamp5
    ## F2S4_171122_009_G01      Lamp5
    ## F2S4_171122_009_H01      Lamp5
    ## F2S4_171122_010_A01        Sst
    ## F2S4_171122_010_B01        Sst
    ## F2S4_171122_010_C01      Lamp5
    ## F2S4_171122_010_D01      Lamp5
    ## F2S4_171122_010_E01      Lamp5
    ## F2S4_171122_010_F01      Lamp5
    ## F2S4_171122_010_G01      Lamp5
    ## F2S4_171122_010_H01      Lamp5
    ## F2S4_171122_011_A01      Lamp5
    ## F2S4_171122_011_D01      Lamp5
    ## F2S4_171122_011_E01      Lamp5
    ## F2S4_171122_011_F01      Pvalb
    ## F2S4_171122_011_G01      Lamp5
    ## F2S4_171122_011_H01      Lamp5
    ## F2S4_171122_012_A01      Lamp5
    ## F2S4_171122_012_B01      Lamp5
    ## F2S4_171122_012_C01      Lamp5
    ## F2S4_171122_012_D01      Lamp5
    ## F2S4_171122_012_E01      Lamp5
    ## F2S4_171122_012_F01      Lamp5
    ## F2S4_171122_012_G01      Lamp5
    ## F2S4_171122_012_H01      Lamp5
    ## F2S4_171122_013_A01      Lamp5
    ## F2S4_171122_013_B01      Lamp5
    ## F2S4_171122_013_C01      Lamp5
    ## F2S4_171122_013_D01      Lamp5
    ## F2S4_171122_013_E01      Lamp5
    ## F2S4_171122_013_F01      Lamp5
    ## F2S4_171122_013_G01      Lamp5
    ## F2S4_171122_013_H01      Lamp5
    ## F2S4_171128_001_A01        Vip
    ## F2S4_171128_001_D01        Vip
    ## F2S4_171128_001_G01        Vip
    ## F2S4_171128_001_H01        Vip
    ## F2S4_171128_002_C01      Meis2
    ## F2S4_171128_002_E01      Pvalb
    ## F2S4_171128_003_C01      Pvalb
    ## F2S4_171128_003_F01      Pvalb
    ## F2S4_171128_004_A01        Vip
    ## F2S4_171128_004_E01      Pvalb
    ## F2S4_171128_004_G01      L6 IT
    ## F2S4_171128_004_H01      Pvalb
    ## F2S4_171128_005_A01      Pvalb
    ## F2S4_171128_005_B01      Pvalb
    ## F2S4_171128_005_C01        Sst
    ## F2S4_171128_005_D01       Sncg
    ## F2S4_171128_005_H01       Sncg
    ## F2S4_171128_006_C01       Sncg
    ## F2S4_171128_009_A01      Pvalb
    ## F2S4_171128_009_B01        Sst
    ## F2S4_171128_009_C01      L6 IT
    ## F2S4_171129_001_D01        SMC
    ## F2S4_171129_001_H01      L5 IT
    ## F2S4_171129_002_C01    L2/3 IT
    ## F2S4_171129_002_D01    L2/3 IT
    ## F2S4_171129_002_E01         L4
    ## F2S4_171129_002_F01      L5 IT
    ## F2S4_171129_002_G01        SMC
    ## F2S4_171129_003_A01        SMC
    ## F2S4_171129_004_B01        SMC
    ## F2S4_171129_004_H01        SMC
    ## F2S4_171129_006_E01    L2/3 IT
    ## F2S4_171129_006_G01    L2/3 IT
    ## F2S4_171129_007_A01      L6 IT
    ## F2S4_171129_007_B01        SMC
    ## F2S4_171129_007_C01        SMC
    ## F2S4_171129_007_D01      Lamp5
    ## F2S4_171129_007_E01        Vip
    ## F2S4_171129_007_F01         L4
    ## F2S4_171129_007_G01        SMC
    ## F2S4_171129_007_H01      Astro
    ## F2S4_171129_009_B01      L6 CT
    ## F2S4_171129_009_G01      L6 IT
    ## F2S4_171129_010_A01      L5 PT
    ## F2S4_171129_010_B01      Astro
    ## F2S4_171129_010_C01      L6 CT
    ## F2S4_171201_004_A01      Pvalb
    ## F2S4_171201_004_B01      Pvalb
    ## F2S4_171201_004_C01      Pvalb
    ## F2S4_171201_004_D01      Pvalb
    ## F2S4_171201_004_E01      Pvalb
    ## F2S4_171201_004_F01      Pvalb
    ## F2S4_171201_004_G01      Pvalb
    ## F2S4_171201_004_H01      Pvalb
    ## F2S4_171201_005_A01      Pvalb
    ## F2S4_171201_005_B01      Pvalb
    ## F2S4_171201_005_C01      Pvalb
    ## F2S4_171201_005_D01      Pvalb
    ## F2S4_171201_005_E01      Pvalb
    ## F2S4_171201_005_F01      Pvalb
    ## F2S4_171201_005_G01      Pvalb
    ## F2S4_171201_005_H01      Pvalb
    ## F2S4_171201_006_A01      Pvalb
    ## F2S4_171201_006_B01      Pvalb
    ## F2S4_171201_006_C01      Pvalb
    ## F2S4_171201_006_D01      Pvalb
    ## F2S4_171201_006_E01      Pvalb
    ## F2S4_171201_006_F01      Pvalb
    ## F2S4_171201_006_G01      Pvalb
    ## F2S4_171201_006_H01      Pvalb
    ## F2S4_171201_007_A01      Pvalb
    ## F2S4_171201_007_B01      Pvalb
    ## F2S4_171201_007_C01      Pvalb
    ## F2S4_171201_007_D01      Pvalb
    ## F2S4_171201_007_E01      Pvalb
    ## F2S4_171201_007_F01      Pvalb
    ## F2S4_171201_007_G01      Pvalb
    ## F2S4_171201_007_H01      Pvalb
    ## F2S4_171201_008_A01      Pvalb
    ## F2S4_171201_008_B01      Pvalb
    ## F2S4_171201_008_C01      Pvalb
    ## F2S4_171201_008_D01      Pvalb
    ## F2S4_171201_008_E01      Pvalb
    ## F2S4_171201_008_F01      Pvalb
    ## F2S4_171201_008_G01      Pvalb
    ## F2S4_171201_008_H01      Pvalb
    ## F2S4_171201_009_A01      Pvalb
    ## F2S4_171201_009_B01      Pvalb
    ## F2S4_171201_009_C01      Pvalb
    ## F2S4_171201_009_D01      Pvalb
    ## F2S4_171201_009_E01      Pvalb
    ## F2S4_171201_009_F01      Pvalb
    ## F2S4_171201_009_H01      Pvalb
    ## F2S4_171201_010_A01      Pvalb
    ## F2S4_171201_010_B01      Pvalb
    ## F2S4_171201_010_C01      Pvalb
    ## F2S4_171201_010_D01      Pvalb
    ## F2S4_171201_010_E01      Pvalb
    ## F2S4_171201_010_F01      Pvalb
    ## F2S4_171201_010_G01      Pvalb
    ## F2S4_171201_010_H01      Pvalb
    ## F2S4_171201_011_G01      Pvalb
    ## F2S4_171201_011_H01      Pvalb
    ## F2S4_171215_003_A01    L2/3 IT
    ## F2S4_171215_003_B01    L2/3 IT
    ## F2S4_171215_003_C01    L2/3 IT
    ## F2S4_171215_003_D01    L2/3 IT
    ## F2S4_171215_004_A01    L2/3 IT
    ## F2S4_171215_004_B01    L2/3 IT
    ## F2S4_171215_004_C01    L2/3 IT
    ## F2S4_171215_004_D01    L2/3 IT
    ## F2S4_171215_004_E01    L2/3 IT
    ## F2S4_171215_004_F01    L2/3 IT
    ## F2S4_171215_004_G01    L2/3 IT
    ## F2S4_171215_004_H01    L2/3 IT
    ## F2S4_171215_005_A01    L2/3 IT
    ## F2S4_171215_005_B01    L2/3 IT
    ## F2S4_171215_005_C01    L2/3 IT
    ## F2S4_171215_005_D01    L2/3 IT
    ## F2S4_171215_005_E01    L2/3 IT
    ## F2S4_171215_005_G01    L2/3 IT
    ## F2S4_171215_005_H01    L2/3 IT
    ## F2S4_171215_006_A01    L2/3 IT
    ## F2S4_171215_006_B01    L2/3 IT
    ## F2S4_171215_006_C01    L2/3 IT
    ## F2S4_171215_006_D01    L2/3 IT
    ## F2S4_171215_006_E01    L2/3 IT
    ## F2S4_171215_006_F01    L2/3 IT
    ## F2S4_171215_006_G01    L2/3 IT
    ## F2S4_171215_006_H01    L2/3 IT
    ## F2S4_171215_007_E01    L2/3 IT
    ## F2S4_171215_007_F01    L2/3 IT
    ## F2S4_171215_007_G01    L2/3 IT
    ## F2S4_171215_007_H01    L2/3 IT
    ## F2S4_180110_004_A01      L5 IT
    ## F2S4_180110_004_B01      L5 PT
    ## F2S4_180110_004_C01      L5 IT
    ## F2S4_180110_004_D01      L5 PT
    ## F2S4_180110_004_E01      L5 IT
    ## F2S4_180110_004_F01      L5 IT
    ## F2S4_180110_004_G01      L6 IT
    ## F2S4_180110_004_H01      L6 IT
    ## F2S4_180110_005_A01      L5 IT
    ## F2S4_180110_005_B01      L5 IT
    ## F2S4_180110_005_C01      L5 IT
    ## F2S4_180110_005_D01        Sst
    ## F2S4_180110_005_E01      L5 PT
    ## F2S4_180110_005_F01      L5 PT
    ## F2S4_180110_005_G01      L5 IT
    ## F2S4_180110_005_H01      L6 IT
    ## F2S4_180110_006_A01      L5 IT
    ## F2S4_180110_006_B01      L5 PT
    ## F2S4_180110_006_C01      L6 IT
    ## F2S4_180110_006_D01      L5 IT
    ## F2S4_180110_006_E01      L5 IT
    ## F2S4_180110_006_F01      L5 PT
    ## F2S4_180110_006_G01      L5 PT
    ## F2S4_180110_006_H01      L5 IT
    ## F2S4_180110_010_B01      L5 IT
    ## F2S4_180110_010_C01      L5 PT
    ## F2S4_180110_010_D01      L5 IT
    ## F2S4_180110_010_E01        Sst
    ## F2S4_180110_010_F01      L5 PT
    ## F2S4_180110_010_G01      L5 IT
    ## F2S4_180110_010_H01      L5 PT
    ## F2S4_180110_011_A01      L5 PT
    ## F2S4_180110_011_B01      L5 PT
    ## F2S4_180110_011_C01      L5 IT
    ## F2S4_180110_011_D01      L5 IT
    ## F2S4_180110_011_E01      L6 IT
    ## F2S4_180110_011_F01      L5 PT
    ## F2S4_180110_011_G01      L5 PT
    ## F2S4_180110_011_H01      L5 PT
    ## F2S4_180110_012_A01      L5 PT
    ## F2S4_180110_012_B01      L5 PT
    ## F2S4_180110_012_C01      L5 PT
    ## F2S4_180110_012_D01      L5 IT
    ## F2S4_180110_012_E01      L5 IT
    ## F2S4_180110_012_F01      L5 IT
    ## F2S4_180110_012_G01      L5 IT
    ## F2S4_180110_012_H01      L5 IT
    ## F2S4_180110_041_A01        Sst
    ## F2S4_180110_041_B01      L5 IT
    ## F2S4_180110_041_D01      L5 IT
    ## F2S4_180110_041_E01      L5 IT
    ## F2S4_180110_041_G01      L5 PT
    ## F2S4_180110_042_A01      L5 IT
    ## F2S4_180110_042_B01      L5 IT
    ## F2S4_180110_042_E01      L5 PT
    ## F2S4_180110_042_F01      L5 PT
    ## F2S4_180110_042_G01      L5 PT
    ## F2S4_180110_043_A01      L5 IT
    ## F2S4_180110_043_B01      L5 PT
    ## F2S4_180110_043_C01      L5 IT
    ## F2S4_180110_043_E01      L5 PT
    ## F2S4_180110_043_F01      L5 PT
    ## F2S4_180110_043_H01      L5 PT
    ## F2S4_180110_047_A01      L5 PT
    ## F2S4_180110_047_B01      L5 IT
    ## F2S4_180110_047_D01        Vip
    ## F2S4_180110_047_E01      L5 PT
    ## F2S4_180110_047_F01      L5 PT
    ## F2S4_180110_047_G01      L5 IT
    ## F2S4_180110_047_H01      L5 PT
    ## F2S4_180110_048_A01      L5 IT
    ## F2S4_180110_048_B01      L5 PT
    ## F2S4_180110_048_C01      L5 PT
    ## F2S4_180110_048_F01      L5 PT
    ## F2S4_180110_048_G01      L5 IT
    ## F2S4_180110_048_H01      L5 PT
    ## F2S4_180110_049_A01      L5 PT
    ## F2S4_180110_049_B01      L5 PT
    ## F2S4_180110_049_C01      L5 PT
    ## F2S4_180110_049_D01      L5 IT
    ## F2S4_180110_049_E01      L5 IT
    ## F2S4_180110_049_F01      L5 PT
    ## F2S4_180110_049_G01      L5 IT
    ## F2S4_180110_053_A01      L5 PT
    ## F2S4_180110_053_B01      L5 IT
    ## F2S4_180110_053_D01      L5 IT
    ## F2S4_180110_053_E01      L5 IT
    ## F2S4_180110_053_F01      L5 PT
    ## F2S4_180110_053_G01      L5 IT
    ## F2S4_180110_053_H01      L5 PT
    ## F2S4_180110_054_A01      L5 PT
    ## F2S4_180110_054_B01      L5 IT
    ## F2S4_180110_054_C01      L5 PT
    ## F2S4_180110_054_D01      L5 PT
    ## F2S4_180110_054_E01      L5 PT
    ## F2S4_180110_054_F01      L5 PT
    ## F2S4_180110_054_G01      L5 IT
    ## F2S4_180110_054_H01      L5 IT
    ## F2S4_180110_057_A01      L5 PT
    ## F2S4_180110_057_B01      L5 PT
    ## F2S4_180110_057_C01      L5 PT
    ## F2S4_180110_057_E01      L5 IT
    ## F2S4_180110_057_F01      L5 PT
    ## F2S4_180110_057_G01      L5 IT
    ## F2S4_180110_057_H01      L5 PT
    ## F2S4_180110_058_A01      L5 PT
    ## F2S4_180110_058_B01      L5 PT
    ## F2S4_180110_058_C01      L5 IT
    ## F2S4_180110_058_D01      L5 IT
    ## F2S4_180110_058_E01        Sst
    ## F2S4_180110_058_F01      L5 IT
    ## F2S4_180110_058_G01      L5 IT
    ## F2S4_180110_058_H01      L5 IT
    ## F2S4_180110_065_A01      L5 IT
    ## F2S4_180110_065_B01      L5 IT
    ## F2S4_180110_065_C01      L5 PT
    ## F2S4_180110_065_D01      L5 IT
    ## F2S4_180110_065_E01      L5 PT
    ## F2S4_180110_065_F01      L5 IT
    ## F2S4_180110_065_H01      L5 IT
    ## F2S4_180110_066_A01      L5 IT
    ## F2S4_180110_066_B01      L5 PT
    ## F2S4_180110_066_C01      L5 PT
    ## F2S4_180110_066_D01      L5 PT
    ## F2S4_180110_066_E01      L5 IT
    ## F2S4_180110_066_F01      L5 PT
    ## F2S4_180110_066_G01      L5 IT
    ## F2S4_180110_066_H01      L5 PT
    ## F2S4_180110_067_A01      L5 PT
    ## F2S4_180110_067_B01      L5 PT
    ## F2S4_180110_067_C01        Sst
    ## F2S4_180110_067_D01      L5 IT
    ## F2S4_180110_067_E01      L5 PT
    ## F2S4_180110_067_F01      L5 IT
    ## F2S4_180110_067_G01      L5 IT
    ## F2S4_180110_067_H01      L5 PT
    ## F2S4_180110_071_A01      L5 PT
    ## F2S4_180110_071_B01      L5 PT
    ## F2S4_180110_071_C01      L5 PT
    ## F2S4_180110_071_D01      L5 PT
    ## F2S4_180110_071_E01      L5 PT
    ## F2S4_180110_071_G01      L5 IT
    ## F2S4_180110_071_H01      L5 IT
    ## F2S4_180110_072_A01      L5 IT
    ## F2S4_180110_072_B01      L5 IT
    ## F2S4_180110_072_C01      L5 IT
    ## F2S4_180110_072_D01      L5 PT
    ## F2S4_180110_072_E01      L5 IT
    ## F2S4_180110_072_F01        Sst
    ## F2S4_180110_072_H01      L5 PT
    ## F2S4_180110_073_A01      L5 IT
    ## F2S4_180110_073_B01      L5 IT
    ## F2S4_180110_073_C01      L5 IT
    ## F2S4_180110_073_D01      L5 PT
    ## F2S4_180110_073_F01      L5 PT
    ## F2S4_180110_073_G01      L5 PT
    ## F2S4_180110_073_H01      L5 IT
    ## F2S4_180122_021_A01         NP
    ## F2S4_180122_021_B01         NP
    ## F2S4_180122_021_C01         NP
    ## F2S4_180122_021_D01         NP
    ## F2S4_180122_021_E01         NP
    ## F2S4_180122_021_F01         NP
    ## F2S4_180122_021_G01         NP
    ## F2S4_180122_021_H01         NP
    ## F2S4_180122_022_B01         NP
    ## F2S4_180122_022_C01        Sst
    ## F2S4_180122_022_D01         NP
    ## F2S4_180122_022_E01         NP
    ## F2S4_180122_022_F01         NP
    ## F2S4_180122_022_G01         NP
    ## F2S4_180122_022_H01         NP
    ## F2S4_180122_023_A01         NP
    ## F2S4_180122_023_B01         NP
    ## F2S4_180122_023_C01         NP
    ## F2S4_180122_023_D01         NP
    ## F2S4_180122_023_E01        Sst
    ## F2S4_180122_023_F01         NP
    ## F2S4_180122_023_G01         NP
    ## F2S4_180122_023_H01         NP
    ## F2S4_180122_024_A01        Sst
    ## F2S4_180122_024_B01       Sncg
    ## F2S4_180122_024_C01        Sst
    ## F2S4_180122_024_D01         CR
    ## F2S4_180122_024_E01         NP
    ## F2S4_180122_024_F01        Sst
    ## F2S4_180122_024_G01         NP
    ## F2S4_180122_024_H01         NP
    ## F2S4_180122_025_A01         NP
    ## F2S4_180122_025_C01         NP
    ## F2S4_180122_025_D01         NP
    ## F2S4_180122_025_E01         NP
    ## F2S4_180122_025_F01        Sst
    ## F2S4_180122_025_H01         NP
    ## F2S4_180122_026_A01         NP
    ## F2S4_180122_026_B01      Pvalb
    ## F2S4_180122_026_C01        Vip
    ## F2S4_180122_026_D01         NP
    ## F2S4_180122_026_E01         NP
    ## F2S4_180122_026_F01      Pvalb
    ## F2S4_180122_026_G01         NP
    ## F2S4_180122_026_H01      L6 CT
    ## F2S4_180122_027_A01         NP
    ## F2S4_180122_027_C01         NP
    ## F2S4_180122_027_D01        Sst
    ## F2S4_180122_027_E01         NP
    ## F2S4_180122_027_F01         NP
    ## F2S4_180122_027_G01         NP
    ## F2S4_180122_027_H01         NP
    ## F2S4_180226_001_C01    L2/3 IT
    ## F2S4_180228_004_A01         L4
    ## F2S4_180228_004_B01         L4
    ## F2S4_180228_004_C01      L5 IT
    ## F2S4_180228_004_D01         L4
    ## F2S4_180228_004_E01      L5 IT
    ## F2S4_180228_004_F01         L4
    ## F2S4_180228_004_G01         L4
    ## F2S4_180228_004_H01         L4
    ## F2S4_180228_005_A01         L4
    ## F2S4_180228_005_B01      L5 IT
    ## F2S4_180228_005_C01      L5 IT
    ## F2S4_180228_005_D01         L4
    ## F2S4_180228_005_E01         L4
    ## F2S4_180228_005_F01         L4
    ## F2S4_180228_005_G01      L5 IT
    ## F2S4_180228_005_H01         L4
    ## F2S4_180228_006_A01         L4
    ## F2S4_180228_006_B01      L5 IT
    ## F2S4_180228_006_C01      Astro
    ## F2S4_180228_006_D01      L5 IT
    ## F2S4_180228_006_E01      L5 IT
    ## F2S4_180228_006_F01         L4
    ## F2S4_180228_006_G01         L4
    ## F2S4_180228_006_H01      L5 IT
    ## F2S4_180228_007_A01      L5 IT
    ## F2S4_180228_007_B01      L5 IT
    ## F2S4_180228_007_C01      L5 IT
    ## F2S4_180228_007_D01      L5 IT
    ## F2S4_180228_007_E01      L5 IT
    ## F2S4_180228_007_F01         L4
    ## F2S4_180228_007_G01      L5 IT
    ## F2S4_180228_007_H01      L5 IT
    ## F2S4_180228_008_A01      L5 IT
    ## F2S4_180228_008_B01      L5 IT
    ## F2S4_180228_008_C01      L5 IT
    ## F2S4_180228_008_D01      L5 IT
    ## F2S4_180228_008_E01      L5 IT
    ## F2S4_180228_008_F01      L5 IT
    ## F2S4_180228_008_G01      L5 IT
    ## F2S4_180228_008_H01      L5 IT
    ## F2S4_180228_009_A01      L5 IT
    ## F2S4_180228_009_C01      L5 IT
    ## F2S4_180228_009_D01      L5 IT
    ## F2S4_180228_009_E01      L5 IT
    ## F2S4_180228_009_F01      L5 IT
    ## F2S4_180228_009_G01      L5 IT
    ## F2S4_180228_009_H01      L5 IT
    ## F2S4_180228_016_A01         L4
    ## F2S4_180228_016_B01         L4
    ## F2S4_180228_016_C01         L4
    ## F2S4_180228_016_D01         L4
    ## F2S4_180228_016_E01         L4
    ## F2S4_180228_016_F01      L5 IT
    ## F2S4_180228_016_G01         L4
    ## F2S4_180228_016_H01         L4
    ## F2S4_180228_017_A01         L4
    ## F2S4_180228_017_B01         L4
    ## F2S4_180228_017_C01      L5 IT
    ## F2S4_180228_017_D01         L4
    ## F2S4_180228_017_E01         L4
    ## F2S4_180228_017_F01         L4
    ## F2S4_180228_017_G01         L4
    ## F2S4_180228_017_H01         L4
    ## F2S4_180228_018_A01         L4
    ## F2S4_180228_018_B01         L4
    ## F2S4_180228_018_C01         L4
    ## F2S4_180228_018_D01         L4
    ## F2S4_180228_018_E01         L4
    ## F2S4_180228_018_F01         L4
    ## F2S4_180228_018_G01         L4
    ## F2S4_180228_018_H01         L4
    ## F2S4_180228_019_A01      L5 IT
    ## F2S4_180228_019_B01      Astro
    ## F2S4_180228_019_C01      L5 IT
    ## F2S4_180228_019_D01      L5 IT
    ## F2S4_180228_019_E01      L5 IT
    ## F2S4_180228_019_F01      L5 IT
    ## F2S4_180228_019_G01      L5 IT
    ## F2S4_180228_019_H01      L5 IT
    ## F2S4_180228_020_A01      L5 IT
    ## F2S4_180228_020_B01      L5 IT
    ## F2S4_180228_020_C01      L5 IT
    ## F2S4_180228_020_D01      L5 IT
    ## F2S4_180228_020_E01      L5 IT
    ## F2S4_180228_020_F01      L5 IT
    ## F2S4_180228_020_G01      L5 IT
    ## F2S4_180228_020_H01      L5 IT
    ## F2S4_180228_021_A01      L5 IT
    ## F2S4_180228_021_B01         L4
    ## F2S4_180228_021_C01      L5 IT
    ## F2S4_180228_021_D01      L5 IT
    ## F2S4_180228_021_E01      L5 IT
    ## F2S4_180228_021_F01      L5 IT
    ## F2S4_180228_021_G01      L5 IT
    ## F2S4_180228_021_H01      L6 IT
    ## F2S4_180312_005_A01        Sst
    ## F2S4_180312_005_B01        Sst
    ## F2S4_180312_005_C01        Sst
    ## F2S4_180312_005_D01        Sst
    ## F2S4_180312_005_E01        Sst
    ## F2S4_180312_005_F01        Sst
    ## F2S4_180312_005_G01        Sst
    ## F2S4_180312_005_H01        Sst
    ## F2S4_180312_006_A01        Sst
    ## F2S4_180312_006_B01        Sst
    ## F2S4_180312_006_C01        Sst
    ## F2S4_180312_006_D01        Sst
    ## F2S4_180312_006_E01        Sst
    ## F2S4_180312_006_F01        Sst
    ## F2S4_180312_006_G01        Sst
    ## F2S4_180312_006_H01        Sst
    ## F2S4_180312_007_A01        Sst
    ## F2S4_180312_007_B01        Sst
    ## F2S4_180312_007_C01        Sst
    ## F2S4_180312_007_D01        Sst
    ## F2S4_180312_007_E01        Sst
    ## F2S4_180312_007_F01        Sst
    ## F2S4_180312_007_G01        Sst
    ## F2S4_180312_007_H01        Sst
    ## F2S4_180312_008_A01        Sst
    ## F2S4_180312_008_B01        Sst
    ## F2S4_180312_008_C01      Pvalb
    ## F2S4_180312_008_D01        Sst
    ## F2S4_180312_008_E01        Sst
    ## F2S4_180312_008_F01        Sst
    ## F2S4_180312_008_G01        Sst
    ## F2S4_180312_008_H01        Sst
    ## F2S4_180312_009_A01        Sst
    ## F2S4_180312_009_B01        Sst
    ## F2S4_180312_009_C01        Sst
    ## F2S4_180312_009_D01        Sst
    ## F2S4_180312_009_E01        Sst
    ## F2S4_180312_009_F01        Sst
    ## F2S4_180312_009_G01        Sst
    ## F2S4_180312_009_H01        Sst
    ## F2S4_180312_010_B01        Sst
    ## F2S4_180312_010_C01        Sst
    ## F2S4_180312_010_D01        Sst
    ## F2S4_180312_010_E01        Sst
    ## F2S4_180312_010_F01        Sst
    ## F2S4_180312_010_G01        Sst
    ## F2S4_180312_010_H01        Sst
    ## F2S4_180312_011_B01        Sst
    ## F2S4_180312_011_D01        Sst
    ## F2S4_180312_011_E01        Sst
    ## F2S4_180312_011_F01        Sst
    ## F2S4_180312_011_G01        Sst
    ## FYS4_171003_504_A01        Vip
    ## FYS4_171003_504_B01      Lamp5
    ## FYS4_171003_504_C01      Lamp5
    ## FYS4_171003_504_D01        Vip
    ## FYS4_171003_504_E01        Vip
    ## FYS4_171003_504_F01        Vip
    ## FYS4_171003_504_G01      Lamp5
    ## FYS4_171003_504_H01        Sst
    ## FYS4_171003_508_A01        Sst
    ## FYS4_171003_508_B01      L6 IT
    ## FYS4_171003_508_C01        Sst
    ## FYS4_171003_508_D01      L6 IT
    ## FYS4_171003_508_E01        Vip
    ## FYS4_171003_508_F01        Sst
    ## FYS4_171003_508_G01      L6 IT
    ## FYS4_171003_508_H01      L6 IT
    ## FYS4_171003_509_A01      L6 IT
    ## FYS4_171003_509_B01      L6 IT
    ## FYS4_171003_509_C01        Sst
    ## FYS4_171003_509_D01      L6 IT
    ## FYS4_171003_509_E01      Pvalb
    ## FYS4_171003_509_F01      L6 IT
    ## FYS4_171003_509_G01      L6 IT
    ## FYS4_171003_509_H01        Sst
    ## FYS4_171004_102_A01         NP
    ## FYS4_171004_102_B01      Pvalb
    ## FYS4_171004_102_C01      Pvalb
    ## FYS4_171004_102_D01      L5 PT
    ## FYS4_171004_102_E01      L5 PT
    ## FYS4_171004_102_F01        Sst
    ## FYS4_171004_102_G01        Sst
    ## FYS4_171004_102_H01        Vip
    ## FYS4_171004_103_A01      Lamp5
    ## FYS4_171004_103_B01       Sncg
    ## FYS4_171004_103_C01      Pvalb
    ## FYS4_171004_103_D01        Sst
    ## FYS4_171004_103_E01      L5 PT
    ## FYS4_171004_103_F01      L5 PT
    ## FYS4_171004_103_G01        Sst
    ## FYS4_171004_103_H01      L5 PT
    ## FYS4_171004_104_A01       Sncg
    ## FYS4_171004_104_B01      Pvalb
    ## FYS4_171004_104_C01      L5 PT
    ## FYS4_171004_104_D01        Sst
    ## FYS4_171004_104_F01      L5 PT
    ## FYS4_171004_104_G01        Sst
    ## FYS4_171004_104_H01      Pvalb

In this tutorial, we apply the ‘anchor’-based integration workflow, that
enables the probabilistic transfer of annotations from a reference to a
query set.

``` r
options(future.globals.maxSize = 8000 * 1024^2)
allen_reference <- SCTransform(allen_reference)
```

    ## Running SCTransform on assay: RNA

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 34608 by 14249

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 5000 cells

    ## Found 33 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 34608 genes

    ## Computing corrected count matrix for 34608 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 52.29345 secs

    ## Determine variable features

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## Set default assay to SCT

``` r
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
```

    ## Normalizing query using reference SCT model

    ## Performing PCA on the provided reference using 2735 features as input.

    ## Projecting cell embeddings

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 266 anchors

``` r
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = cortex[["pca"]], dims = 1:30)
```

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Predicting cell labels

``` r
cortex[["predictions"]] <- predictions.assay
rm(allen_reference)
```

Now we get prediction scores for each spot for each class. Of particular
interest in the frontal cortex region are the laminar excitatory
neurons. Here we can distinguish between distinct sequential layers of
these neuronal subtypes, for example:

![](Spatialtranscriptomics_files/figure-gfm/deconvolution%20results-1.png)<!-- -->

## Working with multiple slices

Until now we are only focusing on the anterior section of the brain, now
we load the other half of the brain and perform the same quality control
and normalization

    ## Validating object structure

    ## Updating object slots

    ## Ensuring keys are in the proper structure
    ## Ensuring keys are in the proper structure

    ## Ensuring feature names don't have underscores or pipes

    ## Updating slots in Spatial

    ## Updating slots in posterior1

    ## Validating object structure for Assay5 'Spatial'

    ## Validating object structure for VisiumV2 'posterior1'

    ## Object representation is consistent with the most current Seurat version

Then we can merge multiple slice in the same Seurat object:

Now plot the gene expression onto the merge tissue sections

![](Spatialtranscriptomics_files/figure-gfm/Gene%20expression%20onto%20merge%20sections-1.png)<!-- -->

Running the joint dimensional reduction and clustering on the underlying
RNA expression data.

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## 14:49:23 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:49:23 Read 5659 rows and found 30 numeric columns

    ## 14:49:23 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:49:23 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:49:23 Writing NN index file to temp file /var/folders/kq/_ntfn5zx1b199xp4vyd82t840000gn/T//RtmpqlyX5i/file82f7285eb92
    ## 14:49:23 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:49:24 Annoy recall = 100%
    ## 14:49:24 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:49:25 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 14:49:25 Commencing optimization for 500 epochs, with 226224 positive edges
    ## 14:49:29 Optimization finished

Now we can visualize the join dimensionality reduction and clustering
results in UMAP plot and onto the tissue slides:

![](Spatialtranscriptomics_files/figure-gfm/joint%20UMAP-1.png)<!-- -->

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.
    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](Spatialtranscriptomics_files/figure-gfm/Clustering%20onto%20merge%20sections-1.png)<!-- -->

You can check if the clustering labels here match the reference brain
region labels in [Allen Brain
atlas](http://atlas.brain-map.org/atlas?atlas=2&plate=100883900#atlas=2&plate=100883818&resolution=10.47&x=8709.999694824219&y=4040&zoom=-3)

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
    ## time zone: Europe/Amsterdam
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.4               patchwork_1.2.0          
    ## [3] ggplot2_3.5.1             stxBrain.SeuratData_0.1.2
    ## [5] SeuratData_0.2.2.9001     Seurat_5.1.0             
    ## [7] SeuratObject_5.0.2        sp_2.1-4                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          jsonlite_1.8.8             
    ##   [3] magrittr_2.0.3              spatstat.utils_3.1-0       
    ##   [5] ggbeeswarm_0.7.2            farver_2.1.2               
    ##   [7] rmarkdown_2.28              zlibbioc_1.50.0            
    ##   [9] vctrs_0.6.5                 ROCR_1.0-11                
    ##  [11] DelayedMatrixStats_1.26.0   spatstat.explore_3.3-2     
    ##  [13] S4Arrays_1.4.1              htmltools_0.5.8.1          
    ##  [15] SparseArray_1.4.8           sctransform_0.4.1          
    ##  [17] parallelly_1.38.0           KernSmooth_2.23-24         
    ##  [19] htmlwidgets_1.6.4           ica_1.0-3                  
    ##  [21] plyr_1.8.9                  plotly_4.10.4              
    ##  [23] zoo_1.8-12                  igraph_2.0.3               
    ##  [25] mime_0.12                   lifecycle_1.0.4            
    ##  [27] pkgconfig_2.0.3             Matrix_1.7-0               
    ##  [29] R6_2.5.1                    fastmap_1.2.0              
    ##  [31] GenomeInfoDbData_1.2.12     MatrixGenerics_1.16.0      
    ##  [33] fitdistrplus_1.2-1          future_1.34.0              
    ##  [35] shiny_1.9.1                 digest_0.6.37              
    ##  [37] colorspace_2.1-1            S4Vectors_0.42.1           
    ##  [39] tensor_1.5                  RSpectra_0.16-2            
    ##  [41] irlba_2.3.5.1               GenomicRanges_1.56.1       
    ##  [43] labeling_0.4.3              progressr_0.14.0           
    ##  [45] fansi_1.0.6                 spatstat.sparse_3.1-0      
    ##  [47] httr_1.4.7                  polyclip_1.10-7            
    ##  [49] abind_1.4-8                 compiler_4.4.1             
    ##  [51] withr_3.0.1                 fastDummies_1.7.4          
    ##  [53] highr_0.11                  MASS_7.3-61                
    ##  [55] DelayedArray_0.30.1         rappdirs_0.3.3             
    ##  [57] tools_4.4.1                 vipor_0.4.7                
    ##  [59] lmtest_0.9-40               beeswarm_0.4.0             
    ##  [61] httpuv_1.6.15               future.apply_1.11.2        
    ##  [63] goftest_1.2-3               glmGamPoi_1.16.0           
    ##  [65] glue_1.7.0                  nlme_3.1-166               
    ##  [67] promises_1.3.0              grid_4.4.1                 
    ##  [69] Rtsne_0.17                  cluster_2.1.6              
    ##  [71] reshape2_1.4.4              generics_0.1.3             
    ##  [73] gtable_0.3.5                spatstat.data_3.1-2        
    ##  [75] tidyr_1.3.1                 data.table_1.16.0          
    ##  [77] XVector_0.44.0              utf8_1.2.4                 
    ##  [79] BiocGenerics_0.50.0         spatstat.geom_3.3-2        
    ##  [81] RcppAnnoy_0.0.22            ggrepel_0.9.6              
    ##  [83] RANN_2.6.2                  pillar_1.9.0               
    ##  [85] stringr_1.5.1               limma_3.60.4               
    ##  [87] spam_2.10-0                 RcppHNSW_0.6.0             
    ##  [89] later_1.3.2                 splines_4.4.1              
    ##  [91] lattice_0.22-6              survival_3.7-0             
    ##  [93] deldir_2.0-4                tidyselect_1.2.1           
    ##  [95] miniUI_0.1.1.1              pbapply_1.7-2              
    ##  [97] knitr_1.48                  gridExtra_2.3              
    ##  [99] IRanges_2.38.1              SummarizedExperiment_1.34.0
    ## [101] scattermore_1.2             stats4_4.4.1               
    ## [103] xfun_0.47                   Biobase_2.64.0             
    ## [105] statmod_1.5.0               matrixStats_1.4.1          
    ## [107] UCSC.utils_1.0.0            stringi_1.8.4              
    ## [109] lazyeval_0.2.2              yaml_2.3.10                
    ## [111] evaluate_0.24.0             codetools_0.2-20           
    ## [113] tibble_3.2.1                cli_3.6.3                  
    ## [115] uwot_0.2.2                  xtable_1.8-4               
    ## [117] reticulate_1.39.0           munsell_0.5.1              
    ## [119] GenomeInfoDb_1.40.1         Rcpp_1.0.13                
    ## [121] globals_0.16.3              spatstat.random_3.3-1      
    ## [123] png_0.1-8                   ggrastr_1.0.2              
    ## [125] spatstat.univar_3.0-1       parallel_4.4.1             
    ## [127] presto_1.0.0                dotCall64_1.1-1            
    ## [129] sparseMatrixStats_1.16.0    listenv_0.9.1              
    ## [131] viridisLite_0.4.2           scales_1.3.0               
    ## [133] ggridges_0.5.6              leiden_0.4.3.1             
    ## [135] purrr_1.0.2                 crayon_1.5.3               
    ## [137] rlang_1.1.4                 cowplot_1.1.3
