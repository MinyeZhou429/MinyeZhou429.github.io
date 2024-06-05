---
layout: page
title: Single Cell Analysis of Breast cancer dataset
description: BMI710 Individual project
img: /assets/img/bc.jpg
importance: 3
category: work
---

**Title**:

Single cell analysis for breast cancer dataset

**Data type**:

Single cell sequencing transcriptomic data

**Summary of Dataset**:

The dataset used in this project is Azizi2018_BreastCancer dataset. Using single-cell RNA sequencing, this study analyzed 45,000 immune cells from eight breast carcinomas as well as matched normal breast tissue, blood, and lymph nodes.

**Method**:

Quality Check and Batch Effect Correction:

The percentage of mitochondrial genes (percent_mito) was calculated as an indicator of cell status, with higher values often suggesting increased cell stress or death. Expression percentages for key markers, EpCAM and PTPRC, were also computed to identify potential contaminants. Additionally, the logarithm of the ratio of the number of detected genes per UMI was derived to assess the overall transcription activity state of cells. Violin plots and density plots were generated for these metrics, along with 'nFeature_RNA' and 'nCount_RNA’, which measure the number of unique genes and the total number of molecules detected within a cell respectively. 

Based on the figure below, one data issue is that many cells have low numbers of genes per cell. Based on the violin plot and density plot for nFeature_RNA, we can see that many cells are centered at 500 features per UMI, which is the common lower bound. This indicates low RNA content, which can be indicative of dying cells. When we remove low-quality cells, half of the cells are expected to be removed based on the distribution. Based on the violin plot for PTPRC percentage, most cells have PTPRC percentage less than 2 with one outlier, which is removes later. Most cells have mitochondrial genes percentage less than 20, cells with mitochondrial genes percentage more than 20 are removed because they are potential dead cells. 


{% include figure.liquid path="assets/img/sc1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Considering the label of replicates does not have meaning, in order to diagnose and correct sample level batch effect, I combined patient ID and replicate ID so that each replicate from each patient represents one sample.  Based on the UMAP(Figure 2a) colored by sample, we can see that the cells are not well-mixed. This indicates sample may create a batch effect. Then, harmony is used to correct sample level batch-effect using variable ‘combined_id’, which is the combination of replicate ID and patient ID. From UMAP generated on batch-corrected PCs(Figure 2b), we can see the cells are mixed together, suggesting less batch effect. 

{% include figure.liquid path="assets/img/sc2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

After processing the data, I performed PCA on the dataset to reduce the complexity of the data while retaining the most significant features that explain the most variance among cells. These PCs are corrected by Harmony mentioned in the previous section about batch effect correction. Then, an elbow plot is generated to visualize the percentage of variance explained by each PC which is crucial for selecting a suitable cutoff to capture the most informative components without including too much noise. Based on the elbow plot(Figure 3), the first 15 PCs explained most variances in this dataset. Finally, I ran UMAP for non-linear dimensionality reduction using the first 15 principal components. UMAP projects the high-dimensional data into a lower-dimensional space, facilitating the visualization of the data. Specifically, UMAP plots are used for diagnosing batch effect and validating the effect of batch effect correction.

{% include figure.liquid path="assets/img/sc3.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

The original paper reported the main immune cell types including monocytes, macrophages, dendritic cells, T cells, B cells, mast cells, and neutrophils. 

First, I applied FindNeighbors() using batch-corrected PC to find 20 nearest neighbors of each cell. This is crucial for the clustering algorithm. Then, I repeatedly run the clustering function at different resolution parameters(0.02, 0.05, 0.2, 0.5, 0.8, 1). Some resolutions produce too many clusters while some produce fewer clusters. When resolution is set to 0.2, there are 8 clusters produced. This seems to be an optimal clustering result. The reason is that we expect 7 clusters, allowing some cell subtype, 8 seems to be a reasonable number. Also, at 0.2 resolution, clusters are well separated while not over-separated into many small clusters(Figure 5). Then, I use FindAllMarkers() to identify differentially expressed genes across clusters. The log fold change threshold is set to be 0.5. Wilcoxon rank sum test is used to compare the median expression levels of a gene between two clusters. only.pos is set to be TRUE to include genes that are positively differentially expressed. Wilcoxon test is used here because it’s flexible and it does not assume normality. From the list of differentially expressed genes, I filter to find the top 10 markers (based on log fold change) for each cluster. Based on Figure 4, we can confirm the difference in the expression of these genes across clusters. They are used as markers for cell type annotation. 

{% include figure.liquid path="assets/img/sc4.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
{% include figure.liquid path="assets/img/sc5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

The main resource used for cell type annotation is CellMarker2.0.  The cell types are determined by expression level and cell marker significance based on CellMarker2.0. For each cluster, at least two genes are used for cell type annotation. There are 8 cell types identified in this dataset: T cells, Monocytes, NK cells, B cells, Fibroblast, Dendritic cells, Mast cells, and Plasmacytoid dendritic cells. Compared to the cell type identified by the original paper, neutrophils are missing potentially due to different cluster resolution and marker genes used for each cluster. The result is summarized below.




