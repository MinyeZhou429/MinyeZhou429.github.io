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

Clustering and Cell Type Annotation:

The original paper reported the main immune cell types including monocytes, macrophages, dendritic cells, T cells, B cells, mast cells, and neutrophils. 

First, I applied FindNeighbors() using batch-corrected PC to find 20 nearest neighbors of each cell. This is crucial for the clustering algorithm. Then, I repeatedly run the clustering function at different resolution parameters(0.02, 0.05, 0.2, 0.5, 0.8, 1). Some resolutions produce too many clusters while some produce fewer clusters. When resolution is set to 0.2, there are 8 clusters produced. This seems to be an optimal clustering result. The reason is that we expect 7 clusters, allowing some cell subtype, 8 seems to be a reasonable number. Also, at 0.2 resolution, clusters are well separated while not over-separated into many small clusters(Figure 5). Then, I use FindAllMarkers() to identify differentially expressed genes across clusters. The log fold change threshold is set to be 0.5. Wilcoxon rank sum test is used to compare the median expression levels of a gene between two clusters. only.pos is set to be TRUE to include genes that are positively differentially expressed. Wilcoxon test is used here because it’s flexible and it does not assume normality. From the list of differentially expressed genes, I filter to find the top 10 markers (based on log fold change) for each cluster. Based on Figure 4, we can confirm the difference in the expression of these genes across clusters. They are used as markers for cell type annotation. 

{% include figure.liquid path="assets/img/sc4.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
{% include figure.liquid path="assets/img/sc5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

The main resource used for cell type annotation is CellMarker2.0.  The cell types are determined by expression level and cell marker significance based on CellMarker2.0. For each cluster, at least two genes are used for cell type annotation. There are 8 cell types identified in this dataset: T cells, Monocytes, NK cells, B cells, Fibroblast, Dendritic cells, Mast cells, and Plasmacytoid dendritic cells. Compared to the cell type identified by the original paper, neutrophils are missing potentially due to different cluster resolution and marker genes used for each cluster. The result is summarized below.

{% include figure.liquid path="assets/img/sc_table.jpg" title="example image" class="img-fluid rounded z-depth-1" %}


Differential Abundance Analysis for each cell type:

Differential abundance analysis is done on eight cell type. For each cell type, I used boxplot to plot their proportion stratified by tissues using ggplot2. Then, I visually inspected the distribution to identify two groups with different distributions and compared the mean of each group. Based Figure 6, T cells, B cells, Dendritic cells,  seem to be potentially differentially abundant between: 1. normal and tumor tissues; 2. blood and tumor tissue. For monocytes, NK cells, mast cells, cells seem to be potentially differentially abundant between normal and tumor tissue. For Plasmacytoid dendritic cells, cells seem to be differentially abundant between blood and tumor tissues. No trend is observed for fibroblasts because there is a small number of fibroblasts.

Solely relying on visual inspection is not enough. Wilcoxon rank sum exact test is done to compare the abundance of cells of each cell type between the two tissues systematically and statistically. None of p-values is significant for these comparisons. This suggests that there is no differentially abundant cell type. This is potentially due to the choice of resolution in clustering. Since the number of clusters is relatively small, subtle but potentially important differences might be overlooked. 

{% include figure.liquid path="assets/img/sc6.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Differential expression analysis
For top three abundant cell types, I first subsetted the data to include only the cells belonging to that type. Then, I performed differential expression analysis for each cell type across different tissue types using Seurat Findmakers(). In findmarkers(), I specified ident.1 to be tumor tissue and ident.2 to be normal tissue. Only. pos is set to be TRUE because we only care about up-regulated genes. The test used here is DESeq2 as DESeq2 has rigorous statistical methods and normalization which is suitable for capturing small changes between cells in different tissues. A pseudocount 1 is added when conducting DE analysis on monocytes to avoid bugs. The differentially expressed genes are sorted based on log2FC. Positive and large log2FC means genes are upregulated with a large difference between two tissues. Then ranked gene list feeds into fgsea for pathway analysis. Ontology gene sets is used to investigate what pathways are associated with differentially expressed genes between tumor and normal samples for each cell type. Enriched pathways are sorted according to NES.

For T cell, HSPA6, SPP1, and FN1 are the top three over-expressed genes in tumor tissue. Previous study has shown that HSPA6 gene was significantly upregulated and was found to suppress the growth, migration, and invasion of triple-negative breast cancer cells. Higher levels of HSPA6 expression were also associated with longer overall survival in breast cancer patients, highlighting its potential as a tumor suppressor[1]. Increased expression of the SPP1 gene was linked to a higher risk of breast cancer recurrence[2].  Other studies have also reported that FN1 was found to be upregulated in breast cancer tissues compared to normal tissues. Higher levels of FN1 mRNA were associated with poorer clinical outcomes and were effective in predicting the survival of BRCA patients. Additionally, Cox regression analysis identified FN1 as an independent prognostic factor for overall survival in BRCA patients[3]. The top three enriched pathways are: GOBP_COMPLEMENT_ACTIVATION, GOCC_BLOOD_MICROPARTICLE, GOCC_IMMUNOGLOBULIN_COMPLEX. Complement activation refers to any activity initiating the activation of any phase within the complement cascade enables the direct elimination of microbes, clearance of immune complexes, and regulation of additional immune functions. In the context of breast cancer, upregulation of complement activation in T cells within the tumor microenvironment suggests an active immune response against tumor cells. Blood microparticle are small membrane-bound vesicles released from various cell types. It has been recognized that there's an elevation in the presence of cancer in the circulation of platelet-derived microparticles. These microparticles have been suggested to support cancer growth through various mechanisms, such as shielding circulating tumor cells, thus aiding in immune evasion[4]. Immunoglobulin complexes are proteins typically consisting of two identical immunoglobulin heavy chains and two identical immunoglobulin light chains, interconnected by disulfide bonds and occasionally associated with other proteins. Immunoglobulin complexes can be found integrated into the plasma membrane, located in extracellular regions such as mucosal surfaces or various tissues, or circulating within the bloodstream or lymphatic system. In the context of breast cancer, the enrichment of immunoglobulin complex-related pathways in T cells may indicate their involvement in recognizing and responding to tumor-specific antigens. 

For Monocytes, CHI3L1, SPP1, and CCL18 are the top three over-expressed genes in tumor tissue. CHI3L1 promotes tumor growth by influencing the movement and infiltration of cancerous cells by altering the surroundings within the tumor[5]. SPP1 is also differentially expressed in T cells which is discussed in the above session. Study has shown that TAMs produce significant amounts of CCL18 in breast cancer. They found that elevated CCL18 expression in blood or the cancer stroma correlates with metastasis and reduced patient survival[6]. The top three enriched pathways are: HP_MALAR_RASH, GOBP_SEMI_LUNAR_VALVE_DEVELOPMENT, and HP_RAYNAUD_PHENOMENON. 
The first pathway is associated with malar rash, a facial rash that is a characteristic feature of systemic lupus erythematosus. Although it seems less relevant to breast cancer, rash is one of the earliest and most prevalent indicators of inflammatory breast cancer. GOBP_SEMI_LUNAR_VALVE_DEVELOPMENT is associated with semi lunar valve formation. This process is irrelevant to breast cancer. The third pathway is associated with reduced blood flow. In breast cancer, vascular abnormalities and alterations in blood flow dynamics might play a role in tumor progression and metastasis. Overall, these pathways are less relevant potentially due to the pseudocount. 

For NK cells, SPP1, APOE, and IGLC3 are the top three over-expressed genes in tumor tissue. APOE is shown to be associated with breast cancer as the expression of apoA-I or apoE enhances the proliferation, migration, and tumor growth of MCF-7 cells[7].  IGLC3 facilitates antigen binding activity and interaction with immunoglobulin receptors. They participate in various processes, including the activation of immune responses. In the context of breast cancer, upregulation of IGLC3 might activate immune responses against cancer cells. The top three enriched pathways are: GOBP_CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS, GOBP_NEUTRAL_LIPID_BIOSYNTHETIC_PROCESS, and GOBP_POSITIVE_REGULATION_OF_CHOLESTEROL_EFFLUX. The first pathway suggests suggests that NK cells might respond to stimuli from lipoprotein particles. Liproprotein is shown to have higher expression in breast cancer samples, which could potentially trigger the cellular response[8]. The second pathway is associated with the process of neutral lipids formation. It has been previously shown that lipid metabolism is an indicator of breast cancer[9]. GOBP_POSITIVE_REGULATION_OF_CHOLESTEROL_EFFLUX means processes that elevate cholesterol efflux level. Cholesterol level is associated with breast cancer and should be controlled for better outcome[10]. 

Collectively, these findings suggest that different genes and pathways contribute to breast cancers across different cell types. 


