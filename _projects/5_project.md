---
layout: page
title: From Genes to Proteins-A Comprehensive RNA Sequencing Approach to Studying Alzheimer's Disease Pathology
description: BST281 Final group project
img: assets/img/brain.jpg
importance: 3
category: work
---

**Title**:

From Genes to Proteins-A Comprehensive RNA Sequencing Approach to Studying Alzheimer's Disease Pathology

**Members and Contributors**:

Lucy Tian, Minye Zhou, Michelle Dai, Evelyn(Sizhu) Jiang, Yushan Xu


**Data type**:

Bulk RNA-Seq data

**Summary**:

Complate workflow of RNA sequencing analysis. It starts from alignment to Gene Set Enrichment Analysis.

**Introduction and Rationale**:

Alzheimer's disease (AD) is a complex neurodegenerative disorder marked by the abnormal buildup of amyloid-β and neurofibrillary tangles in the brain, leading to cognitive decline and memory loss. Extensive research, including our study, leverages RNA sequencing to explore AD’s genetic underpinnings, revealing key genes like APP, PSEN1, and PSEN2, and risk factors such as the APOEε4 allele. This RNA-seq analysis of AD patients versus controls identifies differentially expressed genes and integrates clustering, protein-protein interactions, and pathway analyses to offer a comprehensive view of the molecular dynamics at play. 

**Methods**:

Data summary:

Our data is downloaded from the Rush Alzheimer’s Disease (AD) Project in the ENCODE consortium. To investigate comprehensively the molecular mechanisms underlying AD pathology at both the transcriptional and protein levels, we selected a dataset containing information on total RNA. We manually selected five control samples and five AD cases with sequencing data from their dorsolateral prefrontal cortex, matched by gender and sex. to reduce the potential influence of covariates and boost our confidence in our downstream analysis results.


Data Preprocessing Pipeline:

1. STAR for Alignment
2. RSEM for Quantification
3. FastQC for Quality Check

Differential expression analysis:

1. DESeq2
2. Limma

Protein protein interaction analysis:

STRING + Cytoscape

Functional analysis:

1. PASCAL to calculate pathway scores based on AD GWAS summary statistics
2. Enrichr to perform enrichment analysis on DE genes using KEGG pathways
3. Prerank perform enrichment analysis on ranked DE gene list using GO


**Results**:

Quality control:

The quality assessment of our RNA-seq data using FastQC and RSeQC tools indicates high-quality sequence data suitable for further analysis. FastQC analysis shows high average quality scores across all bases and minimal adaptor content, suggesting no significant degradation or contamination. Similarly, RSeQC's tin.py shows all samples have a median TIN score above 80, indicating excellent RNA integrity and consistent gene coverage. Alignment results further validate the data quality, with unique mapping rates around 85% and total mapping rates above 99%, contrasting sharply with previous lower rates using different methodologies. These findings collectively affirm the suitability of our data for detailed downstream genomic analysis.


{% include figure.liquid path="assets/img/az_qc.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
{% include figure.liquid path="assets/img/az_qc2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}


Differential Expression Analysis:

In our study using DESeq2 and Limma for differential expression analysis, we identified 71 genes with notable expression differences between AD and control groups. DESeq2 revealed 9 upregulated and 62 downregulated genes, while Limma identified 9 genes with more modest expression changes. The overlap between these methods underscores their reliability. Notably, genes such as ADAMTS1 and CYBA were among the differentially expressed, suggesting roles in neuroplasticity and immune responses in AD, respectively. The predominance of downregulated genes may reflect cellular adaptations to metabolic dysfunction in AD. The variance in findings between DESeq2 and Limma likely stems from their different statistical frameworks and assumptions, with DESeq2 using a negative binomial model suitable for RNA-seq data and Limma assuming a normal distribution, typically more conservative. These results provide a solid basis for further analyses to explore the molecular mechanisms and pathways involved in AD.

{% include figure.liquid path="assets/img/az_degs.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
{% include figure.liquid path="assets/img/az_degs2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Clustering:

Our clustering analysis using Hierarchical Clustering, K-Means, and Spectral Clustering assessed the quality of clusters of differentially expressed genes (DEGs) identified by DESeq2, with Hierarchical Clustering providing the best results based on silhouette scores. This measure evaluates the fit of each data point within its cluster relative to other clusters, where a higher score indicates better clustering. Hierarchical Clustering grouped DEGs into four distinct clusters, with the majority of up-regulated genes in the largest cluster. 



PPI analysis:


PPINs in our study demonstrated scale-free and "small-world" properties, typical of biological networks. Most nodes had few connections, while a few hub nodes were highly interconnected, influencing fundamental biological processes like synaptic activities and protein localization, particularly relevant to AD. We identified key hub genes such as TYROBP and CD163, which are crucial in AD pathogenesis, influencing microglial activity and inflammatory responses.

Our PPIN analysis revealed a community structure with each cluster associated with specific biological functions. For instance, one cluster was linked to smooth muscle contraction and another to the unfolded protein response, both contributing to neuroinflammation and immunity. These clusters suggest distinct pathways that could be exploited for therapeutic targets.

Further, using the guilt-by-association principle, we predicted that certain genes involved in the cell cycle could be implicated in AD, highlighting the importance of further validation studies. This analysis underscores the potential of PPINs to identify diagnostic biomarkers and therapeutic targets for AD, as well as providing a deeper understanding of the disease's molecular basis.

{% include figure.liquid path="assets/img/az_degs2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
