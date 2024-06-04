---
layout: page
title: Gene regulatory network
description: with background image
img: assets/img/subnetwork.jpg
importance: 1
category: work
related_publications: false
---

Title: Integrating RNA-Seq and ATAC-Seq Data to Infer Gene Regulatory Networks

**Summary**: Using dual-pipeline approach to build gene regulatory network

**Data Type**: Bulk-RNA Transcriptomic data and ATAC-Seq data

**Tools**: Scenic+ pipeline, HOMER, T-Gene

**Introduction and rationale**:
Gene regulation is a process controlling gene’s expression through many mechanisms. Gene regulation is important for proper cellular function, development, and adaption to different environments. Genes can be regulated at different stages. This project focuses on transcription factor gene regulatory events. How can genes be regulated through transcription factors? How can we use ATAC-Seq and RNA-seq data to infer gene regulatory network? 

**State of art**:
Historically, if we want to build gene regulatory network, we can only rely on experimental data.  With the advancements in sequencing techniques and computational biology, we can obtain gene expression data. So we can build GRN based on co-expression relationships. Initially, the networks are undirected. Undirected networks can provide very limited information as we don’t know who regulates who. To improve this, techniques like GENIE3 and GRNBoost2 help by focusing on known transcription factors and their target genes. These algorithms make the relashionships directed. However, these methods still have their downsides. The main reason is that transcriptomics data do not directly capture many underlying regulatory mechanisms, such chromatin accesibility. Now we have chip-seq and atac-seq, so new methods got developed to incorporate new information. 

**Method and pipeline**:
{% include figure.liquid loading="eager" path="assets/img/scenicplus_pipeline.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Data descriptions**:
Source: The Cancer Genome Atlas(TCGA)

Cohort: Breast cancer cohort
{% include figure.liquid loading="eager" path="assets/img/data.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Results**:
SCENIC+ builds a grn with 322 tfs. Each tf has around 30 target genes. I chose five well studies tfs in breast cancer and draw this subnetwork. Red indicates transcription factors, blue indicates target genes. In this network, we can see there are potential co-regulate events. For example, GATA3 and ESR1 both point to TBC1D9. GATA3 frameshift mutation promotes tumor growth in human luminal breast cancer cells. Mutation in ESR1 promotes ET-resistant and metastatic. The function of TBC1D9 in breast cancer is less characterized. Connecting the dots, we can hypothesize that TBC1D9 could potentially play a specific role in mediating responses that involve both the immune system's regulatory mechanisms and hormone-driven growth. 

{% include figure.liquid loading="eager" path="assets/img/subnetwork.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Among over 10 thousand TF-gene pairs. Most of TF-gene are positive correlated in term of expression. On the left is a scatter plot showing the correlation coefficients of each tf-gene pair. X-axis shows correlation coefficient, y-axis shows log p-value. Based on this plot, we can see that many tf-gene pairs have a significant positive correlation. On the right, I plot the number of positive correlated pairs versus negative correlated pair. Among all significant correlations, the number of positive correlared tf-gene pair is doubled compared to negative correlated pairs. These means many tf play a role in up-regulating target genes in this breast cancer cohort. This makes sense because based on previous research we know that some tfs upregulate target gene to promote tumor growth.

{% include figure.liquid loading="eager" path="assets/img/upregulation.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Differential expression analysis of transcription factors across different subtype may show some potential different activities across subtypes. Here is a plot showing the differential expression of tfs in basal subtype. For example, FOXC2 is one of the de’ed tf in basal subtype. Literature has shown that the expression of FOXC2 is significantly correlated with the highly aggressive basal-like subtype breast cancers. 

{% include figure.liquid loading="eager" path="assets/img/degs1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

The following plot showing the DE’ed tfs in luminal-a subtype. ESR2 belongs estrogen receptor family. This aligns with our knowledge that luminal b subtype are estrogen receptor-positive. CEBPA and ZKSCAN2 are less chategrized in previous studies. 

{% include figure.liquid loading="eager" path="assets/img/degs2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

After some biological interpretation of network built by SCENIC+. Now I want to discuss the performance of two pipelines. Here is a venn plot showing the number of TFs found by cistarget and homer. Homer found 310 TFs, cistaregt found 213 TFs. They have 109 overlaps. There are multiple reasons that could explain this discrepancy. For example, cistarget used a combined motif database. Homer uses its own database. Overall, the overlap provides some credentials to both tools. 

{% include figure.liquid loading="eager" path="assets/img/overlap.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Another aspect is to see how many target genes are found by both tools. I chose target genes of FOXA1 as example. This upset plot shows the similarities between scenic+ predicted gene list, t-gene predicted gene list, and gene list provided by encode chip-seq study. he leftmost bar represents the number of genes validated by the ENCODE project, which is a substantially larger set with over 10 thousand genes.The two smaller bars represent the number of genes predicted as targets by T-Gene (184 genes) and SCENIC+ (213 genes), respectively. The horizontal bars to the right indicate the intersection sizes between the sets. For instance, the rightmost bar corresponds to 94 genes that are common between all three sets: T-Gene predictions, Scenicplus predictions, and ENCODE validations.The bar with 131 indicates the number of genes shared between Scenicplus predictions and ENCODE validated genes, but not predicted by T-Gene.The bar with 90 shows the number of genes shared between T-Gene predictions and ENCODE validated genes, but not predicted by Scenicplus.The bar with 82 indicates the number of genes shared between T-Gene and Scenicplus predictions but not validated by ENCODE.
The main-takeaway is that about half of genes are commonly validated by ENCODE and predicted by both T-Gene and Scenicplus, which could suggest some confidence in these predictions.There are some genes that both T-Gene and Scenicplus agree upon but are not validated by ENCODE, which might be of interest for further investigation.

{% include figure.liquid loading="eager" path="assets/img/overlap2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Limitation and future directions**:

The main limitation is that bulk data cannot provide cell-level heterogeneity, which is very important as TFs can act differently in different cell types. 

Future direction would be find consensus regulons and conduct mutation analysis to see whether there are non-coding mutations in TF binding sites.

{% include figure.liquid loading="eager" path="assets/img/future_direction.jpg" title="example image" class="img-fluid rounded z-depth-1" %}



