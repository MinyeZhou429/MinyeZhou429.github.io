---
layout: page
title: Leukemia Diagnostics
description: Undergraduate senior design project
img: assets/img/leukemia.jpg
importance: 2
category: work
giscus_comments: true
---

**Title**:

Assessment of Dimension Reduction and Machine Learning Methods for Clinical Flow Cytometry Data Analysis


**Summary**: 

Test the performance of t-SNE, opt-SNE, UMAP, and other methods like kmeans and random forest on distinguishing cancer cells.

**Data Type**:

Flow cytometry data

**Introduction and rationale**:

Leukemia is a fatal disease that affects a big population in our society as the data shown here. And currently there is still a lot of challenges in diagnosis, prognosis, and optimize treatment. There are four types of leukemia: ALL, CLL, AML, and CML. In this project, we will only focus on CLL and ALL leukemia. We chose to use CLL and ALL because we have a better understanding of their cancer heterogeneity. 

{% include figure.liquid path="assets/img/leukemia_dis.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

There are many assays used for leukemia diagnosis, for example CBC, Nuclear Morphology, Cytogenetics, and Immunohistochemistry(IHC). The most widely used assay is flow cytometry. Flow cytometry has simple protoco.Base on flow cytometry data, doctors can plot significant cancer markers and plot them on 2d plots and select cancer cells from 2d plots. This is the process called manual gating. Leukemia diagnosis and prognosis required accurate and complete identification of cancer cells. 
Manual gating is subjective as you can see from the graph, the border of cancer cells versus non-cancer can be abrupt which results in inaccurate proportion of cancer cells.

{% include figure.liquid path="assets/img/rationale.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Study design**:

Feature selection is one of the most important step of ML. With this in mind, we began with some popular dimensional reduction tools such as t-sne, opt-sne, and UMAP. These nonlinear embedding methods can reduce a multi dimensional dataset into 2 dimension. Additionally, we also took another route in clustering method such as K-means and HDBSCAN. These tools will cluster n samples into a defined number of clusters. Different machine learning models were used to provide sample-level diagnosis. Finally, we constructed a pipelines that could successfully provide tumor burden estimation.

{% include figure.liquid path="assets/img/pipeline.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Data Description**:

The project used both CLL and ALL dataset. We got our CLL dataset from UC San Diego and ALL from stanford university. It is worth noting that only a portion of the CLL has tumor burden label which we used as ground truth to validate our tumor burden estimation. At the same time, CLL dataset contains samples from patients’ peripheral blood+bone marrow and ALL contains samples purely from bone marrow. 2 panels are design for CLL with 16 different markers measured and 11 markers are measured for ALL.

{% include figure.liquid path="assets/img/leukemia_data.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Results**:

opt-SNE and t-SNE both fail on identifying the whole set of cancer cells. In the figure below, significant cancer markers for cells in the potential cancer cluster are plotted on 2D plots. The right-up region is the expert gating region. However, there are many cells out of the expert gating region. This indicates that the potential cluster contains both cancer and non-cancer cells. We tested on different potential clusters and similar results were obtained.

{% include figure.liquid path="assets/img/tsne1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

With the increase in file size, runtime of t-SNE follows a rapidly increasing trend while opt-SNE has a more gentle slope of increasing. This means that opt-SNE is more efficient than t-SNE when encountering large datasets. It is worth mentioning that t-SNE fails to output results when file size reaches a length of 1000000 lines. Moreover, the cross-map of potential cancer clusters found by opt-SNE on the original t-SNE plot shows that while t-SNE has a longer runtime and more iterations, it does not necessarily perform better than opt-SNE on clustering.

{% include figure.liquid path="assets/img/runtime.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
{% include figure.liquid path="assets/img/tsne2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

While we found t-SNE and opt-SNE do not work very well, we tried a different dimension reduction method called UMAP. For CLL cohort, we used two million cells as the training dataset and one million cells as the testing dataset, half of it is from cancer patients and half of it is from healthy patients. As shown in the figure on the left, the red data points are from cancer patients, while the blue data points are from healthy patients. We would compare this two plots and found the non-overlapped regions as the potential cancer clusters. In this example, we will use the cluster in the box to show how do we determine if this cluster is a cancer cluster. We plotted all the cells from that cluster onto 2D plots with different markers. As we can see in figure C, the positions of the cells are consistent with the typical CLL cancer cell phenotype. Therefore, we can determine this cluster is a cancer cluster. In addition, we found a very small p-value that shows there is a significant difference in cancer cells between healthy samples and cancer samples from the cancer clusters we determined. Lastly, we tried to plot the cancer cells we found in UMAP onto opt-SNE, shown in figure D. We found the cancer cells appear in the none cancer regions of opt-SNE. This further demonstrates that opt-SNE failed to find cancer cells.

{% include figure.liquid path="assets/img/umap1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Next, we used the same UMAP pipeline on ALL data. The red data points are from cancer patients, while the blue data points are from healthy patients. We went through the same process as CLL. However, the positions of the cells are not consistent with the typical ALL cancer cell phenotype. Therefore, we selected a few other clusters of cells but are not sure they are the complete and accurate set of cancer cells

{% include figure.liquid path="assets/img/umap2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

We just showed that how the process can be done through visualization. However, this whole process can be done automatically. We first tried HDBSCAN. HDBSCAN is a clustering method that is build into UMAP package. The first step is to create a UMAP; then apply HDBSCAN; next, we would compare cluster proportion and identify cell-based biomarkers; lastly we would use biomarkers to classify samples and estimate tumor burden. We can also use linear svc with UMAP transformation to reach the same goal. We tested this method manually with a small dataset of 7 samples and a 1% cut off. And, we reached a sample level prediction with 100% accuracy. In the future, we hope to test on a larger dataset

{% include figure.liquid path="assets/img/hdbscan.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Then, we are also curious whether a simple k-means clustering on the original data space can be used to do diagnosis and tumor burden estimate without using UMAP or t-SNE Here, we perform clustering on both CLL and ALL samples. Specifically, we select several files from each leukemia type as the template and perform k-means to identify the cluster centroids. 2 different templates are built for ALL because of its high cancer cell heterogeneity. Because of the high run-time of k-means on such large dataset, we separated the template from training to identify cancer cluster. Here, we selected around 20 samples each for CLL and ALL to identify cancer clusters. We use Wilcoxon Rank Sum test to select clusters that are significantly amplified in cancer clusters and deprived in healthy samples. These would be our cancer clusters. Finally, to evaluate the cancer clusters we identified from k-means, we used around 100 samples for CLL and ALL to calculate the tumor burden. In this case, we will compare the tumor burden output from our pipeline to the doctor’s label and evaluate the pearson’s r correlation. 

{% include figure.liquid path="assets/img/burden1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Our initial tumor burden estimation reached a correlation of around 0.94. However, the false positive rates is high. We adjusted our p-values using the FDR correction and significantly reduced the false positive rates. 

{% include figure.liquid path="assets/img/burden2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Overall ALL has higher heterogeneity of cancer so it is more challenging to estimate ALL tumor burden than in CLL Two templates of ALL similar correlation of around 0.85, which indicates that our pipeline has robust performance

{% include figure.liquid path="assets/img/burden3.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

We also took the output of our pipeline to train a Random Forest model. We used a training and testing ratio of 7:3 and all of the datasets have achieved very high accuracy and f1 scores.
{% include figure.liquid path="assets/img/other_methods.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Conclusion**:

1. UMAP is superior to t-SNE and opt-SNE in terms of identification of CLL cells

2. Our pipeline using k-means clustering could output accurate tumor burden estimation for CLL

3. Constructing a random forest model using the clusters from our pipeline helped sample level diagnostics for both CLL and ALL

