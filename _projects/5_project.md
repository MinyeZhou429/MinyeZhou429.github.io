---
layout: page
title: Advancing thoracic disease classification with CNNs
description: BMI707 group project
img: assets/img/bc.jpg
importance: 6
category: work
---

**Title**:

Advancing Thoracic Disease Identification with CNNs 

**Members and Contributors**:

Minye Zhou, Michelle Dai, Evelyn Jiang

**Data type**:

Medical Images

**Summary of Dataset**:

We used the ChestX-Ray14 dataset. This dataset contains 112,120 frontal-view X-ray images from 30,805 distinct patients, each potentially labeled with one or more of fourteen identified diseases. These labels were extracted from corresponding radiological reports using natural language processing techniques, which are expected to have an accuracy greater than 90%. The dataset identifies fourteen prevalent thoracic conditions using number 0 to 14: No Finding: 0, Atelectasis: 1, Cardiomegaly: 2, Effusion: 3, Infiltration: 4, Mass: 5, Nodule: 6, Pneumonia: 7, Pneumothorax: 8, Consolidation: 9, Edema: 10, Emphysema: 11, Fibrosis: 12, Pleural_Thickening: 13, Hernia: 14. 

We split the data into train and test sets, containing 8000 and 2000 examples respectively. We found that among 8000 training samples, 1187 samples have multiple labels. For 2000 testing samples, 304 samples have multiple labels (more than one thoracic condition). We excluded those multi-labeled samples. Additionally, we allocated 25% of the data from the training set to serve as a validation set to monitor and tune the model's performance during training.


**Introduction**:

Chest X-ray imaging stands as a cornerstone in radiology. With advanced deep learning technologies, significant strides have been made in the automated analysis of medical images. Particularly, advancements in disease classification through these methods have seen considerable improvements over traditional approaches. Our study harnesses the "ChestX-ray14" dataset, with its detailed labeling, to investigate disease classification and pathology localization within X-ray images.

The primary objective of our research is to evaluate the precision of various Convolutional Neural Network (CNN) architectures in identifying specific thoracic diseases from X-ray images. This involves a comprehensive comparison of the performance of different CNN models to assess their efficacy in classifying a broad spectrum of thoracic conditions. Furthermore, our study seeks to adapt and finetune pre-trained ImageNet models for this task.


**Methods**:



{% include figure.liquid path="assets/img/sc6.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
