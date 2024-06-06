---
layout: page
title: Advancing thoracic disease classification with CNNs
description: BMI707 group project
img: assets/img/cnn.jpg
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
{% include figure.liquid path="assets/img/cnn_pipeline.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Data preprocessing:

The datasets have been downloaded and formatted into appropriately pre-processed floating point tensors utilizing `Keras` packages. In particular, we employed the `class ImageDataGenerator` to retrieve and decode from PNG images (1024x1024) to RGB pixel grids, further transforming them into floating point tensors. Then, the pixel values were normalized from their original range of 0 to 255 to a new scale of [0, 1], which is optimal for neural network processing.

The ImageDataGenerator yields batches of 224x224 RGB images and 15 categorical labels. The batch size is 80. The resulting data batch shape is (80,224,224,3). Labels batch shape is (80,15)

Base Convolutional Network:

We first build a base convolutional network with four convolutional layers with each layer followed by a max pooling layer, and two dense layers. The activation function for intermediate layers is ReLU. We used softmax as an activation function for the top layer because softmax could give our model partial credits, thus increasing training effectiveness. The model is summarized in figure below(left). 

Since we have 15 classes, we use categorical cross entropy as a loss function where the target class is a label in a one-hot encoded vector. RMSprop with a learning rate of 2e-5 is used to optimize the loss. We used an accuracy metric to assess model performance. The model is trained for 30 epochs. For each epoch, 63 steps are used to yield from the training generator before declaring one epoch finished and starting the next epoch. 22 steps are used to yield from the validation generator.

Pre-trained Convolutional Network:

Then, we applied transfer learning on a pre-trained VGG19 model from ImageNet, with which we 
froze all layers except for the top classification layers to adapt the network to the specific features of our dataset. The conv_base VGG19 model is summarized in Appendix Figure 2. 

After the VGG19 base, the model structure includes one flattened layer, a fully connected layer with 512 units and ReLU activation, a dropout layer with a rate of 0.5 to help prevent overfitting, as well as the final output layer with 15 units corresponding to the number of disease classes, using softmax activation to output probabilities of each class. Similar to our custom CNN, we utilized the same optimizer and loss function for this multi-class classification task. The model is trained for 10 epochs.

In our second transfer learning model, we attempted a similar structure as above on a pre-trained Xception model on ImageNet. After the Xception base, the model structure includes one flattened layer, two fully connected layers with 512 units, and ReLU activation, a dropout layer with a rate of 0.5. In the first dense layer after the base, an elastic net kernel regularizer has been added to prevent further overfitting (l1 = 1e-4, l2 = 1e-4). The output layer, optimizer, loss function, and number of training epochs are the same as our first transfer learning model.


{% include figure.liquid path="assets/img/cnn_models.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

The evaluation metrics of the deep learning models:

AUC, Accuracy, F1 score, precision, and recall are used to evaluate the performance of our models. Accuracy is the simplest and most intuitive performance measure. It is defined as the ratio of correctly predicted observations to the total observations. Accuracy itself is not sufficient since it may be misleading for an imbalanced dataset. Precision measures the accuracy of positive predictions. Formulated as the ratio of true positives to the sum of true and false positives. Recall measures the ability of a model to find all the true positives within a dataset. Defined as the ratio of true positives to the sum of true positives and false negatives. Together, these metrics help in evaluating the model's performance when dealing with imbalanced datasets, offering insights into the trade-offs between identifying all relevant instances and maintaining a high level of prediction precision. The F1 Score is the harmonic mean of precision and recall, taking both false positives and false negatives into account. AUC measures separability. It tells how much the model is capable of distinguishing between classes. The higher the AUC, the better the model is at predicting 0s as 0s and 1s as 1s. By measuring the area under the curve, the AUC metric captures the overall performance of the model across all classification thresholds, thus providing a robust indication of its effectiveness. 

**Results**:

Deep learning models:

On the training set, the Base CNN has an accuracy of 0.6855 and a comparatively low loss of 1.2976. It achieved a high AUC of 0.9144. The precision and recall were closely matched at 0.6864 and 0.6786, respectively, leading to a balanced F1 score of 0.6825. VGG19 has a slightly lower accuracy of 0.6749 and a higher loss of 4.5783. Although its AUC dropped to 0.8546, it maintained a reasonable precision of 0.6845 and a recall of 0.6721, resulting in an F1 score of 0.6782. The Xception model, on the other hand, matched the Base CNN in accuracy at 0.6855 but suffered from a higher loss of 4.1506. Its AUC of 0.9013 was commendable, yet a notable disparity between a lower precision of 0.47 and a higher recall of 0.69 led to a reduced F1 score of 0.56, highlighting issues with generating too many false positives.

On the testing set, the Base CNN maintained consistent performance with a slight decrease in accuracy to 0.6810, a loss of 1.3701, and an AUC similar to training at 0.9142. The precision and recall remained well-balanced at 0.6805 and 0.6757, respectively, with an F1 score of 0.6781, showing that the model performs similarly well under testing conditions. VGG19 also achieved an accuracy of 0.6810 but with a higher loss of 1.7703. Its AUC decreased to 0.8933, yet it managed to slightly improve both precision and recall to 0.6810, leading to a higher F1 score of 0.6809 compared to training, indicating slight enhancements or adaptability in test scenarios. Xception maintained an accuracy of 0.6810 but with an increased loss of 4.3375 and an AUC of 0.8850. The precision dropped further to 0.46 while recall was at 0.68, resulting in a lower F1 score of 0.55.


{% include figure.liquid path="assets/img/cnn_deep.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

{% include figure.liquid path="assets/img/cnn_performance.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Non-deep learning models:

On the test set, the KMeans clustering approach has a low accuracy of 0.125, complemented by a precision of 0.57 and a recall of 0.12, resulting in an F1 score of 0.17. These metrics indicate that while KMeans was good at identifying specific classes when it predicted correctly, its overall ability to accurately label the dataset was significantly lacking, as evidenced by the low recall and, consequently, the poor F1 score. In contrast, the KNN (n=3) algorithm performed better, achieving an accuracy of 0.6445. Despite a lower precision of 0.49, the recall was considerably higher at 0.64, which led to a more balanced F1 score of 0.55. On the testing set, the performance trends observed were similar. KMeans clustering again showed limited effectiveness, with an accuracy of only 0.1047, precision at 0.54, and recall at 0.10, culminating in an F1 score of 0.14. On the other hand, KNN (n=3) has better performance in the training set, with an accuracy of 0.7144. It maintained a precision of 0.52 and an improved recall of 0.71, leading to a higher F1 score of 0.64. 

{% include figure.liquid path="assets/img/cnn_nondeep.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Limitations of the study:

In this study, we have eliminated multi-labeled samples that could result in the loss of valuable information useful for more detailed classifications or for understanding the underlying complexities of the data. This could limit the depth of insights our model can provide.

While the Base CNN shows promising results, the generalization capability of these pre-trained models to new, unseen data remains uncertain. The models' performance metrics are only as good as the data they are trained on, and if the training data lacks diversity or is not representative of the broader population, the model might not perform as well in real-world settings.

Challenges:

Even with strategies like stratified sampling and weight class, managing class imbalance can still pose significant challenges. While these methods help in creating a more balanced dataset for training models, they may only guarantee that the model will generalize well to real-world data, where the class distribution may still be highly skewed.
