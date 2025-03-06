# PDAC-Survival-Prediction-Using-Deep-Learning-and-Tile-Analysis

## Context and Objective

Pancreatic Ductal Adenocarcinoma (PDAC) is an aggressive and highly fatal form of pancreatic cancer. It is the second most common digestive cancer globally and the fourth leading cause of cancer-related death. The treatment decisions for PDAC are currently based on clinical criteria such as RECIST (Response Evaluation Criteria in Solid Tumors).

This project aims to use deep learning techniques to predict patient survival in PDAC by analyzing digitized tissue slides, along with RNA-seq and Exome data. We combine Convolutional Neural Networks (CNNs) to classify tumor and normal tissue regions from the samples and apply a Multiple Instance Learning (MIL) approach to improve survival prediction.

## Data

The data used in this study comes from two cohorts:
- **Private Cohort from Besançon** (n=206): Non-metastatic patients treated with surgery.
- **TCGA Cohort** (n=166): Non-metastatic patients treated with surgery.

The types of data included are:
- Whole tissue slides stained with Hematoxylin and Eosin (H&E) or Hematoxylin, Eosin, and Safran (HES).
- RNA-seq and Exome data.

## Methodology

### 1. Tile Creation
We used **Groovy QuPath** to divide tissue slides into smaller tiles, allowing us to analyze specific tissue regions, such as tumor and normal areas.

### 2. Image Normalization
The images were normalized using the **Reinhard** method to standardize color distributions across slides, making it easier for the model to learn from them.

![Image of aciduino on protoboard](https://github.com/dinaOuahbi/DenseNet-_-prediction-of-responders-and-non-responders-patients/blob/main/filtration_norm_slide.PNG)

### 3. CNN Model Construction
Two separate CNN models were created:
- **CNN1**: Stroma vs. Other tissue classification.
- **CNN2**: Normal tissue vs. Tumor tissue classification.

### 4. Predictions and Model Combination
The results from both CNN models were combined to predict patient survival outcomes.

### 5. Independent Tile Prediction
Tiles were predicted independently to assess their contribution to survival prediction.

### 6. Grad-CAM Analysis
We used **Grad-CAM** to visualize the importance of specific regions in the images, helping to identify which parts of the tiles influenced the survival predictions the most.

### 7. Survival Model with Intermediate Layers
We integrated intermediate layers from the CNN model into a survival model to assess the impact of learned features on patient survival.

### 8. MIL (Multiple Instance Learning) Approach
The MIL approach was used to treat the tiles as "instances" and predict survival based on a set of heterogeneous instances, which enhanced prediction accuracy on complex data.

## Conclusion
This study demonstrates that deep learning can be effectively applied to predict the prognosis of PDAC and assist in making better decisions for adjuvant treatments. This model can provide valuable support in selecting the most appropriate treatments for patients with PDAC.
https://www.mdpi.com/2227-9059/12/12/2754/pdf

