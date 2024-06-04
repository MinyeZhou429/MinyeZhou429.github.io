---
layout: page
title: Fluids Resuscitation of Patients with Congestive Heart Failure and/or End Stage Renal Disease
description: MIT HST.953 Group project
img: assets/img/mit_intro.jpg
importance: 3
category: work
---

**Title**: 

Model assessments for Fluids Resuscitation of Patients with Congestive Heart Failure and/or End Stage Renal Disease

**Contributor and group members**:

Students: Yichen Wang, Minye Zhou, Amal Cheema, Sud Perera, Jin Kim

Advisors: Max Shen, Sicheng Hao, Xiaoli Liu, Leo Celi

**Data type**:

Clinical data

**Introduction**:

According to the CDC, sepsis affects more than 1.7 million patients in the US each year. Among them 270,000 die each year. It currently accounts for more than half of hospital deaths in the US, and it is the leading cost in the entire US healthcare system, with an estimated $60B spent on sepsis each year.

Septic patients are treated in the intensive care unit (ICU) as they require close monitoring. Alongside broad spectrum antibiotics, early and often aggressive resuscitation with IV fluids is one of the traditional cornerstones of sepsis management to stabilize sepsis-induced tissue and organ hypoperfusion. 

However, some patients have medical conditions such as congestive heart failure (CHF) and/or end stage renal disease (ESRD) and are at risk of fluid overload secondary to sepsis treatments that can exacerbate organ failure. Uncertainty of how to fluid resuscitate patients with ESRD or CHF leads to variable clinical practices, which demonstrates the need for evidence-based modeling that can inform clinical practice.

{% include figure.liquid loading="eager" path="assets/img/mit_intro.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Goal**:

1. Determine the relationship between fluid resuscitation volume and outcomes in patient populations experiencing sepsis with or without CHF or ESRD
2. Construct models to determine the optimal quantity of fluid resuscitation for these patients subgroups to ultimately inform medical management

**Strategy**:

1. Extract data from two large, publicly available  ICU databases: MIMIC-IV and eICU
2. Create a unified dataset with all patients that meet inclusion criteria
3. Iteratively refine our dataset, models, and results as we gain deeper understanding of our data and our “business” (CRISP-DM)

**Method**:

Data source and Study population:

MIMIC-IV and eICU

Inclusion criteria:

1. Admitted from the ED
2. Met the Sepsis3 criteria during the first day
   
Exclusion criteria:

1. CKD stages I through IV
2. ICU readmission within the same hospital stay

Data sanity check: 

1. Assessed missingness and cohort sizes
2. Illustrated need for treatment limit ranges (.e.g, clerical errors in lab value entry)

{% include figure.liquid loading="eager" path="assets/img/mit_design.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Covariate selection:

1. Based on our subject-matter knowledge
2. Scrutinized the missing percentage
3. Used backward stepwise models as reference

{% include figure.liquid loading="eager" path="assets/img/mit_feature.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Statistical Analysis:

1. Used ANOVA and t-test to compare the fluid volumes among three subgroups
2. Propensity score weighting as foundation for future modeling
3. Fitted a generalized linear model and Random Forest Regression to obtain the generalized propensity score (GPS)
4. Calculated and compared Standardized mean difference (SMD) 

Model evaluation:

1. Splitting dataset
2. Calculate R Squared and RMSE

**Results**:

It was found that the sepsis-only group exhibited the lowest mortality. However, certain covariates were unevenly distributed among the subgroups, influencing the outcomes. For instance, the group with Congestive Heart Failure (CHF) had the oldest median age. Additionally, the group with End-Stage Renal Disease (ESRD) had the highest median first-day Sequential Organ Failure Assessment (SOFA) score.

{% include figure.liquid loading="eager" path="assets/img/res1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

In comparing fluid administration across patients with Congestive Heart Failure (CHF), End-Stage Renal Disease (ESRD), and sepsis-only conditions, significant differences were observed (p<0.05). The highest fluid volume in the first 24 hours was administered to sepsis-only patients, averaging 36.1±27.9 cc/kg with a median of 20.1 cc/kg. In contrast, the lowest fluid volume was observed in patients with ESRD, averaging 17.6±20.6 cc/kg with a median of 10.6 cc/kg.

{% include figure.liquid loading="eager" path="assets/img/res2.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

The selected covariates explained only 17-26% of the variations in fluid volumes. The models exhibited high Root Mean Square Error (RMSE), particularly for sepsis patients diagnosed with CHF and/or ESRD. This unsatisfactory performance indicates a need for further work on covariate selection and model engineering to improve the accuracy and reliability of the predictions.

{% include figure.liquid loading="eager" path="assets/img/res3.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

{% include figure.liquid loading="eager" path="assets/img/res4.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

Random Forest outperformed Multiple Linear Regression. RF was more likely to balance covariates for those with high unweighted Standardized Mean Differences. Despite this, some covariates still exhibited a small-to-medium degree of imbalance (SMD > 0.2) across the four GPS strata. The poorest covariate balance was observed in sepsis patients with ESRD, which could be attributed to the limited subpopulation size of 1,421.

{% include figure.liquid loading="eager" path="assets/img/res5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}

**Conclusion and Discussion**:

Model comparison:  

Based on RMSE and R squared, both GLM and RF models did not have high predictive power. RF model showed a better capability of balancing the covariates across different administered fluid volumes. 

Clinical significance: 

A model that predicts how much fluid to administer to a septic patient with ESRD and/or CHF could have utility in informing clinical practice and reducing patient mortality (if deployed appropriately).







    
