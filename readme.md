# FDA Analysis: Data & Source Code

This repository contains the dataset and source code for the empirical study on **Unit 1**. It includes the complete pipeline from raw data preprocessing and sample extraction to clustering and Huber-loss-based regression analysis.

## Repository Structure

```text
├── raw_data/                      # Original and preprocessed datasets
├── extraction_sample/             # Samples extracted based on load trends
├── cluster_result/                # Clustering outputs utilized in the paper
├── regression_result/             # Regression analysis results
├── fulldata_with_dryness.xlsx     # Dataset incorporating calculated dryness metrics
├── steam.ipynb                    # Jupyter Notebook for physical dryness constraints
├── varcluster_coal.R              # Script for variable/data clustering pipeline
├── varselect_coal.R               # Script for executing regression models
├── temporal cross-validation      # Temporal cross-validation
└── varselect_fixload_huberloss.R  # Wrapper functions for Huber loss & feature selection
```
