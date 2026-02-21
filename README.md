# Systemic-Metabolic-Analysis

This repository contains the code used for the analysis described in the manuscript titled:
"Smoking reshapes systemic metabolic coordination across lung cancer subtypes revealed by total-body PET/CT and tumor habitat imaging".

**Overview**

This repository includes scripts for conducting both population-level metabolic network analysis and individual-level metabolic network analysis, as well as tumor habitat imaging analysis for clustering individual and population-level tumor habitats. These analyses were applied to data from lung cancer patients and healthy controls using total-body 18F-FDG PET/CT.

The main components of the code are:

· Population network analysis: This script processes static PET data from the population cohort to create a systemic metabolic network analysis that reveals subgroup-specific metabolic alterations and the impact of smoking.

· Individual network analysis: This script constructs metabolic networks at the individual level, leveraging dynamic PET data to capture personalized metabolic alterations across lung cancer subtypes.

· Habitat analysis-SLIC: A script for spatially clustering tumor habitats using Simple Linear Iterative Clustering (SLIC) algorithm. This step is essential for mapping the tumor microenvironment. 

· Habitat analysis-Leiden bootstrap consensus clustering: This script is used for clustering tumor habitats based on metabolic features and tumor microenvironment characteristics. It utilizes the Leiden algorithm and bootstrap resampling for robust cluster identification.



**How to Run**

1. Clone or download this repository to your local machine.

2. Ensure you have the required software installed (MATLAB, Python).

3. Follow the instructions for each file:

      · Population network analysis: Use this script to generate population-level network analysis from static PET images.

      · Individual network analysis: This script should be run for individual patient PET data to calculate organ-specific metabolic abnormalities.

      · Habitat analysis-SLIC: This script can be used to segment and cluster tumor habitats based on PET and CT images.

      · Habitat analysis-Leiden bootstrap consensus clustering: Run this script after preparing the necessary PET/CT tumor data.



