# Minimum-Cut-Algorithm-example
This repository contains example code and datasets for the paper: "Optimization-Based Framework with Flux Balance Analysis (FBA) and Metabolic Pathway Analysis (MPA) for Identifying Metabolic Objective Functions"

The repository provides case study data, metabolic reaction content, flux bounds for E. coli, iCAC802, iJL680, hybrid Cac, and hybrid Clj metabolic models, as well as supplemental MATLAB and Python codes for running graph-based analysis simulations described in the article.

# 📂 Repository Contents
The repository is organized into five main script folders, each corresponding to case studies and supporting analyses described in the manuscript.

Script 1 – Case Study 2 (TIObjFind CoI)
Files:
1. Cac_Clj_MuliFBA_example.mat
2. multiOBJ_FBA.m
3. README.rtf

Script 2 – Min-Cut Graph & Coefficient of Importance (CoI)
Files:
1. FBA_V.mat
2. mc_plot_coculture_system.m
3. This script generates the flux-dependent graph using the minimum cut algorithm and calculates CoI values.

Script 3 – Case Study 2 (Comparison of Objectives)
Files:
1. singleOBJ_FBA.m
2. Cac_Clj_MuliFBA_example.mat
This script demonstrates the co-culture system case study, including comparison between ObjFind results/ single-reaction objectives.

Script 4 – E. coli Case Study (Supplementary Information Example)
Data reference: Orth, J. D.; Fleming, R. M. T.; Palsson, B. Ø. Reconstruction and Use of Microbial Metabolic Networks: The Core Escherichia Coli Metabolic Model as an Educational Guide. EcoSal Plus 2010, 4 (1), ecosalplus.10.2.1.
This script demonstrates flux-dependent graph construction, minimum cut application, and Sankey diagram visualization for the E. coli case study.

Script 5 – Cac Case Study with 13C Data (Supplementary Information Example)
Data reference: Au, J., Choi, J., Jones, S. W., Venkataramanan, K. P., & Antoniewicz, M. R. (2014). Parallel labeling experiments validate Clostridium acetobutylicum metabolic network model for 13C metabolic flux analysis. Metabolic Engineering, 26, 23–33.
This script illustrates the case study of 13C-validated metabolic flux aligns with the use of model.

# 📊 Data & Results
E. coli case study (Glucose condition):
Sankey diagrams generated using
node_names.csv, node_values_lim.csv, sankey_plot_Ecoli.py, and DATA SET S1.xlsx.

Co-culture system case study:
Flux balance simulations using
singleOBJ_FBA.m, multiOBJ_FBA.m, and Cac_Clj_MultiFBA_example.mat.

Multi-objective FBA weighting:
Use FBA_V.mat and mc_plot_coculture.m to generate flux-dependent graphs and edge densities.
Final normalized weights are provided in DATA SET S1.xlsx.

# 🛠️ Installation & Requirements
To run the MATLAB-based FBA and graph analysis:
Install the COBRA Toolbox
Install a supported solver such as GUROBI (academic license available):

# 🛠️ How to Run
Download the relevant script folder.
Load required .mat data files before running scripts.
For Sankey diagrams and plotting, install Python dependencies (e.g., matplotlib, plotly).
Run MATLAB scripts for FBA analysis and graph construction.
