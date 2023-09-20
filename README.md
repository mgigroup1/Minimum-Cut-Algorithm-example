# Minimum-Cut-Algorithm-example
example code for the paper "An Integrated Approach to Graph-based Analysis of Metabolic Networks Using Minimum Cut Algorithm and Flux Balance Analysis"

This includes all case study data; metabolic reaction content; flux bounds of E coli, iCAC802, iJL680, hybrid Cac, and hybrid Clj metabolic models; 
and supplemental MATLAB codes for running the graph-based analysis simulations described in this article. 

# Details
1. The Sankey diagrams for the E. coli case study, Gluc condition: download the files of node_names.csv; node_values_lim.csv; sankey_plot_Ecoli.py; DATA SET S1.xlsx
2. The predicted profiles for the co-culture system case study: download the files of singleOBJ_FBA.m, multiOBJ_FBA.m, Cac_Clj_MultiFBA_example.mat
3. The Ecoli_case_study.m file demonstrates the case study of Ecoli, including the construction of the flux-dependent graph & application of the min-cut algorithm & criteria of edge density. Download ECOLI_example.mat file and load data before running.
