READMEWritten 8/8/22 by Ryan Otto
#####################################################################################################
The following code was used to perform all analysis and generate all figures in:

Otto R. M., Turska-Nowak A., Brown P. M., Reynolds K. A. "A continuous epistasis model for predicting growth rate given combinatorial variation in gene expression and environment.

Copyright (C) 2022 Kimberly A. Reynolds
This collection of python notebooks is free software distributed under the BSD 3-clause license, please see the file LICENSE for details
######################################################################################################
Code was written using Python 3.9.12 with Jupyter Notebook, and uses pickle (v4.0), numpy (v1.21.2), scipy (v1.7.1), matplotlib (v3.5.2), pandas (v1.4.2), scikit-learn (v1.0.1), seaborn (v0.11.2), biopython (v1.79), and regex (v2021.8.21).

It is presented as a series of seven consecutive Jupyter notebooks, many of which use data structures generated by previous notebooks. This code computes CRISPRi repression intensity for single sgRNA perturbations and growth rate effects for up to fourth-order CRISPRi and environmental perturbations. The code then trains and evaluate the continuous epistasis model as a methodology to predict changes to E. coli growth rate following continuous, high-order transcriptional and environmental perturbations.

If you use anaconda for package management, creating the following environment using the command line code below will ensure that each required package is installed and uses a version compatible with this our analysis. Replace <ENVIRONMENT_NAME> with your desired environment name.

conda create -n <ENVIRONMENT_NAME> python=3.9.12 numpy=1.21.2 scipy=1.7.1 matplotlib=3.5 pandas=1.4.2 scikit-learn=1.0.1 seaborn=0.11.2 biopython=1.78 openpyxl=3.0.1 regex=2022.7.9 jupyter nbconvert==5.4.1
######################################################################################################
Contents:

.gitignore				File stopping git from committing files that will be generated 
					locally by this repo.

1_qPCR_measurements.ipynb		Code that quantifies gene repression intensity from qPCR data.

2_Pairwise_HiSeq_QC.ipynb		Basic quality control run on sgRNA counts for the pairwise 
					library.

3_Pairwse_Growth_Rates.ipynb		Code quantifying relative growth rates following pairwise 
					CRISPRi.

4_Continuous_Epistasis_Model.ipynb	Code training the continuous epistasis model on pairwise data.

5_Third_Order_CRISPRi.ipynb		Code evaluating the continuous epistasis model on third-order 
					CRISPRi perturbations.

6_Env_Growth_Rates.ipynb		Code quantifying growth rates for gene-by-environment data.

7_Env_Model_Fitting.ipynb		Code training the continuous epistasis model on gene-by-
					environment data.

Figures/				Folder containing figure panels generated by this code.

HPC_Analysis/				Copies of code run on our HPC to count sgRNA abundances from 
					HiSeq data.

input_files/				All input data required for this analysis.

intermediate_files			Data files generated by this analysis.

LICENSE					BSD 3-clause license.

plot_defaults.py/			Script to standardize figure formats between notebooks.

README.txt				This file.

Supplementary_Tables.xlsx		Excel file containing all Supplementary Tables. Data in these 
					tables can be updated using the notebooks provided.