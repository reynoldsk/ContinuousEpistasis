READMEWritten 8/8/22 by Ryan Otto
Updated 4/20/23 by Ryan Otto
Updated 10/12/23 by Ryan Otto
#########################################################################################################################
The following files are all data input used in this publication.

Otto R. M., Turska-Nowak A., Brown P. M., Reynolds K. A. "A continuous epistasis model for predicting growth rate given combinatorial variation in gene expression and environment.

Copyright (C) 2023 Kimberly A. Reynolds
This collection of python notebooks is free software distributed under the BSD 3-clause license, please see the file LICENSE for details
#########################################################################################################################
Contents:
20200923_glu_growth_rates_TableS5.csv		Table taken directly from Mathis et al., 2021 containing growth rate 							measurements.
bootstrapped_files/				Files generated on one bootstrapping iteration for both the pairwise 							CRISPRi library and the gene-by-environment data to estimate variability 						in model predictions and accuracy between separate single-perturbation 							experiments.
date_to_gene.txt				File relating qPCR experiment dates to the genes measured on those days.
date_to_strainID.txt				File relating qPCR experiment dates to the strains measured on those 							days.
env_norm_ind_GR_A.pickle			File containing normalized replicate growth rates of E. coli in various 						media conditions (part 1).
env_norm_ind_GR_B.pickle			File containing normalized replicate growth rates of E. coli in various 						media conditions (part 2).
env_raw_GR_A.pickle				File containing raw mean growth rates of E. coli in various media 							conditions (part 1).
env_raw_GR_B.pickle				File containing raw mean growth rates of E. coli in various media 							conditions (part 2).
no_filter_counts/				Pairwise CRISPRi sgRNA counts before filtering.
pairwise_SID_ecoli.pickle			File containing pairwise sequence identity distributions between all 							Enterobacterales species with orthologs to the 19 genes tested in the 							pathway CRISPRi experiment.
pathway_counts/					Folder containing all NGS counts for the pathway CRISPRi experiment.
pathway_lib.fasta				File containing sgRNA sequences for all sgRNAs in the pathway CRISPRi 							experiment.
optimization_files/				Files generated beforehand used to optimize subsampling percentage and 							regularization of the continuous epistasis model.
q30_filter/					Pairwise CRISPRi sgRNA counts after filtering.
qPCR/						Folder containing all qPCR results.
README.txt					This file.
strainID_to_sgrna.txt				File relating strain ID to the sgRNA that strain carried.
third_order_counts.pickle			File containing sgRNA counts from the third-order library.
turbidostat_GR.pickle				File containing growth rate information for the pairwise and third-order 						libraries.