## Summary

This repo contains code to support analyses presented in the manuscript: "Multivariate analysis of metabolomic data to identify biological pathways modified by a clinical intervention"

Code contributed by Laura Corbin (LJC) and Rachel Wood (RMW).

## Files

### experimental
Code relating to the analysis of metabolomics data from the By-Band-Sleeve trial.

#### clarity_analysis (RMW)
01\_cross\_validation.R - main CLARITY analysis. Also runs cross-validation procedure.  
02\_make\_plots.R - code to make Figure 4 - CLARITY plots for experimental data.

#### data_prep (LJC)
01\_prepare\_data.R - Code to prepare data for LMM and CLARITY analysis.

#### linear_analysis (LJC)
00\_functions.R - Functions called by 01_association_analysis.R.  
01\_association\_analysis.R - Code to run an LMM analysis of data.  

#### results_process
00\_functions.R - Functions called by 01 script. Edited versions of functions from bubblehatmap package. 
01\_explore\_hits.R - exploration of clarity results including heatmap generation. 
01a\_feature\_labels.R - relabelling of features to work with bubbleheatmap package.  
02\_compare\_results.R - merging of clariity and lmm results to compare.  
03\_tables\_and\_figures.r - code to make figures and tables for manuscript.  


### simulation (RMW)
Code relating to the generation and analysis of simulated data.

00\_lmm\_funs.R - LMM function. 
00\_sim\_funs.R - function to simulate data.  
01\_sim\_example.R - analysis of simluated data by LMM and CLARITY.  
02\_sim\_sensitivity\_analysis.R - sensitivity analyses.  
03\_plot\_sensitivity\_analysis.R - plots of sensitivity analyses. 