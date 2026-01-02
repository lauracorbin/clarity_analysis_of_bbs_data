### script to analyze metabolite change from baseline to 36m in BBS ###
### initially created by Maddy Smith ###
### run from script dir 

####################
###### SET UP ######
####################

## read in parameter file
source("parameter_files/00_parameters.R")

## read in function file
source("00_functions.R")

## capture output
sink(file=paste0(log_output_dir, "01_association_analysis.txt"),split=TRUE, append=F)

## record todays date
today = Sys.Date()
today = gsub("-","_",today)

## load libraries
library(tidyr)
library(dplyr)
library(lme4)
library(nlme)
library(predictmeans)
library(partR2)
library(rsq)

###############################################################################
############################### READ IN DATA ##################################
###############################################################################

# read in filtered normalised data 
nmr_data <- read.table(filtered_data_file,h=T, sep=",")
  
# read in manifest that links metabolite data to study ID and timepoint
manifest <- read.table(nmr_manifest_file,h=T, sep=",")

# read in latest clinical data (from primary dataaset)
clinic_long <- read.csv(paste0(clinical_data_input_dir, "01_clinical_data_bbsmain_long.csv"), header = T)

## restrict clinical data to pilot patients
clinic_long <- clinic_long[clinic_long$studyid %in% manifest$Study.Identifier,]
dim(clinic_long)

###############################################################################
############################ PROCESS AND MERGE ################################
###############################################################################

# merge 3 data sources together - required data only
# restrict clinical to vars needed in model
essential_clinical <- clinic_long[,c("sample","studyid","timeinfreezer","age_at_rand","sex","Timepoint", "opdate","op_to_end", "baseline_to_op")]

# add sample id column to data to merge on
nmr_data$sample_id <- rownames(nmr_data)

# make id to merge on
manifest$timepoint <- gsub(" ","",manifest$Visit.Identifier)
manifest$sample <- paste0(manifest$Study.Identifier, "_", manifest$timepoint)

# bring sample id into clinical
essential_clinical$sample_id <- manifest$sample_id[match(essential_clinical$sample,manifest$sample)]

# merge metabolite and clinical
nmr_data_with_covar <- merge(essential_clinical,nmr_data,by=c("sample_id"))
dim(nmr_data_with_covar)

###############################################################################
############################### PROCESS DATA ##################################
###############################################################################

# make Timepoint a factor and set baseline as ref
nmr_data_with_covar$Timepoint <- as.factor(nmr_data_with_covar$Timepoint)
nmr_data_with_covar$Timepoint <- relevel(nmr_data_with_covar$Timepoint, "Baseline")

# count individuals
length(unique(nmr_data_with_covar$studyid))

###############################################################################
############################### SAMPLE SUMMARY ################################
###############################################################################

# restrict to unique individuals
nmr_data_with_cover_unique <- nmr_data_with_covar[!duplicated(nmr_data_with_covar[,c('studyid')]),]
dim(nmr_data_with_cover_unique)

# summarise time to surgery
summary(nmr_data_with_cover_unique$baseline_to_op)

# summarise surgery to end
summary(nmr_data_with_cover_unique$op_to_end)

###############################################################################
################################ RUN MODELS ###################################
###############################################################################

######################
# LINEAR MIXED MODEL #
######################

# run linear mixed model - with covars (TIF, age, sex)
set.seed(17082025)

# find start of nmr data
w <- which(colnames(nmr_data_with_covar)=="nonhdlc")
columns_with_metabolite_data <- c(w:ncol(nmr_data_with_covar))
length(columns_with_metabolite_data)

mixed_results <- metab.mixed.model(dtst=nmr_data_with_covar,columns_with_metabolite_data=columns_with_metabolite_data,
                                      covariables=c("timeinfreezer", "age_at_rand", "sex"), dataset_name=lmm_nmr)

# add metabolite name into data frame
mixed_results$metabolite_name <- rownames(mixed_results)
mixed_results <- mixed_results[,c(ncol(mixed_results),1:(ncol(mixed_results)-1))]

# add adjusted p-value
mixed_results$p_adjusted <- p.adjust(mixed_results$P_mixed, method = c("BH"))

# export nmr results
write.table(mixed_results, file=paste0(results_output_dir, "01_nmr_results_", today, ".csv"), sep=",", row.names=F)

###############################################################################
######################### WRITE OUT SAMPLE LIST ###############################
###############################################################################

# write out sample list after removing those with no op date or op after sampling for use in data prep for clarity
write.table(nmr_data_with_covar[,c(1:7)], file=paste0(clarity_input_dir, "sample_inclusion_info_", today, ".csv"), sep=",", row.names=F)

# write out metab list for keeping
write.table(names(nmr_data_with_covar)[w:ncol(nmr_data_with_covar)], file=paste0(clarity_input_dir, "feature_inclusion_info_", today, ".txt"),quote=F,row.names=F,col.names=F)

####################
######## END #######
####################

# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# clear workspace
rm(list = ls())
