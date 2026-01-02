# code to take By-Band-Sleeve NMR data post-filtering and generate covariance matrices for baseline and endpoint.
# covar matrices will be used by Dan Lawson to look for multivariate signatures of change associated with the weight loss intervention.

##############################################################################################################################
### SET UP ###################################################################################################################
##############################################################################################################################

## read in parameters 
source("parameter_files/00_parameters.R")

## load libraries 
library(dplyr)
library(tidyverse)


## record todays date
today = Sys.Date()
today = gsub("-","_",today)
today

## capture output
sink(file=paste0(log_output_dir, "01_prepare_data_", today, ".log"),split=TRUE)

##############################################################################################################################
### READ IN DATA #############################################################################################################
##############################################################################################################################

# read in QC'd metabolite data
nmr_data <- read.table(data_file,h=T)

# read in manifest that links metabolite data to study ID and timepoint
manifest <- read.table(manifest_file,h=T, sep=",")

# read in data restricted to final filtered dataset (after removal of those without an op date or will op after 36mth sampling)
#filtered <- read.table(sample_input_file,h=T, sep=",")

# read in feature annotations to identify derived vars
nmr_feature_metadata <- as.data.frame(read_delim(metadata_file,delim="\t"))

# read in latest clinical data (from primary dataaset)
clinic_long <- read.csv(paste0(clinical_data_input_dir, "01_clinical_data_bbsmain_long.csv"), header = T)

## restrict clinical data to pilot patients
clinic_long <- clinic_long[clinic_long$studyid %in% manifest$Study.Identifier,]
dim(clinic_long)

###############################################################################
############################ PROCESS AND MERGE ################################
###############################################################################

# merge 3 data sources together - required data only
# restrict clinical to vars needed to filter
essential_clinical <- clinic_long[,c("sample","studyid", "opdate","op_to_end", "baseline_to_op")]

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

##############################################################################################################################
### EXCLUDE INELIGIBLE #######################################################################################################
##############################################################################################################################

# remove participants who don't have an op date (haven't had surgery), and who have opdate after visit 6 (negative op_to_end)
nmr_data_with_covar <- nmr_data_with_covar[!is.na(nmr_data_with_covar$opdate),]
dim(nmr_data_with_covar)

nmr_data_with_covar <- nmr_data_with_covar[(nmr_data_with_covar$op_to_end > 0),]
dim(nmr_data_with_covar)

# drop clinical vars and return to original nmr_data structure
row.names(nmr_data_with_covar) <- nmr_data_with_covar$sample_id
nmr_data <- nmr_data_with_covar[,c(7:ncol(nmr_data_with_covar))]
dim(nmr_data)

##############################################################################################################################
### RESTRICT FEATURES TO NON-DERIVED #########################################################################################
##############################################################################################################################

# extract list of non-derived features
non_derived <- nmr_feature_metadata[which(nmr_feature_metadata$derived_features == "no"),c("feature_names")]

# drop apobapoa1 (wasn't identified as ratio but is one) 
non_derived <- non_derived[!non_derived %in% c("apobapoa1")]

# restrict dataset to those non-derived
nmr_data <- nmr_data[,which(names(nmr_data) %in% non_derived)]
dim(nmr_data)

##############################################################################################################################
### SPLIT BY TIMEPOINT #######################################################################################################
##############################################################################################################################

baseline_ids <- manifest$sample_id[manifest$Visit.Identifier == "Baseline"]
endpoint_ids <- manifest$sample_id[manifest$Visit.Identifier == "36 months"]

# standardise
nmr_data_z <- nmr_data %>% mutate_all(~(scale(.) %>% as.vector))

# split by timepoint
baseline_data_z <- nmr_data_z[which(row.names(nmr_data) %in% baseline_ids),]
endpoint_data_z <- nmr_data_z[which(row.names(nmr_data) %in% endpoint_ids),]

##############################################################################################################################
### GENERATE COVARIANCE MATRICES #############################################################################################
##############################################################################################################################

# covar after standardising/z-scoring
baseline_covar_z <- cov(baseline_data_z, use="complete.obs")
endpoint_covar_z <- cov(endpoint_data_z, use="complete.obs")

##############################################################################################################################
### SAVE OUT #################################################################################################################
##############################################################################################################################

write.table(baseline_covar_z,paste0(output_dir,"baseline_z_covar_matrix.csv"),sep=",",col.names=T,row.names=T)
write.table(endpoint_covar_z,paste0(output_dir,"endpoint_z_covar_matrix.csv"),sep=",",col.names=T,row.names=T)

write.table(nmr_data_z,paste0(output_dir,"standardized_after_filtering.csv"),sep=",",col.names=T,row.names=T)
            
##############################################################################################################################
### FINISH ###################################################################################################################
##############################################################################################################################

# capture session info
print("Session information:")
sessionInfo()

# save output
sink()

# clear working space
rm(list = ls())

