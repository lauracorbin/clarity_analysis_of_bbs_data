# code to compare clarity results to LMM results

##############################################################################################################################
### SET UP ###################################################################################################################
##############################################################################################################################

## read in parameters 
source("parameter_files/00_parameters.R")

## load libraries 
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(bubbleHeatmap)
library(Hmisc)
library(tidyverse)
library(Cairo)


## generate new feature labels
source("01a_feature_labels.R")

## capture output
sink(file=paste0(log_output_dir, "02_compare_results.log"),split=TRUE)

## record todays date
today = Sys.Date()
today = gsub("-","_",today)
today

##############################################################################################################################
### READ IN DATA #############################################################################################################
##############################################################################################################################

# read in filtered normalised data 
nmr_data <- read.table(filtered_data_file,h=T, sep=",")

# read in manifest that links metabolite data to study ID and timepoint
manifest <- read.table(manifest_file,h=T, sep=",")

# read in top 50 
clarity_results <- read.table(paste0(results_dir, "clarity/20251114/clarity_results.txt"), sep="\t", h=T)

# read in LMM results
lmm_out <- read.table(paste0(lmm_dir, "01_nmr_results_2026_01_02.csv"), sep=",", h=T)

# feature metadata (ids to match data)
nmr_feature_metadata <- as.data.frame(read_delim(metadata_file,delim="\t"))

##############################################################################################################################
### PREP DATA #######################################################################################################
##############################################################################################################################

# add feature annotations to lmm results
lmm_out <- merge(lmm_out,nmr_feature_metadata,by.x="metabolite_name", by.y="feature_names")

##############################################################################################################################
### COMPARE WITH LMM RESULTS #############################################################################################
##############################################################################################################################

# add p-value rank to LMM output
lmm_out <- lmm_out[order(lmm_out$P_mixed),]
lmm_out$p.rank <- rank(lmm_out$P_mixed)

# merge results
all_results <- merge(clarity_results,lmm_out, by = "metabolite_name", all.y=T)
dim(all_results)
names(all_results)[1] <- "feature.name"
names(all_results)[2] <- "persistence.p.value"
names(all_results)[3] <- "persistence.value"

# order by persistence p
all_results <- all_results[order(all_results$persistence.p.value),]

# persistence p<0.05
all_results$clarity_flag <- 0
all_results <- all_results %>% mutate(`clarity_flag` = replace(`clarity_flag`, all_results$persistence.p.value < 0.05 , 1))
table(all_results$clarity_flag)

all_results$lmm_flag <- 0
all_results <- all_results %>% mutate(`lmm_flag` = replace(`lmm_flag`, all_results$p_adjusted < 0.05 , 1))
table(all_results$lmm_flag)

all_results$combo_flag <- 0
all_results <- all_results %>% mutate(`combo_flag` = replace(`combo_flag`, all_results$p_adjusted < 0.05 & all_results$persistence.p.value < 0.05, 1))
table(all_results$combo_flag)

table(all_results$clarity_flag,all_results$lmm_flag)

persistence_vars_05 <- all_results[which(all_results$persistence.p.value < 0.05),c("feature.name")]
persistence_vars_01 <- all_results[which(all_results$persistence.p.value < 0.01),c("feature.name")]

# look for overlap in signal
print("metabolites with clarity p < 0.05 and flagged in LMM:")
nrow(all_results[which(all_results$p_adjusted <0.05 & all_results$persistence.p.value <0.05),])
all_results[which(all_results$p_adjusted <0.05 & all_results$persistence.p.value <0.05),]

print("metabolites in with persistence p<0.05 and LMM unadj P < 0.05:")
nrow(all_results[which(all_results$P_mixed <0.05 & all_results$persistence.p.value <0.05),])
all_results[which(all_results$P_mixed <0.05 & all_results$persistence.p.value <0.05),]

print("metabolites in with persistence p<0.01 and LMM unadj P < 0.05:")
length(which(all_results$persistence.p.value <0.01 & all_results$P_mixed <0.05))
all_results$raw.label[which(all_results$persistence.p.value <0.01 & all_results$P_mixed <0.05)]

print("metabolites with clarity p < 0.05 and not flagged in LMM:")
nrow(all_results[which(all_results$p_adjusted >0.05 & all_results$persistence.p.value <0.05),])
all_results[which(all_results$p_adjusted >0.05 & all_results$persistence.p.value <0.05),]

print("metabolites with clarity p < 0.05 and not flagged in LMM:")
nrow(all_results[which(all_results$P_mixed >0.05 & all_results$persistence.p.value <0.05),])
all_results[which(all_results$P_mixed >0.05 & all_results$persistence.p.value <0.05),]

##############################################################################################################################
### SAVE OUT #################################################################################################################
##############################################################################################################################

columns_to_omit <- c("p.rank")
write.table(all_results[,!(names(all_results) %in% columns_to_omit)], file=paste0(results_dir, "combined/all_results.txt"),row.names=F,quote=F,sep="\t")
            
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
