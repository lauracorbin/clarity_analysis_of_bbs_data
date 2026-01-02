library(dplyr)
library(Clarity)
library(readr)

source("parameters/00_parameters.R")


#---------- Read in data ----------
nmr_data <- read.table(nmr_file,
                       sep='\t',h=T)
manifest <- read.table(manifest_file,
                       h=T, sep=",")


features <- read.csv(feature_file, header = F)
keep_features <- features$V1

filtered <- read.table(filtered_file,
                       h=T, sep=",")


#---------- Prepare data for analysis ----------
dim(nmr_data)
nmr_data <- nmr_data[which(row.names(nmr_data) %in% filtered$sample_id),]
dim(nmr_data)

nmr_data <- nmr_data %>% mutate_all(~(scale(.) %>% as.vector))

baseline_ids <- manifest$sample_id[manifest$Visit.Identifier == "Baseline"]
endpoint_ids <- manifest$sample_id[manifest$Visit.Identifier == "36 months"]


baseline_data <- nmr_data[which(row.names(nmr_data) %in% baseline_ids),]
endpoint_data <- nmr_data[which(row.names(nmr_data) %in% endpoint_ids),]

#---------- Do clarity analysis ----------

# create inputs  
baseline_covar <- cov(baseline_data, use="complete.obs")
baseline_covar <- baseline_covar[keep_features, keep_features]


endpoint_covar <- cov(endpoint_data, use = "complete.obs")
endpoint_covar <- endpoint_covar[keep_features, keep_features]


baseline_data <- baseline_data[, keep_features]
endpoint_data <- endpoint_data[ , keep_features]

# drop any rows with na vals 
baseline_data <- baseline_data[rowSums(is.na(baseline_data)) ==0,]

endpoint_data <- endpoint_data[rowSums(is.na(endpoint_data)) == 0,]

# transpose for Clarity_Compare
baseline_datat <- t(baseline_data)
endpoint_datat <- t(endpoint_data)
covt <- function(mat){
  return(cov(t(mat)))
}

# compute p-values
base_scan <- Clarity_Scan(baseline_covar, kmax = 60)

compare_from_base <- Clarity_Compare(base_scan, D = baseline_datat, Dnew = endpoint_datat,
                                     nbs = 1000, distfn = covt, H0="scale")

comp_from_base_edge <- Clarity_Compare(base_scan, D = baseline_datat, Dnew = endpoint_datat,
                                       nbs = 1000, distfn = covt, H0="scale", k = 25)

##
hide<-function(){ # for comparison purposes to see what our old residuals were
tmp_scale=Clarity_Compare(base_scan,Ynew=endpoint_covar, H0="scale")
tmp_val=Clarity_Compare(base_scan,Ynew=endpoint_covar)
pers_scale=Clarity_Persistence(tmp_scale)
pers_val=Clarity_Persistence(tmp_val)
test=as.matrix(read.csv(paste0(output_dir,"persistences/endpoint_from_baseline_pers.csv"),header=T,row.names=1))
range(pers_val-test) # [1] -4.560796e-13  3.452794e-13
range(pers_scale-test) # [1] -26.01833695  -0.00013293
range(pers_val/mean(pers_val) - pers_scale/mean(pers_scale)) # much smaller though still non-zero
Clarity_Chart(pers_val, scalefun=I, cex.axis=0.7, las=2,
              mar=c(5,10,2,2), main="Endpoint from Baseline - Val")
Clarity_Chart(pers_scale, scalefun=I, cex.axis=0.7, las=2,
              mar=c(5,10,2,2), main="Endpoint from Baseline - Scale")
}
#---------- Save data ---------- 

# save bootstrap samples
samples <- compare_from_base$Ylist
for (i in 1:1000) {
  dir <- paste0('clinical_analysis/bootstrapped_samples/', i, '/')
  dir.create(dir, showWarnings = FALSE)
  
  current_sample <- samples[[i]]
  
  write.table(current_sample$Y1, paste0(dir, "Y1.csv"), sep=",",
              col.names=T, row.names=T)
  write.table(current_sample$Y2, paste0(dir, "Y2.csv"), sep=",",
              col.names=T, row.names=T)
  write.table(current_sample$Y3, paste0(dir, "Y3.csv"), sep=",",
              col.names=T, row.names=T)
}

write.table(compare_from_base$pvals,paste0(output_dir,"pvals/endpoint_from_baseline_pvals.csv"),sep=",",col.names=T,row.names=T)

write.table(comp_from_base_edge$pvals,paste0(output_dir,"pvals/edge_pvals.csv"),sep=",",col.names=T,row.names=T)
