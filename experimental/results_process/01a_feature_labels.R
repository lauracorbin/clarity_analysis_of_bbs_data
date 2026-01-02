# code to explore apob association with other metabs and weight loss, etc

##############################################################################################################################
### SET UP ###################################################################################################################
##############################################################################################################################

## read in parameters 
source("parameter_files/00_parameters.R")

# feature annotations file
feature_annot <- read.table(paste0(input_dir,"feature_annotation_and_exclusions.csv"), sep=",", h=T)

# make feature ID for bubble plot merge
feature_annot$ckb_id1 <- tolower(feature_annot$CSV.column.name)
feature_annot$ckb_id <- paste0("nmr_",feature_annot$ckb_id1)
feature_annot$bbs.id <- gsub("_","",feature_annot$ckb_id1)

# edit those that still do not match
rown <- which(feature_annot$CSV.column.name == "ApoB_by_ApoA1")
feature_annot[rown,10] <- "apobapoa1"
feature_annot[rown,9] <- "nmr_apob_apoa1"

rown <- which(feature_annot$CSV.column.name == "PUFA_by_MUFA")
feature_annot[rown,10] <- "pufamufa"
feature_annot[rown,9] <- "nmr_pufa_mufa" # don't think this exists in template file

rown <- which(feature_annot$CSV.column.name == "Omega_6_by_Omega_3")
feature_annot[rown,10] <- "omega6omega3"
feature_annot[rown,9] <- "nmr_omega6_omega3" # don't think this exists in template file

# matched
matched <- feature_annot[feature_annot$ckb_id %in% ntngale249$ckb_id,]
dim(matched)

# not matched
unmatched <- feature_annot[!feature_annot$ckb_id %in% ntngale249$ckb_id,]
dim(unmatched)

# add to matched by cross matching metabolite names to those in bubble package data
#myData <- merge_template(ntngale249, "ckb_id")

# try matching on full name
unmatched$ckb_id <- ntngale249$ckb_id[match(tolower(unmatched$Biomarker.name),tolower(ntngale249$description))]
length(which(!unmatched$ckb_id %in% ntngale249$ckb_id))
new_match1 <- unmatched[unmatched$ckb_id %in% ntngale249$ckb_id,]
dim(new_match1)
matched <- rbind(matched,new_match1)
dim(matched)

unmatched <- unmatched[!unmatched$ckb_id %in% ntngale249$ckb_id,]
dim(unmatched)

unmatched_mydata <- ntngale249[!ntngale249$ckb_id %in% matched$ckb_id,]

# replace pct with per for % vars
unmatched$ckb_id <- gsub("pct", "per", unmatched$ckb_id1)
unmatched$ckb_id <- paste0("nmr_",unmatched$ckb_id)

new_match2 <- unmatched[unmatched$ckb_id %in% ntngale249$ckb_id,]
dim(new_match2)
matched <- rbind(matched,new_match2)
dim(matched)

unmatched <- unmatched[!unmatched$ckb_id %in% ntngale249$ckb_id,]
dim(unmatched)

unmatched_mydata <- ntngale249[!ntngale249$ckb_id %in% matched$ckb_id,]

# edit those that still do not match
rown <- which(unmatched$CSV.column.name == "DHA_pct")
unmatched[rown,10] <- "dhapct"
unmatched[rown,9] <- "nmr_dha_fa"

rown <- which(unmatched$CSV.column.name == "Omega_3_pct")
unmatched[rown,10] <- "omega3pct"
unmatched[rown,9] <- "nmr_faw3_fa"

rown <- which(unmatched$CSV.column.name == "Omega_6_pct")
unmatched[rown,10] <- "omega6pct"
unmatched[rown,9] <- "nmr_faw6_fa"

rown <- which(unmatched$CSV.column.name == "LA_pct")
unmatched[rown,10] <- "lapct"
unmatched[rown,9] <- "nmr_la_fa"

rown <- which(unmatched$CSV.column.name == "MUFA_pct")
unmatched[rown,10] <- "mufapct"
unmatched[rown,9] <- "nmr_mufa_fa"

rown <- which(unmatched$CSV.column.name == "PUFA_pct")
unmatched[rown,10] <- "pufapct"
unmatched[rown,9] <- "nmr_pufa_fa"

rown <- which(unmatched$CSV.column.name == "SFA_pct")
unmatched[rown,10] <- "sfapct"
unmatched[rown,9] <- "nmr_sfa_fa"

rown <- which(unmatched$CSV.column.name == "TG_by_PG")
unmatched[rown,10] <- "tgpg"
unmatched[rown,9] <- "nmr_tg_pg"

rown <- which(unmatched$CSV.column.name == "PUFA_by_MUFA")
unmatched[rown,10] <- "pufamufa"
unmatched[rown,9] <- "nmr_21"

rown <- which(unmatched$CSV.column.name == "Omega_6_by_Omega_3")
unmatched[rown,10] <- "omega6omega3"
unmatched[rown,9] <- "nmr_22"

new_match3 <- unmatched[unmatched$ckb_id %in% ntngale249$ckb_id,]
dim(new_match3)
matched <- rbind(matched,new_match3)
dim(matched)

unmatched <- unmatched[!unmatched$ckb_id %in% ntngale249$ckb_id,]
dim(unmatched)
unmatched$ckb_id <- "nmr_glycerol" # only one left - glycerol

#unmatched_mydata <- ntngale249[!ntngale249$ckb_id %in% matched$ckb_id,]

# put back together
final_annotation <- matched
final_annotation <- rbind(matched,unmatched)

# add in glycerol
new_row <- data.frame(
  ckb_id = "nmr_glycerol",
  ukbb_id = 99999,
  description = "Glycerol",
  rowName = "Glycerol",
  colName = "Glycolysis",
  plotGroup = "Glycolysis"
)

ntngale250 <- rbind(ntngale249,new_row)

