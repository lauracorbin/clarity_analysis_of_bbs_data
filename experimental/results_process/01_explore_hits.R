# code to explore clarity results 

##############################################################################################################################
### SET UP ###################################################################################################################
##############################################################################################################################

# clear working space
rm(list = ls())

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
library(grid)
library(viridis)

## generate new feature labels
source("01a_feature_labels.R")

## read in function file
source("00_functions.R")

## capture output
sink(file=paste0(log_output_dir, "01_explore_hits.log"),split=TRUE)

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

# read in results
pval <- read.table(paste0(results_dir, "clarity/20251114/endpoint_from_baseline_pvals.csv"), sep=",", h=T)
persistence <- read.table(paste0(results_dir, "clarity/20251114/endpoint_from_baseline_pers.csv"), sep=",", h=T)

# read in LMM results
lmm_out <- read.table(paste0(lmm_dir, "01_nmr_results_2026_01_02.csv"), sep=",", h=T)

# contemporary feature annotations
feature_annot <- read.table(paste0(input_dir,"feature_annotation_and_exclusions.csv"), sep=",", h=T)

# feature metadata (ids to match data)
nmr_feature_metadata <- as.data.frame(read_delim(metadata_file,delim="\t"))

##############################################################################################################################
### PREP DATA #######################################################################################################
##############################################################################################################################

# split by timepoint
baseline_ids <- manifest$sample_id[manifest$Visit.Identifier == "Baseline"]
endpoint_ids <- manifest$sample_id[manifest$Visit.Identifier == "36 months"]

baseline_data_z <- nmr_data[which(row.names(nmr_data) %in% baseline_ids),]
endpoint_data_z <- nmr_data[which(row.names(nmr_data) %in% endpoint_ids),]

# add feature annotations to lmm results
lmm_out <- merge(lmm_out,nmr_feature_metadata,by.x="metabolite_name", by.y="feature_names")

##############################################################################################################################
### MERGE WITH LMM RESULTS #############################################################################################
##############################################################################################################################

# add p-value rank to LMM output
lmm_out <- lmm_out[order(lmm_out$P_mixed),]
lmm_out$p.rank <- rank(lmm_out$P_mixed)

# extract results for k=25
pval_k25 <- as.data.frame(pval[,c("V10")])
pval_k25$metabolite_name <- row.names(pval)
names(pval_k25)[1] <- "pval_k25"

pers_k25 <- as.data.frame(persistence[,c("V10")])
pers_k25$metabolite_name <- row.names(persistence)
names(pers_k25)[1] <- "pers_k25"

clarity_results <- merge(pval_k25, pers_k25, by="metabolite_name")

# merge results
all_results <- merge(clarity_results,lmm_out, by = "metabolite_name", all.y=T)
names(all_results)[1] <- "feature.name"
names(all_results)[2] <- "persistence.p.value"
names(all_results)[3] <- "persistence.value"

# order by persistence p
all_results <- all_results[order(all_results$persistence.p.value),]

# persistence p<0.05
all_results[which(all_results$persistence.p.value < 0.05),]
dim(all_results[which(all_results$persistence.p.value < 0.05),])

persistence_vars_05 <- all_results[which(all_results$persistence.p.value < 0.05),c("feature.name")]
persistence_vars_01 <- all_results[which(all_results$persistence.p.value < 0.01),c("feature.name")]

print("Metabolites with persistence p<0.05")
all_results$raw.label[all_results$persistence.p.value <0.05]

print("Metabolites with persistence p<0.01")
all_results$raw.label[all_results$persistence.p.value <0.01]


##############################################################################################################################
### CONSIDER CLASSES, ETC #############################################################################################
##############################################################################################################################

persistence_vars_05_data <- all_results[which(all_results$persistence.p.value < 0.05),]

table(persistence_vars_05_data$class)
table(persistence_vars_05_data$subclass)
table(persistence_vars_05_data$label)

##############################################################################################################################
### bubble plots, etc #############################################################################################
##############################################################################################################################

# AIM: take top hits and generate cor matrix with associated p-values then display in bubble plot

# create ID to match bubbleplot template
metabolite_names <- as.data.frame(names(baseline_data_z))
names(metabolite_names)[1] <- "bbs.feature.id"

# calculate cor matrix with rcorr() to also generate p-values
baseline_cormat <- rcorr(as.matrix(baseline_data_z), type=c("pearson"))

endpoint_cormat <- rcorr(as.matrix(endpoint_data_z), type=c("pearson"))

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flattened_baseline_cormat <- flattenCorrMatrix(baseline_cormat$r,baseline_cormat$P)
flattened_endpoint_cormat <- flattenCorrMatrix(endpoint_cormat$r,endpoint_cormat$P)

print(sum(abs(flattened_baseline_cormat$cor)))
print(sum(abs(flattened_endpoint_cormat$cor)))

# add BCAAs to plot list as examples of non-persistent
bcaa <- c("ile","leu","val","totalbcaa")

plot_list_combined <- c(persistence_vars_05,bcaa)

# list to save graphs
plot_list = vector("list", length = 2*(length(plot_list_combined)))
plot_count = 0

# capture correlations metric columns
all_results$total_r_baseline <- 0
all_results$total_r_endpoint <- 0
baseline_r_col <- which(names(all_results) == "total_r_baseline")
endpoint_r_col <- which(names(all_results) == "total_r_endpoint")
feature.name.col <- which(names(all_results) == "feature.name")

# list to cpature histograms of corr
plot_list_baseline <- list()
plot_list_endpoint <- list()

# loop over metabolites
for(i in 1:length(all_results$feature.name)){
#for(i in 1:2){
  metab <- all_results[i,feature.name.col]
  print(metab)
  flattened_baseline_cormat_hit <- flattened_baseline_cormat[flattened_baseline_cormat$row == metab | flattened_baseline_cormat$column == metab ,]
  dim(flattened_baseline_cormat_hit)
  flattened_endpoint_cormat_hit <- flattened_endpoint_cormat[flattened_endpoint_cormat$row == metab | flattened_endpoint_cormat$column == metab ,]
  dim(flattened_endpoint_cormat_hit)

  flattened_baseline_cormat_hit$metabolite <- NA
  for (j in 1:nrow(flattened_baseline_cormat_hit)) {
    if (flattened_baseline_cormat_hit[j,1] == metab) {
      flattened_baseline_cormat_hit[j,5] = as.character(flattened_baseline_cormat_hit[j,2])}
    if (flattened_baseline_cormat_hit[j,1] != metab) {
      flattened_baseline_cormat_hit[j,5] = as.character(flattened_baseline_cormat_hit[j,1])}
    # replace pval = 0 with very small number
    if (flattened_baseline_cormat_hit[j,4] == 0) {
      flattened_baseline_cormat_hit[j,4] = 1.0e-17
      }
    }

  flattened_endpoint_cormat_hit$metabolite <- NA
  for (j in 1:nrow(flattened_endpoint_cormat_hit)) {
    if (flattened_endpoint_cormat_hit[j,1] == metab) {
      flattened_endpoint_cormat_hit[j,5] = as.character(flattened_endpoint_cormat_hit[j,2])}
    if (flattened_endpoint_cormat_hit[j,1] != metab) {
      flattened_endpoint_cormat_hit[j,5] = as.character(flattened_endpoint_cormat_hit[j,1])}
    # replace pval = 0 with very small number
    if (flattened_endpoint_cormat_hit[j,4] == 0) {
      flattened_endpoint_cormat_hit[j,4] = 1.0e-17
      }
    }

  # look at overall corr
  print(sum(abs(flattened_baseline_cormat_hit$cor)))
  all_results[i,baseline_r_col] <- sum(abs(flattened_baseline_cormat_hit$cor))
  print(sum(abs(flattened_endpoint_cormat_hit$cor)))
  all_results[i,endpoint_r_col] <- sum(abs(flattened_endpoint_cormat_hit$cor))

  # histogram of r
  plot_list_baseline[[metab]] <- ggplot(data = flattened_baseline_cormat_hit, aes(x = cor)) +
                                  geom_histogram() +
                                  xlim(-1,1) +
                                  ggtitle(metab)

  plot_list_endpoint[[metab]] <- ggplot(data = flattened_endpoint_cormat_hit, aes(x = cor)) +
                                  geom_histogram() +
                                  xlim(-1,1) +
                                  ggtitle(metab)

  # graph hits and selected others
  if (metab %in% plot_list_combined){
    plot_count = plot_count + 1
    # merge with annotations
    flattened_baseline_cormat_hit$ckb_id <- final_annotation$ckb_id[match(flattened_baseline_cormat_hit$metabolite,final_annotation$bbs.id)]
    flattened_endpoint_cormat_hit$ckb_id <- final_annotation$ckb_id[match(flattened_endpoint_cormat_hit$metabolite,final_annotation$bbs.id)]

    # build input data
    # merge with template but only keep the non-derived that have been analysed
    myData_baseline <- merge(flattened_baseline_cormat_hit,ntngale250,all.x = TRUE,by="ckb_id",sort=FALSE)
    names(myData_baseline)[4] <- "estimate"
    names(myData_baseline)[5] <- "pval"
    myData_baseline$negLog10P <- -log10(myData_baseline$pval)

    myData_endpoint <- merge(flattened_endpoint_cormat_hit,ntngale250,all.x = TRUE,by="ckb_id",sort=FALSE)
    names(myData_endpoint)[4] <- "estimate"
    names(myData_endpoint)[5] <- "pval"
    myData_endpoint$negLog10P <- -log10(myData_endpoint$pval)

    # plot baseline
    filename <- paste0(fig_dir,plot_count,"_",metab,"_bubbleplot_baseline.pdf")
    print("baseline format data")
    gridData <- MyformatData(myData_baseline, colorValue="estimate", sizeValue = "negLog10P",
                           nightingale = TRUE)
    print("baseline create plot input")
    treeInput <- MymultiPlotInput(colorList=gridData$colorList,
                                sizeList=gridData$sizeList,
                                nightingale=TRUE, legendHeight=8)
    print("baseline create heatmap list")
    treeList <- bubbleHeatmapList(treeInput)
    print("baseline create plot")
    metabTree <- MymetabFigurePlot(treeList)
    cairo_pdf(filename,width=11, height =8)
    grid.draw(metabTree)
    dev.off()

    # plot endpoint
    filename <- paste0(fig_dir,plot_count,"_",metab,"_bubbleplot_endpoint.pdf")
    gridData <- MyformatData(myData_endpoint, colorValue="estimate", sizeValue = "negLog10P",
                             nightingale = TRUE)
    treeInput <- MymultiPlotInput(colorList=gridData$colorList,
                                sizeList=gridData$sizeList,
                                nightingale=TRUE, legendHeight=8)
    treeList <- bubbleHeatmapList(treeInput)
    metabTree <- MymetabFigurePlot(treeList)
    cairo_pdf(filename,width=11, height =8)
    grid.draw(metabTree)
    dev.off()
  }
}

# export histograms
filename = paste0(fig_dir,"Online_resource_4.pdf")
outputfig <- ggpubr::ggarrange(plotlist=plot_list_baseline,
                               ncol=2, nrow=3,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)

filename = paste0(fig_dir,"Online_resource_5.pdf")
outputfig <- ggpubr::ggarrange(plotlist=plot_list_endpoint,
                               ncol=2, nrow=3,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)

##############################################################################################################################
### hist of r diff #############################################################################################
##############################################################################################################################

clarity_results$total_r_baseline <- all_results$total_r_baseline[match(clarity_results$metabolite_name,all_results$feature.name)]
clarity_results$total_r_endpoint <- all_results$total_r_endpoint[match(clarity_results$metabolite_name,all_results$feature.name)]
clarity_results$total_r_diff <- clarity_results$total_r_endpoint - clarity_results$total_r_baseline

clarity_results$clarity_flag <- "not persistent"
clarity_results <- clarity_results %>% mutate(`clarity_flag` = replace(`clarity_flag`, clarity_results$pval_k25 < 0.05 , "persistent"))

# update theme/colours
theme_set(theme_bw())
options(ggplot2.discrete.colour= c("#E69F00","#56B4E9","#CC79A7", "#009E73", "#F0E442", "#0072B2",
                                   "#D55E00","#999999", "#000000"))

filename = paste0(fig_dir,"SuppFig2_r_hist.pdf")
pdf(filename)
ggplot(data = clarity_results, aes(x = total_r_diff, fill = as.factor(clarity_flag))) + 
  geom_histogram() +
  theme_set(theme_bw()) + 
  guides(fill=guide_legend(title="Metabolite group")) +
  xlab("Difference in total correlation\n (post-treatment - pre-treatment)") 
  #scale_fill_viridis(discrete=T,option="turbo",begin=0.2, end=0.9)
dev.off()

##############################################################################################################################
### SAVE OUT #################################################################################################################
##############################################################################################################################
  
columns_to_omit <- c("clarity_flag")
write.table(clarity_results[,!(names(clarity_results) %in% columns_to_omit)], paste0(results_dir,"clarity/20251114/clarity_results.txt"), sep="\t", row.names=F, quote=F)

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
