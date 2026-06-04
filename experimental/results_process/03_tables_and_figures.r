# code to generate tables and figures for manuscript

##############################################################################################################################
### SET UP ###################################################################################################################
##############################################################################################################################

# clear working space
rm(list = ls())

## read in parameters 
source("parameter_files/00_parameters.R")

## load libraries 
library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(viridis)

## capture output
sink(file=paste0(log_output_dir, "03_tables_and_figures.log"),split=TRUE)

## record todays date
today = Sys.Date()
today = gsub("-","_",today)
today

##############################################################################################################################
### READ IN RESULTS ##########################################################################################################
##############################################################################################################################

# combined lmm and clarity output
all_results <- read.table(file=paste0(results_dir, "combined/all_results.txt"),h=T, sep="\t")

# edge pvals
edge_p <- read.table(file=paste0(results_dir, "clarity/20251114/edge_pvals.csv"),h=T, sep=",")

##############################################################################################################################
### GENERATE GROUP VAR #######################################################################################################
##############################################################################################################################

# add new labels that are collapsed classes
all_results$Group <- "Other"
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Amino acids", "Amino acids"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Apolipoproteins", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Lipoprotein subclasses", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Lipoprotein particle concentrations", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Lipoprotein particle sizes", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Other lipids", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Glycerides and phospholipids", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Phospholipids", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Total lipids", "Lipid and lipoprotein molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Fatty acids", "Fatty acids"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Cholesterol", "Cholesterol molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Free cholesterol", "Cholesterol molecules"))
all_results <- all_results %>% mutate(`Group` = replace(`Group`, `class` == "Cholesteryl esters", "Cholesterol molecules"))
table(all_results$class, all_results$Group)


##############################################################################################################################
### CLASS DISTRIBUTIONS ######################################################################################################
##############################################################################################################################

## identify associated features for each analysis
all_results$lmm_flag = 0
all_results <- all_results %>% mutate(`lmm_flag` = replace(`lmm_flag`, all_results$p_adjusted < 0.05 , 1))

all_results$clarity_flag = 0
all_results <- all_results %>% mutate(`clarity_flag` = replace(`clarity_flag`, all_results$persistence.p.value < 0.05 , 1))

# separate tested and associated
clarity_associated <- all_results[all_results$clarity_flag == 1,]
clarity_tested <- all_results[!is.na(all_results$persistence.p.value),]
lmm_associated_169 <- all_results[all_results$lmm_flag == 1 & !is.na(all_results$persistence.p.value),]

# extract numbers for stack for tested
table_for_all <- as.data.frame(table(clarity_tested$class))
# add %
table_for_all$class_pct <- (table_for_all$Freq / sum(table_for_all$Freq))*100

# extract numbers for stack for clarity assoc
table_for_assoc <- as.data.frame(table(clarity_associated$class))
# add %
table_for_assoc$class_pct <- (table_for_assoc$Freq / sum(table_for_assoc$Freq))*100

# extract numbers for stack for lmm assoc
table_for_associated_lmm <- as.data.frame(table(lmm_associated_169$class))
# add %
table_for_associated_lmm$class_pct <- (table_for_associated_lmm$Freq / sum(table_for_associated_lmm$Freq))*100

# add label for all versus assoc
table_for_all$outcome <- paste0("All tested","\n","metabolites")
table_for_associated_lmm$outcome <- paste0("LMM","\n", "associated")
table_for_assoc$outcome <- paste0("CLARITY", "\n", "associated")

# join to make graph input
stacked_bar_data_combo1 <- rbind(table_for_all,table_for_associated_lmm)
stacked_bar_data_combo <- rbind(stacked_bar_data_combo1,table_for_assoc )

# update var labels
names(stacked_bar_data_combo)[1] <- "Class"
names(stacked_bar_data_combo)[2] <- "freq_overall"
names(stacked_bar_data_combo)[3] <- "Percentage.in.class"

# re-order outcome factor 
stacked_bar_data_combo$outcome <- as.factor(stacked_bar_data_combo$outcome)
stacked_bar_data_combo <- stacked_bar_data_combo %>%
  mutate(outcome = fct_relevel(outcome, levels(outcome)[c(1,3,2)]))

# plot stacked bar
filename <- paste0(results_dir,"figures/SuppFig1_class_distribution.pdf")
pdf(filename)

ggplot(stacked_bar_data_combo, aes(fill=Class,y=Percentage.in.class,x=outcome)) +
  geom_bar(position="stack",stat="identity",colour="black") +
  theme_set(theme_bw()) +
  theme(axis.title.x = element_blank()) +
  labs(y="Percentage in class") +
  scale_fill_viridis(discrete=T,option="turbo")
invisible(dev.off())


# join to make table output
names(table_for_all)[2:4] <- paste0("tested_",names(table_for_all)[2:4])
names(table_for_assoc)[2:4] <- paste0("clarity_",names(table_for_assoc)[2:4])
names(table_for_associated_lmm)[2:4] <- paste0("lmm_",names(table_for_assoc)[2:4])

table_out1 <- merge(table_for_all, table_for_associated_lmm, by="Var1", all=T)
table_out <- merge(table_out1, table_for_assoc, by="Var1", all=T)
table_out <- table_out[,c(1,2,3,5,6,8,9)]

write.table(table_out,paste0(results_dir,"tables/class_distribution.txt"), row.names=F, quote=F, sep="\t")



##############################################################################################################################
### COMPARE RESULTS ##########################################################################################################
##############################################################################################################################

# update theme/colours
theme_set(theme_bw())
options(ggplot2.discrete.colour= c("#E69F00","#56B4E9","#CC79A7", "#009E73", "#F0E442", "#0072B2",
                                   "#D55E00","#999999", "#000000"))

filename <- paste0(results_dir,"figures/Figure4_scatter_compare.pdf")
pdf(filename)
ggplot(all_results, aes(x=-log10(p_adjusted), y=-log10(persistence.p.value), colour=Group, shape = Group)) + 
  geom_point() +
  geom_text_repel(data=subset(all_results, persistence.p.value < 0.05),
            aes(-log10(p_adjusted),-log10(persistence.p.value),label=feature.name),show.legend = FALSE) +
  geom_text_repel(data=subset(all_results, p_adjusted < 0.05 & persistence.p.value > 0.05),
                  aes(-log10(p_adjusted),-log10(persistence.p.value),label=feature.name),show.legend = FALSE) +
  theme_set(theme_bw()) +
  theme(legend.position = "bottom") + 
  guides(colour=guide_legend(nrow=5)) +
  labs(x=expression(-log[10]*p~"from LMM"),y=expression(-log[10]*p~"from CLARITY")) +
  #scale_colour_viridis(discrete=T,option="turbo") +
  scale_shape_manual(values = 1:length(unique(all_results$Group))) +
  geom_hline(yintercept = -log10(0.05), lty = 'longdash', alpha=0.5) + 
  geom_vline(xintercept = -log10(0.05), lty = 'longdash', alpha=0.5)
invisible(dev.off())

##############################################################################################################################
### cor of r & persistence #############################################################################################
##############################################################################################################################

plot_list <- list()

plot_list[[1]] <- ggplot(data = all_results, aes(x = persistence.value, y = -log10(persistence.p.value),color=Group, shape = Group)) + 
  labs(x="Persistence value ",y=expression(-log[10]*p~"")) +
  geom_point() +
  theme(legend.position = "bottom") + 
  theme_set(theme_bw()) +
  geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  theme(legend.position = "none") +
  #scale_colour_viridis(discrete=T,option="turbo") +
  scale_shape_manual(values = 1:length(unique(all_results$Group)))

plot_list[[2]] <- ggplot(data = all_results, aes(x = total_r_baseline, y = -log10(persistence.p.value),color=Group, shape = Group)) + 
  labs(x="Total correlation at baseline",y=expression(-log[10]*p~"")) +
  geom_point() +
  theme(legend.position = "bottom") + 
  theme_set(theme_bw()) +
  geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  theme(legend.position = "none") +
  #scale_colour_viridis(discrete=T,option="turbo") +
  scale_shape_manual(values = 1:length(unique(all_results$Group))) 

plot_list[[3]] <- ggplot(data = all_results, aes(x = total_r_baseline, y = persistence.value,color=Group, shape = Group)) + 
  labs(x="Total correlation at baseline",y="Persistence value") +
  geom_point() +
  theme(legend.position = "bottom") + 
  theme_set(theme_bw()) +
  guides(colour=guide_legend(nrow=6)) +
  #scale_colour_viridis(discrete=T,option="turbo") +
  scale_shape_manual(values = 1:length(unique(all_results$Group))) 


filename = paste0(fig_dir,"SuppFig3_three_plots.pdf")
outputfig <- ggpubr::ggarrange(plotlist=plot_list,
                               ncol=1, nrow=3,
                               common.legend = T,
                               legend="right",
                               heights = 2,
                               widths = 3)
ggpubr::ggexport(outputfig, filename=filename)

##############################################################################################################################
### make version of all results for suppl table ##############################################################################
##############################################################################################################################

# add column - BH adjusted persistence p
all_results[,'persistence.p.value.adjusted'] <- p.adjust(all_results[,'persistence.p.value'], method = 'BH') 

# restrict to columns needed for table
pub_results <- all_results[,c("feature.name", "label", "class", "subclass",
                              "beta_mixed",  "lci_mixed", "uci_mixed", "P_mixed", "p_adjusted",
                              "R2_model", "R2_timepoint", "R2_model_conditional", "R2_timepoint_conditional",   
                              "persistence.value","persistence.p.value",'persistence.p.value.adjusted', 
                              "total_r_baseline", "total_r_endpoint", "total_r_diff"
)]

# rename columns
names(pub_results) <- c("metabolite_short_name", "metabolite_full_name", "metabolite_class", "metabolite_subclass",
                        "lmm_beta", "lmm_lci", "lmm_uci", "lmm_p", "lmm_adj.p",
                        "lmm_total_R2", "lmm_timepoint_R2", "lmm_total_cond.R2", "lmm_timepoint_cond.R2",
                        "clarity_persistence.value.at.k25", "clarity_persistence.p.value.at.k25", "clarity_persistence.p.value.at.k25.adjusted",
                        "total_r_baseline", "total_r_endopint","total_r_diff")

# Identify numeric columns
is_numeric_col <- sapply(pub_results, is.numeric)
is_numeric_col[8] <- FALSE # remove lmm_p from list as that should be scientific

# Apply round to numeric columns except lmm_p (make that scientific)
pub_results[is_numeric_col] <- lapply(pub_results[is_numeric_col], round, digits = 3)

# change lmm_p to scientific notation
pub_results$lmm_p <- formatC(pub_results$lmm_p, format = "e", digits = 3)

# save out
write.table(pub_results,file = paste0(table_dir,"Supplementary_table_1.txt"), quote=F, row.names=F, sep="\t")

##############################################################################################################################
### extract edge pvals for editing bubbleplots ###############################################################################
##############################################################################################################################

# apoa1
apoa1 <- as.data.frame(edge_p[,c("apoa1")])
names(apoa1)[1] <- "pval"
apoa1$metabs <- rownames(edge_p)
apoa1_sig <- apoa1[apoa1$pval < 0.05,]
apoa1_sig <- apoa1_sig[order(apoa1_sig$pval),]
apoa1_sig

#valine
valine <- as.data.frame(edge_p[,c("val")])
names(valine)[1] <- "pval"
valine$metabs <- rownames(edge_p)
valine_sig <- valine[valine$pval < 0.05,]
valine_sig <- valine_sig[order(valine_sig$pval),]
valine_sig

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