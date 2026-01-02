source("parameters/00_parameters.R")

library(readr)
library(Clarity)

#--------- Load and prep data ----------

baseline_covar <- as.matrix(
  read.table(paste0(project_dir, 'data/intermediate_data/baseline_z_covar_matrix.csv'),
             header=T, sep = ","))

endpoint_covar <- as.matrix(
  read.table(paste0(project_dir, 'data/intermediate_data/endpoint_z_covar_matrix.csv'),
             header=T, sep=","))

get_threshold <- function(mat, n) {
  vals <- sort(as.numeric(abs(mat)), decreasing = TRUE)
  threshold <- vals[n]
  return(threshold)
}

features <- read.csv(paste0(project_dir,"data/input_data/feature_inclusion_info_2025_10_08.txt"),
                     header = F)

keep_features <- features$V1


baseline_covar <- baseline_covar[keep_features,keep_features]

endpoint_covar <- endpoint_covar[keep_features,keep_features]
baseline_covar <- scale(baseline_covar)

baseline_covar <- (baseline_covar + t(baseline_covar))/2

endpoint_covar <- scale(endpoint_covar)

endpoint_covar <- (endpoint_covar + t(endpoint_covar))/2


really_as_matrix<-function(A){
  ret=matrix(as.numeric(A),dim(A)[1],dim(A)[2]) # Coerce to standard matrix
  rownames(ret)=rownames(A)
  colnames(ret)=colnames(A)
  return(ret)
}
kmax <- 60

baseline_scan <- Clarity_Scan(baseline_covar, kmax=kmax)
endpoint_pred <- Clarity_Predict(endpoint_covar, baseline_scan)


endpoint_res <- Clarity_Extract(endpoint_pred, summary = I, k = 25)
endpoint_res <- really_as_matrix(endpoint_res)


endpoint_scan <- Clarity_Scan(endpoint_covar, kmax = kmax)
baseline_pred <- Clarity_Predict(baseline_covar, endpoint_scan)

baseline_res <- Clarity_Extract(baseline_pred, summary = I, k = 25)
baseline_res <- really_as_matrix(baseline_res)


#--------- Make Figure 5 ----------
get_threshold <- function(mat, n) {
  vals <- sort(as.numeric(abs(mat)), decreasing = TRUE)
  threshold <- vals[n]
  return(threshold)
}

base_thresh <- get_threshold(baseline_res, 6000)
end_thresh <- get_threshold(endpoint_res, 6000)

max <- max((abs(range(cbind(baseline_covar, endpoint_covar)))))

my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

pvals=read.csv(paste0(output_dir,"pvals/endpoint_from_baseline_pvals.csv"),header=T,row.names=1)
edge_pvals=read.csv(paste0(output_dir,"pvals/edge_pvals.csv"),header=T,row.names=1)
c_pred=Clarity_Persistence(endpoint_pred)
zlim_pers=c(0,sqrt(max(c_pred)))

mypar=c(5,5,2,1)
png("figures/Figure_5.png",height=600,width=1200)
par(mfrow=c(2,3))

Clarity_Chart(baseline_covar,scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(-1,1)*max,
              cols=my_palette)
#col.axis.Y=grouping_colors,col.axis.X=grouping_colors)
mtext("(a) Raw Baseline Covariance",side=1, padj = 1.5, cex = 1.5)

Clarity_Chart(baseline_scan$scan[[25]]$Yresid,scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(-1,1)*sqrt(max(abs(range(endpoint_res)))),
#              zlim=c(-1,1)*max(abs(range(baseline_scan$scan[[25]]$Yresid))),
              #col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
              cols=my_palette)
mtext("(b) Baseline Residuals",side=1, padj = 1.5, cex = 1.5)


plot(baseline_scan,rotate=TRUE,cex.axis.X = 1,zlim=zlim_pers,
     cex.axis.Y=0,las=1,mar=mypar)
mtext("k=", line=1,side = 1, adj = -0.05, cex = 0.8)
#col.axis.Y=grouping_colors)
abline(v=25)
mtext("(c) Baseline Persistences", side = 1, padj = 2, cex = 1.5)

Clarity_Chart(endpoint_covar,scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(-1,1)*max,
              cols=my_palette)
mtext("(d) Raw Endpoint Covariance", side = 1, padj = 2, cex = 1.5)

Clarity_Chart(endpoint_res,scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(-1,1)*max(abs(range(endpoint_res))),signif = edge_pvals<0.05,
              #col.axis.Y=grouping_colors,col.axis.X=grouping_colors,
              cols=my_palette)
mtext("(e) Endpoint Residuals",side=1, padj = 2, cex = 1.5)

plot(endpoint_pred,rotate=TRUE,cex.axis.X = 1,zlim=zlim_pers,
     cex.axis.Y=0,las=1,mar=mypar,signif = pvals<0.05)
mtext("k=", line=1,side = 1, adj = -0.05, cex = 0.8)
abline(v=25)
mtext("(f) Endpoint Persistences", side = 1, padj = 2, cex = 1.5)
dev.off()

