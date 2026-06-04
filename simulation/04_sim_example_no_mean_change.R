#------------------- Set up session -------------------
source('simulations/00_sim_funs.R')
source('simulations/00_lmm_funs.R')

library(Clarity)
library(readr)
library(gridExtra)

library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
options(ggplot2.discrete.colour= c("#E69F00","#56B4E9","#CC79A7", "#009E73", "#F0E442", "#0072B2",
                                   "#D55E00","#999999", "#000000"))
#------------------- Set parameters -------------------

d <- 30

n <-  200
nU <- 2
nC <- 4
nX <- 200

#------------------- Simulate data using sim_data function -------------------
set.seed(1234)
out <- sim_data(
  n = n, nU = nU, nC = nC, nX = nX, clusters = 7, rho_time = 0.8,
  sigma_x = 0.0375, sigma_c = 0.12, mean_effect = 0
)

#------------------- Extract feature matrices and prepare for analysis -------------------

X1 <- out$X1  
X2 <- out$X2
C <- out$UC[,paste0('C',1:nC)]
anomalies <- out$anomalies


all_data <- rbind(X1, X2)

#------------------- LMM analysis -------------------
lmm_input <- cbind(all_data, 'Timepoint' = rep(c(1,2), each = n), rbind(C,C), studyid = rep(1:n, times = 2))
lmm <- metab.mixed.model(as.data.frame(lmm_input), columns_with_metabolite_data = paste0('X', 1:nX), covariables = paste0('C',1:3))

#extract and adjust p-vals
lmm_p_vals <- lmm$P_mixed
lmm_p_vals <- p.adjust(lmm_p_vals, method = c('BH'))

#------------------- CLARITY analysis -------------------
#Prepare covariance matrices
scaled <- scale(all_data)

X1 <- scaled[1:n,]
X2 <- scaled[(n+1):(2*n),]
X1_cov <- cov(X1)
X2_cov <- cov(X2)

X1_cov <- c_fixdiagonals(X1_cov)
X2_cov <- c_fixdiagonals(X2_cov)

#Do CLARITY pipeline
X1_scan_cov  <- Clarity_Scan(X1_cov, kmax = d,verbose=FALSE)
## Rotate the baseline to best-fit the endpoint
X2_pred_cov <- Clarity_Predict(X2_cov, X1_scan_cov)


pers_cov <- Clarity_Persistence(X2_pred_cov)


covt <- function(mat){
  return(cov(t(mat)))
}

compare <-  Clarity_Compare(X1_scan_cov, D = t(X1), Dnew = t(X2),
                         nbs = 200, distfn = covt, H0="scale") # structure
pvals=compare$pvals[,7]


#------------------- CLARITY pipeline plot -------------------


# CLARITY pipeline plot
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

#Set limits so things are on the same scale
cov_max <- max((abs(range(cbind(X1_cov, X2_cov)))))

X1_res <- X1_scan_cov$scan[[7]]$Yresid
X2_res <- X2_pred_cov$scan[[7]]$Yresid
res_max <- max(abs(cbind(X1_res, X2_res)))


zlim=c(0,1)*max((abs(range(cbind(range(Clarity_Persistence(X1_scan_cov)), 
                                 range(Clarity_Persistence(X2_pred_cov)))))))
zlim=sqrt(zlim)


mypar=c(5,2,2,2)
myline=2
png("figures/Figure_5_no_mean_change.png",
    height=600,width=1200)
par(mfrow=c(2,3))

Clarity_Chart(X1_cov,scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(0,1)*cov_max)

mtext("(a) Raw Baseline Covariance",side=1,line=myline,
      padj = 1.5, cex = 1.1)

Clarity_Chart(X1_scan_cov$scan[[7]]$Yresid,scalefun=I,cex.axis=0,
              las=2,mar=mypar,
              zlim=c(-1,1)*res_max,
              cols=my_palette)
mtext("(b) Baseline Residuals",side=1, line=myline,
      padj = 1.5, cex = 1.2)


plot(X1_scan_cov,rotate=TRUE,cex.axis.X = 1,
     cex.axis.Y=0,las=1,mar=mypar,zlim=zlim)
abline(v=7)
mtext("k=", line=1,side = 1, adj = -0.05, cex = 0.8)
mtext("(c) Baseline Persistences", side = 1,line=myline,
      padj = 2, cex = 1)


Clarity_Chart(X2_cov,scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(0,1)*cov_max)
a1range=range(which(anomalies==1)) + c(-0.5,0.5)
a2range=range(which(anomalies==2))+ c(-0.5,0.5)
a3range=range(which(anomalies==3))+ c(-0.5,0.5)
abline(h=c(a2range[1],a2range[2], a3range[1], a3range[2], a1range[1], a1range[2]), lty = 2)
abline(v=c(a2range[1],a2range[2], a3range[1], a3range[2], a1range[1], a1range[2]), lty = 2)
rect(xleft = c(a1range[1],a2range[1], a3range[1]), xright = c(a1range[2],a2range[2], a3range[2]),
     ybottom = c(a1range[1],a2range[1], a3range[1]), ytop = c(a1range[2],a2range[2], a3range[2]), lwd = 2)
mtext("(d) Raw Endpoint Covariance", side = 1,line=myline,
      padj = 2, cex = 1)


Clarity_Chart(X2_pred_cov$scan[[7]]$Yresid, scalefun=I,cex.axis=0,las=2,mar=mypar,
              zlim=c(-1,1)*res_max,
              cols=my_palette)
abline(h=c(a2range[1],a2range[2], a3range[1], a3range[2], a1range[1], a1range[2]), lty = 2)
abline(v=c(a2range[1],a2range[2], a3range[1], a3range[2], a1range[1], a1range[2]), lty = 2)
rect(xleft = c(a1range[1],a2range[1], a3range[1]), xright = c(a1range[2],a2range[2],  a3range[2]),
     ybottom = c(a1range[1],a2range[1], a3range[1]), ytop = c(a1range[2],a2range[2],  a3range[2]), lwd = 2)

mtext("(e) Endpoint Residuals",side=1, padj = 2, line=myline,cex = 1)

plot(X2_pred_cov,rotate=TRUE,cex.axis.X = 1,cex.axis.Y=0,las=1,
     mar=mypar,zlim=zlim)
abline(v=7)
abline(h=(c(a2range[1],a2range[2], a3range[1], a3range[2], a1range[1], a1range[2])),  lwd = 2, col = 'dark grey')
mtext("k=", line=1,side = 1, adj = -0.05, cex = 0.8)
mtext("(f) Endpoint Persistences", side = 1, padj = 2, line=myline,cex = 1)

dev.off()

#------------------- p-val plot -------------------


clarity_p_vals <- compare$pvals[,7]

status_labels <- c("Unaffected","Whole", "Part", "Remainder")

statuses <- status_labels[anomalies + 1]

plot_data <- tibble(Metabs = colnames(X1), 
                    LMM = -log10(lmm_p_vals),
                    CLARITY = -log10(clarity_p_vals),
                    Anomalies = as.factor(statuses))


testanomalies=as.factor(c(rep(0,21),rep(1,4),rep(2,2),rep(3,2)))
testrocs=sapply(1:100,function(rep){
  set.seed(rep)
  test_pred=rdirichlet(29,c(1,1,1,1))
  colnames(test_pred)=c("0","1","2","3")
  auc_comb <- as.numeric(pROC::multiclass.roc(testanomalies, test_pred)$auc)
  return(auc_comb)
})


#save p-values

write.table(plot_data, 'simulations/results/example_pvals_no_mean_change.csv')

plot_data=tibble(read.table("simulations/results/example_pvals_no_mean_change.csv"))
plot_data$Anomalies[plot_data$Anomalies=="Part"]<- "c2, Part"
plot_data$Anomalies[plot_data$Anomalies=="Remainder"]<- "c3, Remainder"
plot_data$Anomalies[plot_data$Anomalies=="Whole"]<- "c1, Whole"
plot_data$Anomalies=as.factor(plot_data$Anomalies)

## table of plot_data means by Anomalies class
sapply(status_labels, function(x) {
  mean(plot_data$LMM[plot_data$Anomalies == x])
})


legend_labels <- status_labels


png("figures/Figure_6_no_mean_change.png",
    height=600,width=600)
ggplot(plot_data, aes(LMM, CLARITY, colour = Anomalies)) +
  geom_point() + 
  geom_vline(xintercept = -log10(0.05), linetype = 'longdash', alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 'longdash', alpha = 0.5) +
  scale_color_discrete() +
  xlim(min(plot_data$CLARITY), max(plot_data$CLARITY))+
  theme_bw(base_size = 18) +
  labs(x=expression(-log[10]*p~"from LMM"),y=expression(-log[10]*p~"from CLARITY"))+
#  theme(
#    legend.position = c(0.95, 0.05),
#    legend.justification = c("right", "bottom"),
#    legend.background = element_rect(fill = alpha("white", 0.6))
#  ) + 
  ggtitle("Comparison of p-values with no mean change")

dev.off()

