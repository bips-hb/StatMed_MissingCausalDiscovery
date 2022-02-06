library(naniar)
library(mice)
library(CALIBERrfimpute)
library(bnlearn)
library(tpc)
library(micd)
library(ggplot2)

### load raw data ##############################################################
dat <- readRDS("IDEFICS_data.RDS") # we are not allowed to share the data

### visualise missing data pattern #############################################
gg <- vis_miss(dat)
gg +  theme( plot.margin = margin(0, 2, 0, 0.5, "cm") )

cairo_ps(file="VisMiss.eps",
    width=6, height=4)
  gg +  theme( plot.margin = margin(0, 2, 0, 0.5, "cm") )
dev.off()

### specify tiers and additional forbidden edges for tPC #######################
tiers <- rep(c(1,2,3,4,5,6,7,8), c(5,1,1,1,3,2,3,7))
forbEdges <- matrix(0, ncol=23, nrow=23)
forbEdges[1:4, 5] <- TRUE
forbEdges[1:13, 12] <- TRUE

### list-wise deletion #########################################################
datl <- dat[complete.cases(dat), ]
res_lwd <- tpc(suffStat=datl, indepTest=mixCItest, labels=colnames(dat),
               alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
plot(res_lwd)
as(res_lwd, "matrix")[c(7,18,19,21,23)]

### test-wise deletion #########################################################
res_twd <- tpc(suffStat=dat, indepTest=mixCItwd, labels=colnames(dat),
               alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
plot(res_twd)
Rgraphviz::plot(res_twd@graph)
as()

### parametric multiple imputation #############################################
set.seed(28359)
form_mi <- make.formulas.saturated(dat, d=2)
mi <- mice(dat, m=100, formulas=form_mi,
           defaultMethod=c("norm","logreg","polyreg","polr"))
# failed while trying to impute bage (Error in solve.default(xtx + diag(pen)) : 
# system is computationally singular: reciprocal condition number = 1.18701e-18)

### parametric main effects multiple imputation ################################
set.seed(30659)
miw <- mice(dat, m=100, defaultMethod=c("norm","logreg","polyreg","polr"))
suffmiw <- complete(miw, action="all")
res_miw <- tpc(suffStat=suffmiw, indepTest=mixMItest, labels=colnames(dat),
               alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)

## random forests multiple imputation ##########################################
set.seed(30659)
rf <- mice(dat, m=100, method="rf", ntree=100)
suffrf <- complete(rf, action="all")
res_rf <- tpc(suffStat=suffrf, indepTest=mixMItest, labels=colnames(dat),
              alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
plot(res_rf)

## CALIBER multiple imputation #################################################
set.seed(30659)
rfc <- mice(dat, m=100, ntree=100,
            defaultMethod=c("rfcont","rfcat","rfcat","polr"))
suffrfc <- complete(rfc, action="all")
res_rfc <- tpc(suffStat=suffrfc, indepTest=mixMItest, labels=colnames(dat),
               alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
plot(res_rfc)

### mean/mode imputation #######################################################
modes <- lapply(dat, function(i) {
  if (is.factor(i)) {
    wh <- which.max(table(i))
    m <- levels(i)[wh]
    return(m)
  } else {
   return(mean(i, na.rm=TRUE)) 
  }
})
dat2 <- dat
for (i in 1:ncol(dat2)) {
  dat2[is.na(dat2[ ,i]), i] <- modes[[i]]
}

### save results ###############################################################
res_mode <- tpc(suffStat=dat2, indepTest=mixCItest, labels=colnames(dat),
                alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)

save(res_lwd, res_twd, res_rf, res_rfc, res_miw, res_mode,
     miw, rf, rfc,
     file="res.RData")

load("res.RData")

### calculate pairwise Hamming distances
T12 <- hamming(as.bn(res_lwd), as.bn(res_twd))
T13 <- hamming(as.bn(res_lwd), as.bn(res_miw))
T14 <- hamming(as.bn(res_lwd), as.bn(res_rf))
T15 <- hamming(as.bn(res_lwd), as.bn(res_rfc))
T16 <- hamming(as.bn(res_lwd), as.bn(res_mode))

T23 <- hamming(as.bn(res_twd), as.bn(res_miw))
T24 <- hamming(as.bn(res_twd), as.bn(res_rf))
T25 <- hamming(as.bn(res_twd), as.bn(res_rfc))
T26 <- hamming(as.bn(res_twd), as.bn(res_mode))

T34 <- hamming(as.bn(res_miw), as.bn(res_rf))
T35 <- hamming(as.bn(res_miw), as.bn(res_rfc))
T36 <- hamming(as.bn(res_miw), as.bn(res_mode))

T45 <- hamming(as.bn(res_rf), as.bn(res_rfc))
T46 <- hamming(as.bn(res_rf), as.bn(res_mode))

T56 <- hamming(as.bn(res_rfc), as.bn(res_mode))

Ta <- matrix(0, nrow=6, ncol=6)
colnames(Ta) <- c("lwd", "twd", "miw", "rf", "rfc", "mode")
rownames(Ta) <- c("lwd", "twd", "miw", "rf", "rfc", "mode")
Ta[1,2] <- Ta[2,1] <- T12
Ta[1,3] <- Ta[3,1] <- T13
Ta[1,4] <- Ta[4,1] <- T14
Ta[1,5] <- Ta[5,1] <- T15
Ta[1,6] <- Ta[6,1] <- T16
Ta[2,3] <- Ta[3,2] <- T23
Ta[2,4] <- Ta[4,2] <- T24
Ta[2,5] <- Ta[5,2] <- T25
Ta[2,6] <- Ta[6,2] <- T26
Ta[3,4] <- Ta[4,3] <- T34
Ta[3,5] <- Ta[5,3] <- T35
Ta[3,6] <- Ta[6,3] <- T36
Ta[4,5] <- Ta[5,4] <- T45
Ta[4,6] <- Ta[6,4] <- T46
Ta[5,6] <- Ta[6,5] <- T56

### regression diagnostics #####################################################

## randomly choose 10 imputations for diagnostics
set.seed(123)
ind <- sample(100, 10)
# ind <- c(31, 79, 51, 14, 67, 42, 50, 43, 97, 25)

col <- rep("transparent", 100)
col[ind] <- 1:10

## convergence plots
cairo_ps("plot_conv_miw.eps",
    width=10, height=13)
  plot(miw, col=col, layout=c(4,10))
dev.off()

cairo_ps("plot_conv_rf.eps",
    width=10, height=13)
  plot(rf, col=col, layout=c(4,10))
dev.off()

cairo_ps("plot_conv_rfc.eps",
    width=10, height=13)
  plot(rfc, col=col, layout=c(4,10))
dev.off()

## boxplots and density plots
long_miw <- complete(miw, action='long', include=TRUE)
long_miw <- long_miw[long_miw$.imp%in%c(0,ind), ]
long_miw$.imp <- rep(0:10, each=657)
sub_miw <- as.mids(long_miw)

cairo_ps("plot_box_miw.eps",
    width=10, height=13)
  bwplot(sub_miw, layout=c(3,5))
dev.off()

cairo_ps("plot_dens_miw.eps",
    width=10, height=13)
  densityplot(sub_miw, layout=c(3,5))
dev.off()

long_rf <- complete(rf, action='long', include=TRUE)
long_rf <- long_rf[long_rf$.imp%in%c(0,ind), ]
long_rf$.imp <- rep(0:10, each=657)
sub_rf <- as.mids(long_rf)

cairo_ps("plot_box_rf.eps",
    width=10, height=13)
  bwplot(sub_rf, layout=c(3,5))
dev.off()

cairo_ps("plot_dens_rf.eps",
    width=10, height=13)
  densityplot(sub_rf, layout=c(3,5))
dev.off()

long_rfc <- complete(rfc, action='long', include=TRUE)
long_rfc <- long_rfc[long_rfc$.imp%in%c(0,ind), ]
long_rfc$.imp <- rep(0:10, each=657)
sub_rfc <- as.mids(long_rfc)

cairo_ps("plot_box_rfc.eps",
    width=10, height=13)
  bwplot(sub_rfc, layout=c(3,5))
dev.off()

cairo_ps("plot_dens_rfc.eps",
    width=10, height=13)
  densityplot(sub_rfc, layout=c(3,5))
dev.off()


## plots for main paper
cairo_ps("main_box_wb_miw.eps",
    width=3, height=2.5)
bwplot(sub_miw, wb, ylim=c(18,55), main="parametric\nmain-effects")
dev.off()

cairo_ps("main_dens_wb_miw.eps",
    width=3, height=2.8)
densityplot(sub_miw, data=~wb, ylim=c(-0.01, 0.14))
dev.off()


cairo_ps("main_box_wb_rf.eps",
    width=3, height=2.5)
bwplot(sub_rf, wb, ylim=c(18,55), main="random forests\n(rf)")
dev.off()

cairo_ps("main_dens_wb_rf.eps",
    width=3, height=2.8)
densityplot(sub_rf, data=~wb, ylim=c(-0.01, 0.14))
dev.off()


cairo_ps("main_box_wb_rfc.eps",
    width=3, height=2.5)
bwplot(sub_rfc, wb, ylim=c(18,55), main="random forests\n(CALIBER)")
dev.off()

cairo_ps("main_dens_wb_rfc.eps",
    width=3, height=2.8)
densityplot(sub_rfc, data=~wb, ylim=c(-0.01, 0.14))
dev.off()
  



