### large graph simulations
### ECOLI (Gaussian)

rm(list=ls())

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
} 

library(bnlearn) # processing of true graph

library(snow) # for parallel computing
library(Rmpi) # for parallel computing
library(parallel) # for parallel computing

source("gaussCItwd.R")
source("gaussMItest.R")
source("getSuff.R")
source("Tu_MissingValuePC.R")
source("Tu_CITest.R")
source("Tu_myfixes.R")
source("streetparty.R")


# load network
ecoli <- readRDS("ecoli70.rds")

# true adjacencies in moral graph
moralgraph <- moral(ecoli)
trueMoral <- amat(moralgraph)


sim <- function(., network, trueMoral, n, mech) {
  ## network    a bn.fit object
  ## trueMoral  a matrix; true adjacencies in moral graph
  ## n          sample size
  ## mech       type of missingness mechanism: "MCAR", "MAR" or "MNAR"
  
  ### generate data ############################################################
  dat <- rbn(network, n=n)
  
  ### analyse full data ########################################################
  suff <- getSuff(dat, test="gaussCItest", adaptDF=FALSE)
  full_res <- pc(suffStat=suff, indepTest=gaussCItest, alpha=0.05,
                 labels=colnames(dat), skel.method="stable", maj.rule=TRUE,
                 solve.confl=TRUE)
  
  p <- 12

  ### lose data ################################################################
  # missing values are only in the subset of variables used in the other
  # analysis
  dat_small <- dat[ ,c("b1191", "cchB", "eutG", "fixC", "ibpB", "sucA", "tnaA",
                       "yceP", "yfaD", "ygbD", "ygcE", "yjbO")]
  
  if (mech=="MCAR") {
    # randomly delete 18% of data points
    
    miss <- sample(n*p, round(n*p*0.18), replace=FALSE)
    miss <- matrix(1:(n*p) %in% miss, ncol=p)
    dat_small[miss] <- NA
    
  } else if (mech=="MAR") {
    # choose 2 times 3 variables containing missing values
    # desired pattern: for each observation, exactly one of these 3 is missing,
    #   and the probability depends on the other two
    
    v1 <- sample(p, 3)
    v2 <- sample(setdiff(1:p, v1), 3)
  
    # missinness patterns
    pat <- matrix(c(0,1,1,
                    1,0,1,
                    1,1,0),
                  byrow=TRUE, nrow=3)
    
    amp1 <- ampute(data=dat_small[ ,v1], prop=0.999, patterns=pat,
                   freq=c(1/3,1/3,1/3))
    dat_small[ ,v1] <- amp1$amp
    amp2 <- ampute(data=dat_small[ ,v2], prop=0.999, patterns=pat,
                   freq=c(1/3,1/3,1/3))
    dat_small[ ,v2] <- amp2$amp
    
    # in order to obtain a missingness proportion of 18% overall, choose one
    # further variable and delete values completely at random
    v3 <- sample(setdiff(1:p, c(v1,v2)), 1)
    dat_small[sample(n,round(n*p*0.013)),v3] <- NA
    
  } else if (mech=="MNAR") {
    # 12/8*18% of observations contain missing values in 5 selected variables
    # missingness depends on one of these 5, so MNAR
    
    del_index <- order(dat_small$fixC, decreasing=TRUE)[1:(round(n*0.18*12/8))]
    dat_small[del_index,
        c("fixC","eutG","ibpB","sucA","tnaA","yceP","yfaD","ygcE")]<- NA
  }
  
  dat[ ,c("b1191", "cchB", "eutG", "fixC", "ibpB", "sucA", "tnaA",
                       "yceP", "yfaD", "ygbD", "ygcE", "yjbO")] <- dat_small
  
  p <- ncol(dat)
  
  ### list-wise deletion #######################################################
  datl <- dat[complete.cases(dat), ]
  if (nrow(datl) < 10) {
    lwd_res <- NA
  } else {
    
    suffl <- getSuff(datl, test="gaussCItest")
    
    lwd_res <- tryCatch( pc(suffStat=suffl, indepTest=gaussCItest, alpha=0.05,
                            labels=colnames(dat), skel.method="stable",
                            maj.rule=TRUE, solve.confl=TRUE),
                         error=function(e){paste("Error")} )
    
    if ( identical(lwd_res, "Error") ) {lwd_res <- NA}
  }
  
  ### test-wise deletion #######################################################
  twd_res <- tryCatch( pc(suffStat=dat, indepTest=gaussCItwd, alpha=0.05,
                          labels=colnames(dat), skel.method="stable",
                          maj.rule=TRUE, solve.confl=TRUE),
                       error=function(e){paste("Error")} )
  
  if ( identical(twd_res, "Error") ) {twd_res <- NA}
  
  ### Tu 1 method (permutation) ################################################
  tu1_res <- tryCatch( mvpc(suffStat=list(data=dat), indepTest=gaussCItest.td,
                            corrMethod=gaussCItest.permc, alpha=0.05, p=p,
                            labels=colnames(dat), skel.method="stable",
                            maj.rule=TRUE, solve.confl=TRUE),
                       error=function(e){paste("Error")} )
  
  if ( identical(tu1_res, "Error") ) {tu1_res <- NA}
  
  ### Tu 2 method (density ratio weighting) ####################################
  tu2_res <- tryCatch( mvpc(suffStat=list(data=dat), indepTest=gaussCItest.td,
                            corrMethod=gaussCItest.drw, alpha=0.05, p=p,
                            labels=colnames(dat), skel.method="stable",
                            maj.rule=TRUE, solve.confl=TRUE),
                       error=function(e){paste("Error")} )
  
  if ( identical(tu2_res, "Error") ) {tu2_res <- NA}
  
  ### oracle multiple imputation ###############################################
  omi <- tryCatch( mice(dat, m=10, method="norm", print=FALSE,
                        predictorMatrix=trueMoral),
                   error=function(e){paste("Error")} )
  if ( identical(omi, "Error") ){
    omi_res <- NA
  } else {
    omi_res <- pc(suffStat=getSuff(omi, test="gaussMItest"),
                  indepTest=gaussMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### standard multiple imputation #############################################
  mi <- tryCatch( mice(dat, m=10, method="norm", print=FALSE),
                  error=function(e){paste("Error")} )
  if ( identical(mi, "Error") ){
    mi_res <- NA
  } else {
    mi_res <- pc(suffStat=getSuff(mi, test="gaussMItest"),
                 indepTest=gaussMItest,
                 alpha=0.05, labels=colnames(dat), skel.method="stable",
                 maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### multiple imputation using random forests #################################
  rf <- tryCatch( mice(dat, m=10, method="rf", print=FALSE, ntree=100),
                  error=function(e){paste("Error")} )
  if ( identical(rf, "Error") ){
    rf_res <- NA
  } else {
    rf_res <- pc(suffStat=getSuff(rf, test="gaussMItest"),
                 indepTest=gaussMItest,
                 alpha=0.05, labels=colnames(dat), skel.method="stable",
                 maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### multiple imputation using random forests - CALIBER #######################
  rfc <- tryCatch( mice(dat, m=10, method="rfcont", print=FALSE, ntree=100),
                   error=function(e){paste("Error")} )
  if ( identical(rfc, "Error") ){
    rfc_res <- NA
  } else {
    rfc_res <- pc(suffStat=getSuff(rfc, test="gaussMItest"),
                  indepTest=gaussMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### two-step #################################################################
  skel1 <- tryCatch( skeleton(suffStat=dat, indepTest=gaussCItwd, alpha=0.2,
                              labels=colnames(dat), method="stable"),
                       error=function(e){paste("Error")} )
  if ( identical(skel1, "Error") ) {
     pre1_res <- NA
  } else {
    skelm1 <- streetparty(as(skel1, "matrix"))
    mi1 <- tryCatch( mice(dat, m=10, method="norm", print=FALSE,
                          predictorMatrix=skelm1),
                    error=function(e){paste("Error")} )
    if ( identical(mi1, "Error") ){
      pre1_res <- NA
    } else {
      pre1_res <- pc(suffStat=getSuff(mi1, test="gaussMItest"),
                   indepTest=gaussMItest,
                   alpha=0.05, labels=colnames(dat), skel.method="stable",
                   maj.rule=TRUE, solve.confl=TRUE)
    }
  }
  
  ### two-step marginal only ###################################################
  skel2 <- tryCatch( skeleton(suffStat=dat, indepTest=gaussCItwd, alpha=0.2,
                              labels=colnames(dat), method="stable", m.max=0),
                       error=function(e){paste("Error")} )
  if ( identical(skel2, "Error") ) {
     pre2_res <- NA
  } else {
    skelm2 <- streetparty(as(skel2, "matrix"))
    mi2 <- tryCatch( mice(dat, m=10, method="norm", print=FALSE,
                          predictorMatrix=skelm2),
                    error=function(e){paste("Error")} )
    if ( identical(mi2, "Error") ){
      pre2_res <- NA
    } else {
      pre2_res <- pc(suffStat=getSuff(mi2, test="gaussMItest"),
                   indepTest=gaussMItest,
                   alpha=0.05, labels=colnames(dat), skel.method="stable",
                   maj.rule=TRUE, solve.confl=TRUE)
    }
  }
  
  ### two-step marginal only ignoring neighbours of neighbours #################
  skel3 <- tryCatch( skeleton(suffStat=dat, indepTest=gaussCItwd, alpha=0.2,
                              labels=colnames(dat), method="stable", m.max=0),
                       error=function(e){paste("Error")} )
  if ( identical(skel3, "Error") ) {
     pre3_res <- NA
  } else {
    skelm3 <- as(skel3, "matrix")
    mi3 <- tryCatch( mice(dat, m=10, method="norm", print=FALSE,
                          predictorMatrix=skelm3),
                    error=function(e){paste("Error")} )
    if ( identical(mi3, "Error") ){
      pre3_res <- NA
    } else {
      pre3_res <- pc(suffStat=getSuff(mi3, test="gaussMItest"),
                   indepTest=gaussMItest,
                   alpha=0.05, labels=colnames(dat), skel.method="stable",
                   maj.rule=TRUE, solve.confl=TRUE)
    }
  }

  ### mean imputation ##########################################################
  means <- apply(dat, 2, mean, na.rm=TRUE)
  dat2 <- dat
  for (i in 1:ncol(dat2)) {
    dat2[is.na(dat2[ ,i]), i] <- means[i] 
  }
  
  suffm <- getSuff(dat2, test="gaussCItest")
  mean_res <- pc(suffStat=suffm, indepTest=gaussCItest, alpha=0.05,
                 labels=colnames(dat), skel.method="stable", maj.rule=TRUE,
                 solve.confl=TRUE)
  
  ### return results ###########################################################
  res <- list(full_res, lwd_res, twd_res, tu1_res, tu2_res, omi_res, mi_res,
              rf_res, rfc_res, pre1_res, pre2_res, pre3_res, mean_res)
  res <- lapply(res, function(r){
    if(!inherits(r, "pcAlgo")){return(r)}else{return(as(r, "amat"))}
  })
  return(res)
}

seed <- 28359
nrep <- 1000

cl<-makeCluster(239, type="MPI")
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(bnlearn)) # for data generation
clusterCall(cl, function() library(pcalg)) # for pc
clusterCall(cl, function() library(abind)) # cannot be found although it should be loaded by pcalg
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterCall(cl, function() library(CALIBERrfimpute)) # multiple imputation
clusterCall(cl, function() library(ks)) # for Tu et al. method
clusterCall(cl, function() library(weights)) # for Tu et al. method
clusterCall(cl, function() library(DescTools)) # for Tu et al. method
clusterExport(cl, c("gaussCItwd", "gaussMItest", "getSuff", "zStatMI",
                    "log.q1pm", 
                    "detection_prt_m", "get_m_ind", "get_prt_R_ind",
                    "graph2gaps", "mvpc", "skeleton2", "common.neighbor",
                    "compute_weights", "compute.weights.continuous",
                    "cond.PermC", "f_R", "f_weights", "gaussCItest.drw",
                    "gaussCItest.permc", "gaussCItest.td", "gaussCItest.td.ref",
                    "get_ind_r_xys", "get_ind_weights", "get_logi_formula",
                    "get_logidata", "get_prt_i", "get_prt_m_xys",
                    "get_prt_R_ind", "get_rw_pair", "indx_test_wise_deletion",
                    "is.in_prt_m", "iscorr", "kde.weights", "kdrw.weights",
                    "weights_check_table", "perm", "test_wise_deletion",
                    "test_wise_deletion_w", "streetparty"))



start <- Sys.time()
ECOLI_large_100_MCAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=100, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_100_MAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=100, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_100_MNAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=100, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_1000_MCAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=1000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_1000_MAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=1000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_1000_MNAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=1000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_5000_MCAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=5000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_5000_MAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=5000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

start <- Sys.time()
ECOLI_large_5000_MNAR <- parLapply(cl, 1:nrep, sim, network=ecoli,
                                 trueMoral=trueMoral, n=5000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("ECOLI", ls())],
     file="res_ECOLI_large.RData")

stopCluster(cl)
