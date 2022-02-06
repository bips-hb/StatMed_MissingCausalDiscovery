### small graph simulations
### MAGIC (Gaussian)

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


# load network
magic <- readRDS("magic-niab.rds")

# true adjacencies in moral graph
moralgraph <- moral(magic)
trueMoral <- amat(moralgraph)
trueMoral <-
  trueMoral[c("MIL", "G1217", "G257", "G2208", "G1338", "G524", "G1945"),
            c("MIL", "G1217", "G257", "G2208", "G1338", "G524", "G1945")]


sim <- function(., network, trueMoral, n, mech) {
  ## network    a bn.fit object
  ## trueMoral  a matrix; true adjacencies in moral graph
  ## n          sample size
  ## mech       type of missingness mechanism: "MCAR", "MAR" or "MNAR"
  
  ### generate data ############################################################
  dat <- rbn(network, n=n)
  dat <- dat[ ,c("MIL", "G1217", "G257", "G2208", "G1338", "G524", "G1945")]
  
  ### analyse full data ########################################################
  suff <- getSuff(dat, test="gaussCItest", adaptDF=FALSE)
  full_res <- pc(suffStat=suff, indepTest=gaussCItest, alpha=0.05,
                 labels=colnames(dat), skel.method="stable", maj.rule=TRUE,
                 solve.confl=TRUE)
  
  p <- ncol(dat)

  ### lose data ################################################################
  if (mech=="MCAR") {
    # randomly delete 18% of data points
    
    miss <- sample(n*p, round(n*p*0.18), replace=FALSE)
    miss <- matrix(1:(n*p) %in% miss, ncol=p)
    dat[miss] <- NA
    
  } else if (mech=="MAR") {
    # choose 3 variables containing missing values
    # desired pattern: for each observation, exactly one of these 3 is missing,
    #   and the probability depends on the other two
    
    v1 <- sample(p, 3)
  
    # missinness patterns
    pat <- matrix(c(0,1,1,
                    1,0,1,
                    1,1,0),
                  byrow=TRUE, nrow=3)
    
    amp <- ampute(data=dat[ ,v1], prop=0.999, patterns=pat,
                  freq=c(1/3,1/3,1/3))
    dat[ ,v1] <- amp$amp
    
    # in order to obtain a missingness proportion of 18% overall, choose one
    # further variable and delete values completely at random
    v2 <- sample(setdiff(1:p, v1), 1)
    dat[sample(n,round(n*p*0.037)),v2] <- NA
    
  } else if (mech=="MNAR") {
    # 7/5*18% of observations contain missing values in 5 selected variables
    # missingness depends on one of these 5, so MNAR
    
    del_index <- order(dat$MIL, decreasing=TRUE)[1:(round(n*0.18*7/5))]
    dat[del_index, c("MIL","G524","G1338","G1217","G2208")] <- NA
  }
  
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
                            corrMethod=gaussCItest.permc, alpha=0.05, p=8,
                            labels=colnames(dat), skel.method="stable",
                            maj.rule=TRUE, solve.confl=TRUE),
                       error=function(e){paste("Error")} )
  
  if ( identical(tu1_res, "Error") ) {tu1_res <- NA}
  
  ### Tu 2 method (density ratio weighting) ####################################
  tu2_res <- tryCatch( mvpc(suffStat=list(data=dat), indepTest=gaussCItest.td,
                            corrMethod=gaussCItest.drw, alpha=0.05, p=8,
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

  ### mean imputation ##########################################################
  means <- apply(dat, 2, mean, na.rm=TRUE)
  dat2 <- dat
  for (i in 1:ncol(dat2)) {
    dat2[is.na(dat2[ ,i]), i] <- means[i] 
  }
  
  suffm <- getSuff(dat2, test="gaussCItest")
  mode_res <- pc(suffStat=suffm, indepTest=gaussCItest, alpha=0.05,
                 labels=colnames(dat), skel.method="stable", maj.rule=TRUE,
                 solve.confl=TRUE)
  
  ### return results ###########################################################
  res <- list(full_res, lwd_res, twd_res, tu1_res, tu2_res, omi_res, mi_res,
              rf_res, rfc_res, mode_res)
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
                    "test_wise_deletion_w"))



start <- Sys.time()
MAGIC_100_MCAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=100, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_100_MAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=100, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_100_MNAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=100, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_1000_MCAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=1000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_1000_MAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=1000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_1000_MNAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=1000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_5000_MCAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=5000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_5000_MAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=5000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

start <- Sys.time()
MAGIC_5000_MNAR <- parLapply(cl, 1:nrep, sim, network=magic,
                                 trueMoral=trueMoral, n=5000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MAGIC", ls())],
     file="res_MAGIC.RData")

stopCluster(cl)
