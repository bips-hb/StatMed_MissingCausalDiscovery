### small graph simulations
### SACHS (categorial with 3 categories per variable)

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

source("disCItwd.R")
source("disMItest.R")
source("getSuff.R")
source("make_formulas_saturated.R")


# load network
sachs <- readRDS("sachs.rds")

# true adjacencies in moral graph
moralgraph <- moral(sachs)
trueMoral <- amat(moralgraph)


sim <- function(., network, trueMoral, n, mech) {
  ## network    a bn.fit object
  ## trueMoral  a matrix; true adjacencies in moral graph
  ## n          sample size
  ## mech       type of missingness mechanism: "MCAR", "MAR" or "MNAR"
  
  ### generate data ############################################################
  dat <- rbn(network, n=n)
  
  ### analyse full data ########################################################
  datf <- dat
  datf[] <- lapply(datf, factor)
  suff <- getSuff(datf, test="disCItest", adaptDF=FALSE)
  full_res <- pc(suffStat=suff, indepTest=disCItest_new, alpha=0.05,
                 labels=colnames(dat), skel.method="stable", maj.rule=TRUE,
                 solve.confl=TRUE)
  
  p <- ncol(dat)

  ### lose data ################################################################
  if (mech=="MCAR") {
    # randomly delete 18% of data points
    
    miss <- sample(n*p, round(n*p*0.18), replace=FALSE)
    miss <- matrix(1:(n*p) %in% miss, ncol=p)
    dat[miss] <- NA
    dat[] <- lapply(dat, factor)
    
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
    
    amp1 <- ampute(data=dat[ ,v1], prop=0.999, patterns=pat, freq=c(1/3,1/3,1/3))
    dat[ ,v1] <- amp1$amp
    amp2 <- ampute(data=dat[ ,v2], prop=0.999, patterns=pat, freq=c(1/3,1/3,1/3))
    dat[ ,v2] <- amp2$amp
    
  } else if (mech=="MNAR") {
    # 11/7*18% of observations contain missing values in 5 selected variables
    # missingness depends on one of these 5, so MNAR
    
    del_index <- order(dat$PKC, decreasing=TRUE)[1:(round(n*0.18*11/7))]
    dat[del_index, c("PKC","PKA","Raf","Mek","Erk","Akt","Jnk")] <- NA 
  }
  
  dat[] <- lapply(dat, factor)
  
  ### list-wise deletion #######################################################
  datl <- dat[complete.cases(dat), ]
  if (nrow(datl)<5) {
    lwd_res <- NA
  } else {
    
    datl[] <- lapply(datl, factor)
    suffl <- getSuff(datl, test="disCItest", adaptDF=FALSE)
    
    lwd_res <- tryCatch( pc(suffStat=suffl, indepTest=disCItest_new, alpha=0.05,
                            labels=colnames(dat), skel.method="stable",
                            maj.rule=TRUE, solve.confl=TRUE),
                         error=function(e){paste("Error")} )
    
    if ( identical(lwd_res, "Error") ) {lwd_res <- NA}
  }
  
  ### test-wise deletion #######################################################
  datt <- dat
  datt[] <- lapply(datt, factor)
  sufft <- getSuff(datt, test="disCItest", adaptDF=FALSE)
  
  twd_res <- tryCatch( pc(suffStat=sufft, indepTest=disCItwd, alpha=0.05,
                          labels=colnames(dat), skel.method="stable",
                          maj.rule=TRUE, solve.confl=TRUE),
                       error=function(e){paste("Error")} )
  
  if ( identical(twd_res, "Error") ) {twd_res <- NA}
  
  ### oracle multiple imputation ###############################################
  onelevel <- sapply(dat, nlevels)==1
  trueMoral[onelevel, ] <- 0
  trueMoral[ ,onelevel] <- 0
  form_omi <- make.formulas.saturated(dat, predictorMatrix=trueMoral, d=3)
  omi <- tryCatch( mice(dat, m=10, print=FALSE, formulas=form_omi),
                   error=function(e){paste("Error")} )
  if ( identical(omi, "Error") ){
    omi_res <- NA
  } else {
    omi_res <- pc(suffStat=getSuff(omi, test="disMItest"), indepTest=disMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### standard multiple imputation #############################################
  pred <- matrix(1, p, p)
  colnames(pred) <- rownames(pred) <- colnames(dat)
  diag(pred) <- 0
  pred[onelevel, ] <- 0
  pred[ ,onelevel] <- 0
  form_mi <- make.formulas.saturated(dat, predictorMatrix=pred, d=2)
  mi <- tryCatch( mice(dat, m=10, print=FALSE, formulas=form_mi),
                  error=function(e){paste("Error")} )
  if ( identical(mi, "Error") ){
    mi_res <- NA
  } else {
    mi_res <- pc(suffStat=getSuff(mi, test="disMItest"), indepTest=disMItest,
                 alpha=0.05, labels=colnames(dat), skel.method="stable",
                 maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### standard multiple imputation without interaction terms ###################
  miw <- tryCatch( mice(dat, m=10, print=FALSE),
                   error=function(e){paste("Error")} )
  if ( identical(miw, "Error") ){
    miw_res <- NA
  } else {
    miw_res <- pc(suffStat=getSuff(miw, test="disMItest"), indepTest=disMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### multiple imputation using random forests #################################
  rf <- tryCatch( mice(dat, m=10, method="rf", print=FALSE, ntree=100),
                  error=function(e){paste("Error")} )
  if ( identical(rf, "Error") ){
    rf_res <- NA
  } else {
    rf_res <- pc(suffStat=getSuff(rf, test="disMItest"), indepTest=disMItest,
                 alpha=0.05, labels=colnames(dat), skel.method="stable",
                 maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### multiple imputation using random forests - CALIBER #######################
  rfc <- tryCatch( mice(dat, m=10, method="rfcat", print=FALSE, ntree=100),
                   error=function(e){paste("Error")} )
  if ( identical(rfc, "Error") ){
    rfc_res <- NA
  } else {
    rfc_res <- pc(suffStat=getSuff(rfc, test="disMItest"), indepTest=disMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }

  ### mode imputation ##########################################################
  modes <- sapply(dat, function(i) {
    wh <- which.max(table(i))
    m <- levels(i)[wh]
    return(m)
  })
  
  dat2 <- dat
  for (i in 1:ncol(dat2)) {
    dat2[is.na(dat2[ ,i]), i] <- modes[[i]]
  }
  
  datm <- dat2
  datm[] <- lapply(datm, factor)
  suffm <- getSuff(datm, test="disCItest", adaptDF=FALSE)
  mode_res <- pc(suffStat=suffm, indepTest=disCItest_new, alpha=0.05,
                 labels=colnames(dat), skel.method="stable", maj.rule=TRUE,
                 solve.confl=TRUE)
  
  ### return results ###########################################################
  res <- list(full_res, lwd_res, twd_res, omi_res, mi_res, miw_res, rf_res,
              rfc_res, mode_res)
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
clusterCall(cl, function() library(pcalg)) # for gaussCItest
clusterCall(cl, function() library(abind)) # cannot be found although it should be loaded by pcalg
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterCall(cl, function() library(CALIBERrfimpute)) # multiple imputation
clusterExport(cl, c("disCItest_new", "gSquareDis_new", "disCItwd", "disMItest",
                    "getSuff", "make.formulas.saturated",
                    "mice_check_dataform"))



start <- Sys.time()
SACHS_100_MCAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=100, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_100_MAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=100, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_100_MNAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=100, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_1000_MCAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=1000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_1000_MAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=1000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_1000_MNAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=1000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_5000_MCAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=5000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_5000_MAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=5000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

start <- Sys.time()
SACHS_5000_MNAR <- parLapply(cl, 1:nrep, sim, network=sachs,
                                 trueMoral=trueMoral, n=5000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("SACHS", ls())],
     file="res_SACHS.RData")

stopCluster(cl)
