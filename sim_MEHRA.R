### small graph simulations
### MEHRA (mixed)

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

source("mixCItest.R")
source("mixCItwd.R")
source("mixMItest.R")
source("getSuff.R")
source("streetparty.R")
source("make_formulas_saturated.R")


# load network
mehra <- readRDS("mehra.rds")

# true adjacencies in moral graph
moralgraph <- moral(mehra)
trueMoral <- amat(moralgraph)
trueMoral <-
  trueMoral[c("Zone","Type","Year","Region","co","pm10","pm2.5","so2"),
            c("Zone","Type","Year","Region","co","pm10","pm2.5","so2")]


sim <- function(., network, trueMoral, n, mech) {
  ## network    a bn.fit object
  ## trueMoral  a matrix; true adjacencies in moral graph
  ## n          sample size
  ## mech       type of missingness mechanism: "MCAR", "MAR" or "MNAR"
  
  ### generate data ############################################################
  dat <- rbn(network, n=n)
  dat <- dat[ ,c("Zone","Type","Year","Region","co","pm10","pm2.5","so2")]
  
  ### analyse full data ########################################################
  datf <- dat
  datf[] <- lapply(datf, function(x){if (is.factor(x)) {factor(x)} else {x}})
  full_res <- pc(suffStat=datf, indepTest=mixCItest, alpha=0.05,
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
    # choose 4 variables containing missing values
    # desired pattern: for each observation, exactly one of these 4 is
    #   missing, and the probability depends on the other two(three)
    
    v1 <- sample(p, 4)
  
    # missinness patterns
    pat <- matrix(c(0,1,1,1,
                    1,0,1,1,
                    1,1,0,1,
                    1,1,1,0),
                  byrow=TRUE, nrow=4)
    
    amp <- ampute(data=dat[ ,v1], prop=0.999, patterns=pat,
                   freq=c(0.25,0.25,0.25,0.25))
    dat[ ,v1] <- amp$amp
    
    # in order to obtain a missingness proportion of 18% overall, choose one
    # further variable and delete values completely at random
    v2 <- sample(setdiff(1:p, v1), 1)
    dat[sample(n,round(n*p*0.055)),v2] <- NA
    
  } else if (mech=="MNAR") {
    # 8/5*18% of observations contain missing values in 5 selected variables
    # missingness depends on one of these 5, so MNAR
    
    del_index <- order(dat$pm10, decreasing=TRUE)[1:(round(n*0.18*8/5))]
    dat[del_index, c("pm10","Zone","Type","Year","pm2.5")] <- NA 
  }
  
  dat$Zone <- factor(dat$Zone)
  dat$Type <- factor(dat$Type)
  dat$Year <- factor(dat$Year)
  dat$Region <- factor(dat$Region)
  
  ### list-wise deletion #######################################################
  datl <- dat[complete.cases(dat), ]
  if (nrow(datl)<5) {
    lwd_res <- NA
  } else {
    
    lwd_res <- tryCatch( pc(suffStat=datl, indepTest=mixCItest, alpha=0.05,
                            labels=colnames(dat), skel.method="stable",
                            maj.rule=TRUE, solve.confl=TRUE),
                         error=function(e){paste("Error")} )
    
    if ( identical(lwd_res, "Error") ) {lwd_res <- NA}
  }
  
  ### test-wise deletion #######################################################
  twd_res <- tryCatch( pc(suffStat=dat, indepTest=mixCItwd, alpha=0.05,
                          labels=colnames(dat), skel.method="stable",
                          maj.rule=TRUE, solve.confl=TRUE),
                       error=function(e){paste("Error")} )
  
  if ( identical(twd_res, "Error") ) {twd_res <- NA}
  
  ### oracle multiple imputation ###############################################
  onelevel <- sapply(dat, nlevels)==1
  trueMoral[onelevel, ] <- 0
  trueMoral[ ,onelevel] <- 0
  form_omi <- make.formulas.saturated(dat, predictorMatrix=trueMoral, d=2)
  omi <- tryCatch( mice(dat, m=10, print=FALSE, formulas=form_omi,
                        defaultMethod=c("norm","logreg","polyreg","polr")),
                   error=function(e){paste("Error")} )
  if ( identical(omi, "Error") ){
    omi_res <- NA
  } else {
    omi_res <- pc(suffStat=getSuff(omi, test="mixMItest"), indepTest=mixMItest,
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
  mi <- tryCatch( mice(dat, m=10, print=FALSE, formulas=form_mi,
                       defaultMethod=c("norm","logreg","polyreg","polr")),
                  error=function(e){paste("Error")} )
  if ( identical(mi, "Error") ){
    mi_res <- NA
  } else {
    mi_res <- pc(suffStat=getSuff(mi, test="mixMItest"), indepTest=mixMItest,
                 alpha=0.05, labels=colnames(dat), skel.method="stable",
                 maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### standard multiple imputation without interaction terms ###################
  miw <- tryCatch( mice(dat, m=10, print=FALSE,
                        defaultMethod=c("norm","logreg","polyreg","polr")),
                   error=function(e){paste("Error")} )
  if ( identical(miw, "Error") ){
    miw_res <- NA
  } else {
    miw_res <- pc(suffStat=getSuff(miw, test="mixMItest"), indepTest=mixMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### multiple imputation using random forests #################################
  rf <- tryCatch( mice(dat, m=10, method="rf", print=FALSE, ntree=100),
                  error=function(e){paste("Error")} )
  if ( identical(rf, "Error") ){
    rf_res <- NA
  } else {
    rf_res <- pc(suffStat=getSuff(rf, test="mixMItest"), indepTest=mixMItest,
                 alpha=0.05, labels=colnames(dat), skel.method="stable",
                 maj.rule=TRUE, solve.confl=TRUE)
  }
  
  ### multiple imputation using random forests - CALIBER #######################
  rfc <- tryCatch( mice(dat, m=10, print=FALSE, ntree=100,
                        defaultMethod=c("rfcont","rfcat","rfcat","polr")),
                   error=function(e){paste("Error")} )
  if ( identical(rfc, "Error") ){
    rfc_res <- NA
  } else {
    rfc_res <- pc(suffStat=getSuff(rfc, test="mixMItest"), indepTest=mixMItest,
                  alpha=0.05, labels=colnames(dat), skel.method="stable",
                  maj.rule=TRUE, solve.confl=TRUE)
  }

  ### mean/mode imputation #####################################################
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
  
  mode_res <- pc(suffStat=dat2, indepTest=mixCItest, alpha=0.05,
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
clusterCall(cl, function() library(Rfast)) # mixMItest
clusterExport(cl, c("covm", "df_f", "df_h", "dfCG", "evalCell", "evalJoint",
                    "likelihoodCell", "likelihoodJoint", "listsum", "maxCell",
                    "maxJoint", "mixCItest", "mixCItwd", "mixMItest",
                    "multinomialLikelihood", "plus", "getSuff", "streetparty",
                    "make.formulas.saturated", "mice_check_dataform"))



start <- Sys.time()
MEHRA_100_MCAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=100, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_100_MAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=100, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_100_MNAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=100, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_1000_MCAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=1000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_1000_MAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=1000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_1000_MNAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=1000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_5000_MCAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=5000, mech="MCAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_5000_MAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=5000, mech="MAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

start <- Sys.time()
MEHRA_5000_MNAR <- parLapply(cl, 1:nrep, sim, network=mehra,
                                 trueMoral=trueMoral, n=5000, mech="MNAR")
end <- Sys.time()
end-start

save(list=ls()[grepl("MEHRA", ls())],
     file="res_MEHRA.RData")

stopCluster(cl)
