### Run on the cluster:

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

library(snow) # for parallel computing
library(Rmpi) # for parallel computing
library(parallel) # for parallel computing



source("gaussMItest.R")
source("gaussCItwd.R")


sim <- function(., p, q,
            density=c("very_sparse", "sparse","medium", "dense", "very_dense"),
            weights=c("very_weak", "weak", "medium", "strong", "very_strong")) {
  
  prob <- switch(density, "very_sparse"=0.1, "sparse"=0.25, "medium"=0.4,
                 "dense"=0.55, "very_dense"=0.7)
  B <- switch(weights, "very_weak"=0.03, "weak"=0.065, "medium"=0.088,
              "strong"=0.111, "very_strong"=0.145)
  
  ## generate random DAG with specified edge weights
  DAG <- randomDAG(n=p, prob=prob, lB=B, uB=B, V=letters[1:p])
  
  ## generate 500 observations
  n <- 500
  dat <- rmvDAG(n=n, dag=DAG, errDist="normal")
  
  ## estimate graph using PC on the complete data
  pc_full <- pc(suffStat=list(C=cor(dat), n=n), labels=letters[1:p],
                indepTest=gaussCItest, alpha=0.05)
  
  ## delete data (MCAR)
  dat[sample(n*p, n*q*p, replace=FALSE)] <- NA
  dat <- as.data.frame(dat)
  
  ## estimate graph using PC with test-wise deletion on the incomplete data 
  pc_twd <- pc(suffStat=dat, labels=letters[1:p], indepTest=gaussCItwd,
               alpha=0.05)
  
  ## impute and estimate graph using PC on the multiply imputed data
  mi <- mice(dat, method="norm", m=100, print=FALSE)
  mi_suff <- c(lapply(complete(mi, action="all"), cor), n)
  pc_mi <- pc(suffStat=mi_suff, labels=letters[1:p], indepTest=gaussMItest,
              alpha=0.05)
  
  ## construct undirected graphs
  u_full <- ugraph(getGraph(pc_full))
  u_twd <- ugraph(getGraph(pc_twd))
  u_mi <- ugraph(getGraph(pc_mi))
  
  
  ## calculate Hamming distances
  H_full <- hamming(as.bn(u_full), as.bn(DAG))
  H_twd <- hamming(as.bn(u_twd), as.bn(DAG))
  H_mi <- hamming(as.bn(u_mi), as.bn(DAG))
  
  res <- c(H_full, H_twd, H_mi)
  return(res)
}


## run simulations

nrep <- 1000
seed <- 30659

cl<-makeCluster(239, type="MPI")
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(graph))
clusterCall(cl, function() library(pcalg))
clusterCall(cl, function() library(bnlearn)) # Hamming distance
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterExport(cl, c("gaussCItwd", "gaussMItest", "zStatMI","log.q1pm"))

res_vsparse_vweak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="very_sparse", weights="very_weak")
res_vsparse_weak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                              density="very_sparse", weights="weak")
res_vsparse_medium <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                                density="very_sparse", weights="medium")
res_vsparse_strong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                                density="very_sparse", weights="strong")
res_vsparse_vstrong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                                 density="very_sparse", weights="very_strong")

res_sparse_vweak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                              density="sparse", weights="very_weak")
res_sparse_weak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                             density="sparse", weights="weak")
res_sparse_medium <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="sparse", weights="medium")
res_sparse_strong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="sparse", weights="strong")
res_sparse_vstrong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                                density="sparse", weights="very_strong")

res_medium_vweak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                              density="medium", weights="very_weak")
res_medium_weak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                             density="medium", weights="weak")
res_medium_medium <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="medium", weights="medium")
res_medium_strong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="medium", weights="strong")
res_medium_vstrong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                                density="medium", weights="very_strong")

res_dense_vweak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                             density="dense", weights="very_weak")
res_dense_weak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                            density="dense", weights="weak")
res_dense_medium <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                              density="dense", weights="medium")
res_dense_strong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                              density="dense", weights="strong")
res_dense_vstrong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="dense", weights="very_strong")

res_vdense_vweak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                              density="very_dense", weights="very_weak")
res_vdense_weak <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                             density="very_dense", weights="weak")
res_vdense_medium <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="very_dense", weights="medium")
res_vdense_strong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                               density="very_dense", weights="strong")
res_vdense_vstrong <- parSapply(cl, 1:nrep, sim, p=8, q=0.1,
                                density="very_dense", weights="very_strong")

stopCluster(cl)

save(list=ls()[grepl("res", ls())], file="res_TilePlot.RData")
