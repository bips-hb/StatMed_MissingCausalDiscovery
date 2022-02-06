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

source("gaussCItwd.R")
source("gaussMItest.R")
source("getSuff.R")


sim <- function(., scenario, n, q, m) {
  
  ## generate data
  Sigma <- matrix(0.2, 4, 4)
  Sigma[4,1] <- Sigma[1,4] <- 0.9
  Sigma[1,2] <- Sigma[2,1] <- 0.04 ## null hypothesis
  diag(Sigma) <- 1
  
  mu <- c(0,0,0,0)
  dat <- mvrnorm(n, mu, Sigma)
  dat <- data.frame(dat)
  colnames(dat) <- c("X", "Y", "Z", "A")
  Noise <- data.frame(matrix(rnorm(n*99), ncol=99))
  colnames(Noise) <- paste("N", 1:99, sep="")
  
  if (q==0) {
    p_twd <- p_MI <- gaussCItwd(1, 2, 3, suffStat=dat)
    
  } else {
    
    ## lose data
    miss <- sample(n, round(q*n)/100, replace=FALSE)
    if (scenario=="A") {
      dat[miss,3] <- NA
    } else {
      dat[miss,1] <- NA 
    }
    
    
    ####################
    # test-wise deletion
    ####################
  
    p_twd <- gaussCItwd(1, 2, 3, suffStat=dat)
    
    
    #####################
    # multiple imputation
    #####################
    
    if (scenario=="A" | scenario=="B") {
      dat2 <- dat[ ,c("X","Y","Z")]
      mi <- tryCatch( mice(dat2, m=m, method="norm", print=FALSE),
                      error=function(e){paste("Error")} )
    } else if (scenario=="C") {
      mi <- tryCatch( mice(dat, m=m, method="norm", print=FALSE),
                      error=function(e){paste("Error")} )
    } else if (scenario=="D") {
      dat3 <- cbind(dat, Noise)
      mi <- tryCatch( mice(dat3, m=m, method="norm", print=FALSE),
                      error=function(e){paste("Error")} )
    }
    
    if (identical(mi, "Error")) { return(c(NA,NA)) }
  
    p_MI <- tryCatch( gaussMItest(1, 2, 3,
                                  suffStat=getSuff(mi, test="gaussMItest")),
                      error=function(e){NA} )
    
  }
  
  res <- c(p_twd, p_MI)
  names(res) <- c("p_twd", "p_MI")
  
  return(res)
}


seed <- 30989
nrep <- 10000

cl<-makeCluster(119, type="MPI")
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(MASS)) # for data generation
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterExport(cl, c("gaussCItwd", "gaussMItest", "getSuff", "zStatMI",
                    "log.q1pm"))

### run the simulations:
res_A_50_00 <- parSapply(cl, 1:nrep, sim, scenario="A", n=50, q=0, m=100)
res_A_50_10 <- parSapply(cl, 1:nrep, sim, scenario="A", n=50, q=10, m=100)
res_A_50_30 <- parSapply(cl, 1:nrep, sim, scenario="A", n=50, q=30, m=100)
res_A_50_50 <- parSapply(cl, 1:nrep, sim, scenario="A", n=50, q=50, m=100)
res_A_50_70 <- parSapply(cl, 1:nrep, sim, scenario="A", n=50, q=70, m=100)

res_B_50_00 <- parSapply(cl, 1:nrep, sim, scenario="B", n=50, q=0, m=100)
res_B_50_10 <- parSapply(cl, 1:nrep, sim, scenario="B", n=50, q=10, m=100)
res_B_50_30 <- parSapply(cl, 1:nrep, sim, scenario="B", n=50, q=30, m=100)
res_B_50_50 <- parSapply(cl, 1:nrep, sim, scenario="B", n=50, q=50, m=100)
res_B_50_70 <- parSapply(cl, 1:nrep, sim, scenario="B", n=50, q=70, m=100)

res_C_50_00 <- parSapply(cl, 1:nrep, sim, scenario="C", n=50, q=0, m=100)
res_C_50_10 <- parSapply(cl, 1:nrep, sim, scenario="C", n=50, q=10, m=100)
res_C_50_30 <- parSapply(cl, 1:nrep, sim, scenario="C", n=50, q=30, m=100)
res_C_50_50 <- parSapply(cl, 1:nrep, sim, scenario="C", n=50, q=50, m=100)
res_C_50_70 <- parSapply(cl, 1:nrep, sim, scenario="C", n=50, q=70, m=100)

res_D_50_00 <- parSapply(cl, 1:nrep, sim, scenario="D", n=50, q=0, m=100)
res_D_50_10 <- parSapply(cl, 1:nrep, sim, scenario="D", n=50, q=10, m=100)
res_D_50_30 <- parSapply(cl, 1:nrep, sim, scenario="D", n=50, q=30, m=100)
res_D_50_50 <- parSapply(cl, 1:nrep, sim, scenario="D", n=50, q=50, m=100)
res_D_50_70 <- parSapply(cl, 1:nrep, sim, scenario="D", n=50, q=70, m=100)


res_A_500_00 <- parSapply(cl, 1:nrep, sim, scenario="A", n=500, q=0, m=100)
res_A_500_10 <- parSapply(cl, 1:nrep, sim, scenario="A", n=500, q=10, m=100)
res_A_500_30 <- parSapply(cl, 1:nrep, sim, scenario="A", n=500, q=30, m=100)
res_A_500_50 <- parSapply(cl, 1:nrep, sim, scenario="A", n=500, q=50, m=100)
res_A_500_70 <- parSapply(cl, 1:nrep, sim, scenario="A", n=500, q=70, m=100)

res_B_500_00 <- parSapply(cl, 1:nrep, sim, scenario="B", n=500, q=0, m=100)
res_B_500_10 <- parSapply(cl, 1:nrep, sim, scenario="B", n=500, q=10, m=100)
res_B_500_30 <- parSapply(cl, 1:nrep, sim, scenario="B", n=500, q=30, m=100)
res_B_500_50 <- parSapply(cl, 1:nrep, sim, scenario="B", n=500, q=50, m=100)
res_B_500_70 <- parSapply(cl, 1:nrep, sim, scenario="B", n=500, q=70, m=100)

res_C_500_00 <- parSapply(cl, 1:nrep, sim, scenario="C", n=500, q=0, m=100)
res_C_500_10 <- parSapply(cl, 1:nrep, sim, scenario="C", n=500, q=10, m=100)
res_C_500_30 <- parSapply(cl, 1:nrep, sim, scenario="C", n=500, q=30, m=100)
res_C_500_50 <- parSapply(cl, 1:nrep, sim, scenario="C", n=500, q=50, m=100)
res_C_500_70 <- parSapply(cl, 1:nrep, sim, scenario="C", n=500, q=70, m=100)

res_D_500_00 <- parSapply(cl, 1:nrep, sim, scenario="D", n=500, q=0, m=100)
res_D_500_10 <- parSapply(cl, 1:nrep, sim, scenario="D", n=500, q=10, m=100)
res_D_500_30 <- parSapply(cl, 1:nrep, sim, scenario="D", n=500, q=30, m=100)
res_D_500_50 <- parSapply(cl, 1:nrep, sim, scenario="D", n=500, q=50, m=100)
res_D_500_70 <- parSapply(cl, 1:nrep, sim, scenario="D", n=500, q=70, m=100)

stopCluster(cl)

save(list=ls()[grepl("res", ls())], file="res_CurveNull.RData")
