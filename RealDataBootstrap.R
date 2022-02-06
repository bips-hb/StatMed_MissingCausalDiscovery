library(parallel)

dat <- readRDS("IDEFICS_data.RDS") # we are not allowed to share the data

boot <- function(., fulldat) {
  tiers <- rep(c(1,2,3,4,5,6,7,8), c(5,1,1,1,3,2,3,7))
  forbEdges <- matrix(0, ncol=23, nrow=23)
  forbEdges[1:4, 5] <- TRUE
  forbEdges[1:13, 12] <- TRUE
  
  ### draw bootstrap sample
  ind <- sample(1:657, 657, replace=TRUE)
  dat <- fulldat[ind, ]

  ### list-wise deletion
  datl <- dat[complete.cases(dat), ]
  res_lwd <- tpc(suffStat=datl, indepTest=mixCItest, labels=colnames(dat),
                 alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)

  ### test-wise deletion
  res_twd <- tpc(suffStat=dat, indepTest=mixCItwd, labels=colnames(dat),
                 alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
  
  ### parametric main effects multiple imputation
  miw <- tryCatch(
         mice(dat, m=100, defaultMethod=c("norm","logreg","polyreg","polr")),
         error=function(e){paste("Error")} )
  if ( identical(miw, "Error") ){
    res_miw <- NA
  } else {
    suffmiw <- complete(miw, action="all")
    res_miw <- tpc(suffStat=suffmiw, indepTest=mixMItest, labels=colnames(dat),
                   alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
  }

  ## random forests
  rf <- tryCatch( mice(dat, m=100, method="rf", ntree=100),
                  error=function(e){paste("Error")} )
  if ( identical(rf, "Error") ){
    res_rf <- NA
  } else {
    suffrf <- complete(rf, action="all")
    res_rf <- tpc(suffStat=suffrf, indepTest=mixMItest, labels=colnames(dat),
                  alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
  }

  ## CALIBER
  rfc <- tryCatch( mice(dat, m=100, ntree=100,
              defaultMethod=c("rfcont","rfcat","rfcat","polr")),
                   error=function(e){paste("Error")} )
  if ( identical(rfc, "Error") ){
    res_rfc <- NA
  } else {
    suffrfc <- complete(rfc, action="all")
    res_rfc <- tpc(suffStat=suffrfc, indepTest=mixMItest, labels=colnames(dat),
                   alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)
  }

  ## mean/mode imputation
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
  
  res_mode <- tpc(suffStat=dat2, indepTest=mixCItest, labels=colnames(dat),
                  alpha=0.1, tiers=tiers, forbEdges=forbEdges, verbose=TRUE)

  res <- list(res_lwd, res_twd, res_rf, res_rfc, res_miw, res_mode)
  res <- lapply(res, function(r){
    if(!inherits(r, "pcAlgo")){return(r)}else{return(as(r, "amat"))}
  })
  return(res)
}


seed <- 28359
nrep <- 10

cl<-makeCluster(10)
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(tpc))
clusterCall(cl, function() library(micd))
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterCall(cl, function() library(CALIBERrfimpute)) # multiple imputation

start <- Sys.time()
boot_res_1 <- parLapply(cl, 1:nrep, boot, fulldat=fulldat)
end <- Sys.time()
end-start

stopCluster(cl)

save(boot_res_1, file="boot_res_1.RData")


seed <- 30989
nrep <- 10

cl<-makeCluster(10)
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(tpc))
clusterCall(cl, function() library(micd))
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterCall(cl, function() library(CALIBERrfimpute)) # multiple imputation


start <- Sys.time()
boot_res_2 <- parLapply(cl, 1:nrep, boot, fulldat=fulldat)
end <- Sys.time()
end-start

stopCluster(cl)

save(boot_res_2, file="boot_res_2.RData")


seed <- 54321

cl<-makeCluster(6)
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(tpc))
clusterCall(cl, function() library(micd))
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterCall(cl, function() library(CALIBERrfimpute)) # multiple imputation


start <- Sys.time()
boot_res_3 <- parLapply(cl, 1:6, boot, fulldat=fulldat)
end <- Sys.time()
end-start

save(boot_res_3, file="boot_res_3.RData")

stopCluster(cl)


seed <- 12345

cl<-makeCluster(6)
if (!is.null(seed)) {clusterSetRNGStream(cl, seed)}
clusterCall(cl, function() library(tpc))
clusterCall(cl, function() library(micd))
clusterCall(cl, function() library(crayon)) # mice only loads when this one is loaded
clusterCall(cl, function() library(vctrs)) # mice only loads when this one is loaded
clusterCall(cl, function() library(backports)) # mice only loads when this one is loaded
clusterCall(cl, function() library(mice)) # multiple imputation
clusterCall(cl, function() library(CALIBERrfimpute)) # multiple imputation


start <- Sys.time()
boot_res_4 <- parLapply(cl, 1:6, boot, fulldat=fulldat)
end <- Sys.time()
end-start

save(boot_res_4, file="boot_res_4.RData")

start <- Sys.time()
boot_res_5 <- parLapply(cl, 1:6, boot, fulldat=fulldat)
end <- Sys.time()
end-start

save(boot_res_5, file="boot_res_5.RData")

start <- Sys.time()
boot_res_6 <- parLapply(cl, 1:6, boot, fulldat=fulldat)
end <- Sys.time()
end-start

save(boot_res_6, file="boot_res_6.RData")

start <- Sys.time()
boot_res_7 <- parLapply(cl, 1:6, boot, fulldat=fulldat)
end <- Sys.time()
end-start

save(boot_res_7, file="boot_res_7.RData")

stopCluster(cl)



### load results

load("boot_res_1.RData")
load("boot_res_2.RData")
load("boot_res_3.RData")
load("boot_res_4.RData")
load("boot_res_5.RData")
load("boot_res_6.RData")
load("boot_res_7.RData")

boot_res <- c(boot_res_1, boot_res_2, boot_res_3, boot_res_4, boot_res_5,
              boot_res_6, boot_res_7)

### count number of edges
ne <- sapply(boot_res, function(i) {
  sapply(i, function(j) {
    # matrix of undirected edges
    jj <- (j + t(j)) > 0
    # count undirected edges
    sum(jj) /2
  })
})

apply(ne, 1, min)
apply(ne, 1, mean)
apply(ne, 1, max)


### count number of edges adjoining critical nodes (mvpa, sed, sleep, homa)
nce <- sapply(boot_res, function(i) {
  sapply(i, function(j) {
    # matrix of undirected edges
    jj <- (j + t(j)) > 0
    # count edge marks
    total <- sum(jj[c("mvpa","sed","sleep","homa"), ])
    # sub-matrix of critical nodes
    jjj <- jj[c("mvpa","sed","sleep","homa"),c("mvpa","sed","sleep","homa")]
    # count edge marks in sub-matrix and subtract from total count
    within <- sum(jjj) / 2
    total - within
  })
})

apply(nce, 1, min)
apply(nce, 1, mean)
apply(nce, 1, max)
