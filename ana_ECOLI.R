### small graph simulations - analysis
### ECOLI (Gaussian)

library(graph)
library(pcalg)
library(abind)
library(bnlearn)
library(xtable)

# load results
load("res_ECOLI.RData")

# load network
ecoli <- readRDS("ecoli70.rds")
                   
# true CPDAG
cpdag_ecoli <- cpdag(ecoli)
nodes_select <- c("b1191", "cchB", "eutG", "fixC", "ibpB", "sucA", "tnaA",
                  "yceP", "yfaD", "ygbD", "ygcE", "yjbO")
nodes_delete <- setdiff(nodes(cpdag_ecoli), nodes_select)
for (i in nodes_delete) {
  cpdag_ecoli <- remove.node(cpdag_ecoli, i)
}
NELcpdag_ecoli <- as.graphNEL(cpdag_ecoli)



mycompareGraphs <- function (gl, gt) {
  ml <- wgtMatrix(ugraph(gl))
  mt <- wgtMatrix(ugraph(gt))
  p <- dim(ml)[2]
  mt[mt != 0] <- rep(1, sum(mt != 0))
  ml[ml != 0] <- rep(1, sum(ml != 0))
  diffm <- ml - mt
  nmbTrueGaps <- (sum(mt == 0) - p)/2
  fpr <- if (nmbTrueGaps == 0) 
    1
  else (sum(diffm > 0)/2)/nmbTrueGaps
  diffm2 <- mt - ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  tpr <- if (nmbTrueEdges == 0) 
    NA ##
  else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
  trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
  tdr <- if (sum(ml == 1) == 0) {
    if (trueEstEdges == 0) 
      1
    else 0
  }
  else trueEstEdges/(sum(ml == 1)/2)
  c(tpr = tpr, fpr = fpr, tdr = tdr)
}


analyse <- function(amat, g_true){
  if (anyNA(amat)) {return(rep(NA,5))}
  m <- (as(amat, "matrix")>0)*1
  g <- pcalg::getGraph(t(m))
  e <- sum(m) - sum(m*t(m)/2)
  shd <- pcalg::shd(as(g, "graphNEL"), g_true)
  u <- ugraph(pcalg::getGraph(g))
  h <- hamming(as.bn(u), as.bn(g_true))
  com <- mycompareGraphs(u, g_true)
  a <- c(e=e, rec=com[1], pre=com[3], h=h, shd=shd)
  return(a)
}


analyse_all <- function(res) {
  stats <- lapply(res, function(r){
    g_true <- NELcpdag_ecoli
    rr <- sapply(r, analyse, g_true)
    return(rr)
  })
  out <- apply(do.call(abind, list(stats, along=3)), c(1,2), mean, na.rm=TRUE)
  colnames(out) <- c("complete", "List-wise deletion", "-- without correction",
                     "-- permutation correction", "-- density correction",
                     "-- oracle linear models", "-- linear models",
                     "-- random forests (rf)", "-- CALIBER", "Mean imputation")
  return(t(out))
  #return(stats)
}

sum_ECOLI_100_MCAR <- analyse_all(ECOLI_100_MCAR)
sum_ECOLI_100_MAR  <- analyse_all(ECOLI_100_MAR)
sum_ECOLI_100_MNAR <- analyse_all(ECOLI_100_MNAR)

sum_ECOLI_1000_MCAR <- analyse_all(ECOLI_1000_MCAR)
sum_ECOLI_1000_MAR  <- analyse_all(ECOLI_1000_MAR)
sum_ECOLI_1000_MNAR <- analyse_all(ECOLI_1000_MNAR)

sum_ECOLI_5000_MCAR <- analyse_all(ECOLI_5000_MCAR)
sum_ECOLI_5000_MAR  <- analyse_all(ECOLI_5000_MAR)
sum_ECOLI_5000_MNAR <- analyse_all(ECOLI_5000_MNAR)

xtable(sum_ECOLI_100_MCAR, digits=c(0,1,2,2,1,1))
xtable(sum_ECOLI_1000_MCAR, digits=c(0,1,2,2,1,1))
xtable(sum_ECOLI_5000_MCAR, digits=c(0,1,2,2,1,1))

xtable(sum_ECOLI_100_MAR, digits=c(0,1,2,2,1,1))
xtable(sum_ECOLI_1000_MAR, digits=c(0,1,2,2,1,1))
xtable(sum_ECOLI_5000_MAR, digits=c(0,1,2,2,1,1))

xtable(sum_ECOLI_100_MNAR, digits=c(0,1,2,2,1,1))
xtable(sum_ECOLI_1000_MNAR, digits=c(0,1,2,2,1,1))
xtable(sum_ECOLI_5000_MNAR, digits=c(0,1,2,2,1,1))


### Error-Meldungen zählen
error_all <- function(res) {
  is_error <- lapply(res, function(r) {
    is.na(r)
  })
  matr <- do.call(rbind, is_error)
  vec <- apply(matr, 2, mean)
  names(vec) <- c("complete","lwd","twd","Tu1","Tu2","omi","mi","rf",
                  "rfc","mean")
  return(vec)
}

err_ECOLI_100_MCAR <- error_all(ECOLI_100_MCAR)
err_ECOLI_100_MAR  <- error_all(ECOLI_100_MAR)
err_ECOLI_100_MNAR <- error_all(ECOLI_100_MNAR)

err_ECOLI_1000_MCAR <- error_all(ECOLI_1000_MCAR)
err_ECOLI_1000_MAR  <- error_all(ECOLI_1000_MAR)
err_ECOLI_1000_MNAR <- error_all(ECOLI_1000_MNAR)

err_ECOLI_5000_MCAR <- error_all(ECOLI_5000_MCAR)
err_ECOLI_5000_MAR  <- error_all(ECOLI_5000_MAR)
err_ECOLI_5000_MNAR <- error_all(ECOLI_5000_MNAR)

tab <- rbind(err_ECOLI_100_MCAR, err_ECOLI_1000_MCAR, err_ECOLI_5000_MCAR,
      err_ECOLI_100_MAR, err_ECOLI_1000_MAR, err_ECOLI_5000_MAR,
      err_ECOLI_100_MNAR, err_ECOLI_1000_MNAR, err_ECOLI_5000_MNAR)
colnames(tab) <- c("complete", "List-wise del.", "Test-wise del.",
                     "Permutation corr.", "Density corr.",
                     "Oracle linear", "Linear",
                     "Random forests (rf)", "CALIBER", "Mean imputation")
rownames(tab) <- c("100_MCAR","100_MAR","100_MNAR","1000_MCAR","1000_MAR",
                   "1000_MNAR","5000_MCAR","5000_MAR","5000_MNAR")
del <- numeric(9)
for (i in 1:ncol(tab)) {
  if ( sum(tab[ ,i]) == 0 ) {del[i] <- 1}
}
tab <- tab[ ,!del]
tab <- t(tab) * 100
tab <- apply(tab, c(1,2), function(i){if(i==0) {return("--")} else {i}})
tab <- cbind(tab[ ,1:3], " ", tab[ ,4:6], " ", tab[ ,7:9])

xtable(tab)


save(sum_ECOLI_100_MCAR, sum_ECOLI_100_MAR, sum_ECOLI_100_MNAR,
    sum_ECOLI_1000_MCAR, sum_ECOLI_1000_MAR, sum_ECOLI_1000_MNAR,
    sum_ECOLI_5000_MCAR, sum_ECOLI_5000_MAR, sum_ECOLI_5000_MNAR,
     err_ECOLI_100_MCAR, err_ECOLI_100_MAR, err_ECOLI_100_MNAR,
    err_ECOLI_1000_MCAR, err_ECOLI_1000_MAR, err_ECOLI_1000_MNAR,
    err_ECOLI_5000_MCAR, err_ECOLI_5000_MAR, err_ECOLI_5000_MNAR, 
  file="summary_ECOLI.RData")
