### small graph simulations - analysis
### MEHRA (mixed)

library(graph)
library(pcalg)
library(abind)
library(bnlearn)
library(xtable)

# load results
load("res_MEHRA.RData")

# load network
mehra <- readRDS("mehra.rds")
                   
# true CPDAG
cpdag_mehra <- cpdag(mehra)
nodes_select <- c("Zone","Type","Year","Region","co","pm10","pm2.5","so2")
nodes_delete <- setdiff(nodes(cpdag_mehra), nodes_select)
for (i in nodes_delete) {
  cpdag_mehra <- remove.node(cpdag_mehra, i)
}
NELcpdag_mehra <- as.graphNEL(cpdag_mehra)


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
    g_true <- NELcpdag_mehra
    rr <- sapply(r, analyse, g_true)
    return(rr)
  })
  out <- apply(do.call(abind, list(stats, along=3)), c(1,2), mean, na.rm=TRUE)
  colnames(out) <- c("complete", "List-wise deletion", "-- without correction",
                     "-- oracle 2-way parametric", "-- 2-way parametric",
                     "-- main effects parametric", "-- random forests (rf)",
                     "-- CALIBER", "Mode imputation")
  return(t(out))
  #return(stats)
}

sum_MEHRA_100_MCAR <- analyse_all(MEHRA_100_MCAR)
sum_MEHRA_100_MAR  <- analyse_all(MEHRA_100_MAR)
sum_MEHRA_100_MNAR <- analyse_all(MEHRA_100_MNAR)

sum_MEHRA_1000_MCAR <- analyse_all(MEHRA_1000_MCAR)
sum_MEHRA_1000_MAR  <- analyse_all(MEHRA_1000_MAR)
sum_MEHRA_1000_MNAR <- analyse_all(MEHRA_1000_MNAR)

sum_MEHRA_5000_MCAR <- analyse_all(MEHRA_5000_MCAR)
sum_MEHRA_5000_MAR  <- analyse_all(MEHRA_5000_MAR)
sum_MEHRA_5000_MNAR <- analyse_all(MEHRA_5000_MNAR)

xtable(sum_MEHRA_100_MCAR, digits=c(0,1,2,2,1,1))
xtable(sum_MEHRA_1000_MCAR, digits=c(0,1,2,2,1,1))
xtable(sum_MEHRA_5000_MCAR, digits=c(0,1,2,2,1,1))

xtable(sum_MEHRA_100_MAR, digits=c(0,1,2,2,1,1))
xtable(sum_MEHRA_1000_MAR, digits=c(0,1,2,2,1,1))
xtable(sum_MEHRA_5000_MAR, digits=c(0,1,2,2,1,1))

xtable(sum_MEHRA_100_MNAR, digits=c(0,1,2,2,1,1))
xtable(sum_MEHRA_1000_MNAR, digits=c(0,1,2,2,1,1))
xtable(sum_MEHRA_5000_MNAR, digits=c(0,1,2,2,1,1))


### Error-Meldungen zählen
error_all <- function(res) {
  is_error <- lapply(res, function(r) {
    is.na(r)
  })
  matr <- do.call(rbind, is_error)
  vec <- apply(matr, 2, mean)
  names(vec) <- c("complete","lwd","twd","omi","mi","miw","rf","rfc","mode")
  return(vec)
}

err_MEHRA_100_MCAR <- error_all(MEHRA_100_MCAR)
err_MEHRA_100_MAR  <- error_all(MEHRA_100_MAR)
err_MEHRA_100_MNAR <- error_all(MEHRA_100_MNAR)

err_MEHRA_1000_MCAR <- error_all(MEHRA_1000_MCAR)
err_MEHRA_1000_MAR  <- error_all(MEHRA_1000_MAR)
err_MEHRA_1000_MNAR <- error_all(MEHRA_1000_MNAR)

err_MEHRA_5000_MCAR <- error_all(MEHRA_5000_MCAR)
err_MEHRA_5000_MAR  <- error_all(MEHRA_5000_MAR)
err_MEHRA_5000_MNAR <- error_all(MEHRA_5000_MNAR)


tab <- rbind(err_MEHRA_100_MCAR, err_MEHRA_1000_MCAR,
             err_MEHRA_5000_MCAR, err_MEHRA_100_MAR,
             err_MEHRA_1000_MAR, err_MEHRA_5000_MAR,
             err_MEHRA_100_MNAR, err_MEHRA_1000_MNAR,
             err_MEHRA_5000_MNAR)
colnames(tab) <- c("complete", "List-wise del.", "Test-wise del.",
                     "Oracle 2-way", "2-way",
                   "Main effects",
                     "Random forests (rf)", "CALIBER", "Mean/mode")
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


save(sum_MEHRA_100_MCAR, sum_MEHRA_100_MAR, sum_MEHRA_100_MNAR,
    sum_MEHRA_1000_MCAR, sum_MEHRA_1000_MAR, sum_MEHRA_1000_MNAR,
    sum_MEHRA_5000_MCAR, sum_MEHRA_5000_MAR, sum_MEHRA_5000_MNAR,
     err_MEHRA_100_MCAR, err_MEHRA_100_MAR, err_MEHRA_100_MNAR,
    err_MEHRA_1000_MCAR, err_MEHRA_1000_MAR, err_MEHRA_1000_MNAR,
    err_MEHRA_5000_MCAR, err_MEHRA_5000_MAR, err_MEHRA_5000_MNAR, 
  file="summary_MEHRA.RData")
