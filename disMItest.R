disMItest <- function (x, y, S=NULL, suffStat) {
  
  # number of imputations / completed data sets
  M <- length(suffStat)
   
  nl_x <- nlevels(suffStat[[1]][ ,x])
  nl_y <- nlevels(suffStat[[1]][ ,y])
  
  if (nl_x==1) {return(1)}
  if (nl_y==1) {return(1)}
  
  if (length(S)==0) {
    nl_S <- 1
  } else{
    nl_S <- sapply(S, function(s){nlevels(suffStat[[1]][ ,s])})
  }
   
  # degrees of freedom if we had complete data
  k <- ( nl_x - 1 ) * ( nl_y - 1 ) * prod(nl_S)
  if (k==0) {k <- 1}
   
  if (length(S)==0) {
    mres <- lapply(1:M, function(m) {
       data_m <- suffStat[[m]][ ,c(x,y)]
       # observed counts
       o_list_m <- as.matrix(table(data_m[ ,1], data_m[ ,2]))
       # counts expected under marginal independence of x and y
       n <- sum(o_list_m)
       x_marg <- apply(o_list_m, 1, sum)
       y_marg <- apply(o_list_m, 2, sum)
       e_list_m <- x_marg %*% t(y_marg) / n
    
       o_m <- unlist(o_list_m)
       e_m <- unlist(e_list_m)
    
       # zero counts
       zeroes <- sum(o_m==0)
    
       # G^2 statistic for m-th completed data set
       # expected counts of zero are ignored, this is consistent with the pcalg procedure
       G2_m <- 2*sum( o_m * log(o_m/e_m) , na.rm=TRUE)
    
       return(list(o_m, e_m, G2_m, zeroes))
      })
      
   } else {
      mres <- lapply(1:M, function(m) {
         data_m <- suffStat[[m]][ ,c(x,y,S)]
         split_m <- split(data_m, data_m[ ,-c(1,2)])
         # observed counts
         o_list_m <- lapply(split_m, function(i) { as.matrix(table(i[ ,1], i[ ,2])) } )
         # counts expected under conditional independence of x and y given S
         e_list_m <- lapply(o_list_m, function(i) {
            n_i <- sum(i)
            x_marg <- apply(i, 1, sum)
            y_marg <- apply(i, 2, sum)
            xy_prod <- x_marg %*% t(y_marg) / n_i
         } )
      
         o_m <- unlist(o_list_m)
         e_m <- unlist(e_list_m)
      
         # zero counts
         zeroes <- sum(o_m==0)
      
         # G^2 statistic for m-th completed data set
         # expected counts of zero are ignored, this is consistent with the pcalg procedure
         G2_m <- 2*sum( o_m * log(o_m/e_m) , na.rm=TRUE)
      
         return(list(o_m, e_m, G2_m, zeroes))
      } )
   }
   
   # warning about zero counts
   #av_zero <- mean(sapply(mres, function(i){return(i[[4]])}))
   #ncells <- nl_x * nl_x * prod(nl_S)
   #if (verbose) {
    #  cat(paste("On average,", av_zero, "out of", ncells, "cells are empty in the imputed data. Empty cells bias the test result towards concluding independence."))
   #}
   
   # average G^2 statistic
   av_G2 <- mean(sapply(mres, function(i){return(i[[3]])}))
   
   # average observed cell counts
   av_o <- apply(sapply(mres, function(i){return(i[[1]])}), 1, mean)
   
   # average expected cell counts
   av_e <- apply(sapply(mres, function(i){return(i[[2]])}), 1, mean)
   
   # G^2 statistic with average observed and average expected counts for each imputed data set
   G2_av_m <- sapply(mres, function(m) {
      o_m <- m[[1]]
      2*sum( o_m * log(av_o/av_e) , na.rm=TRUE)
      })
   
   G2_av <- mean(G2_av_m)
   
   # r
   r <- (M+1)/(M-1)/k * (av_G2 - G2_av)
   
   # pooled G^2 statistic
   G2_pooled <- G2_av / (k*(r+1))
   
   # df for pooled G2
   t <- k*(M-1)
   if (t > 4) {
      df <- 4 + (t-4) * (1 + (1-2/t) / r)^2
   } else {
      df <- t * (1 + 1/k) * (1 + 1/r)^2 /2
   }
   
   # p-Wert
   pvalue <- pf(G2_pooled, k, df, lower.tail=FALSE)
   
   #return(c(G2=G2_pooled, df=k, p=pvalue))
   return(pvalue)
}

