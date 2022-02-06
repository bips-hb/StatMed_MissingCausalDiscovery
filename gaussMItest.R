gaussMItest <- function (x, y, S, suffStat) {
  # number of imputations
  M <- length(suffStat) - 1
  # sample size
  n <- suffStat[[M+1]]
  suffStat[[M+1]] <- NULL
  
  z <- sapply(suffStat, function(j) {
    zStatMI(x, y, S, C=j, n=n)
  })
  
  # 1. Average of M imputed data sets
  avgz <- mean(z)
  
  # 2. Average of completed-data variance
  W <- 1 / (n - length(S) - 3)
  
  # 3. Between variance
  B <- sum( ( z - avgz )^2 ) / (M-1)
  
  # 4. Total variance
  TV <- W + (1 + 1 / M) * B
  
  # 5. Test statistic
  ts <- avgz / sqrt(TV)
  
  # 6. Degrees of freedom
  df <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2
  
  # 7. pvalue
  pvalue <- 2 * pt(abs(ts), df = df, lower.tail = FALSE)
  
  return(pvalue)
}


zStatMI <- function (x, y, S, C, n) {
  r <- pcalg::pcorOrder(x, y, S, C)
  res <- 0.5 * log.q1pm(r)
  if (is.na(res))
    0
  else res
}

log.q1pm <- function(r) log1p(2*r/(1-r))
