### disCItest ohne Überprüfung, ob jede Variable mindestens zwei Werte hat
### außerdem wird immer getestet, egal, wie klein die Stichprobe ist

### bei twd wird sofort NA ausgegeben, wenn X oder Y nur einen Wert hat. Variablen in S mit nur einem Wert werden ebenfalls entfernt
disCItwd <- function(x, y, S=NULL, suffStat) {
  miss <- apply(suffStat$dm[, c(x, y, S)], 1, anyNA)
  
  if (sum(!miss) < 2) { return(1) }
  
  suffStat$dm <- suffStat$dm[!miss, c(x, y, S)]
  
  if (length(unique(suffStat$dm[ ,1]))==1) {return(NA)}
  if (length(unique(suffStat$dm[ ,2]))==1) {return(NA)}
  
  if (length(S) > 0) {
    S <- 3:ncol(suffStat$dm)
  }
  
  S1 <- apply(suffStat$dm[ ,S, drop=FALSE], 2, function(i){length(unique(i))==1})
  S <- S[!S1]
  if (length(S)==0) {S <- NULL}
  
  gSquareDis_new(x = 1, y = 2, S = S, dm=suffStat$dm, adaptDF=suffStat$adaptDF)
}

disCItest_new <- function(x, y, S=NULL, suffStat) {
    if (is.data.frame(dm <- suffStat$dm)) 
      dm <- data.matrix(dm)
    else stopifnot(is.matrix(dm))
    nlev <- suffStat$nlev
    adaptDF <- suffStat$adaptDF
    gSquareDis_new(x = x, y = y, S = S, dm = dm, nlev = nlev, adaptDF = adaptDF, 
        verbose = FALSE)
}


gSquareDis_new <- function (x, y, S, dm, nlev, adaptDF = FALSE, n.min = -1, 
          verbose = FALSE) {
  stopifnot((n <- nrow(dm)) >= 1, (p <- ncol(dm)) >= 2)
  if (!all(1 <= c(x, y, S) & c(x, y, S) <= p)) 
    stop("x, y, and S must all be in {1,..,p}, p=", 
         p)
  if (any(as.integer(dm) != dm)) 
    stop("'dm' must be discrete, with values in {0,1,..}")
  if (!any(dm == 0)) 
    stop("'dm' must have values in {0,1,..} with at least one '0' value")
  if (verbose) 
    cat("Edge ", x, "--", y, " with subset S =", 
        S, "\n")
  lenS <- length(S)
  if (missing(nlev) || is.null(nlev)) 
    nlev <- vapply(seq_len(p), function(j) length(levels(factor(dm[, 
                                                                   j]))), 1L)
  else stopifnot(is.numeric(nlev), length(nlev) == p, !is.na(nlev))
#  if (!all(nlev >= 2)) 
#    stop("Each variable, i.e., column of 'dm', must have at least two different values")
  nl.x <- nlev[x]
  nl.y <- nlev[y]
  nl.S <- nlev[S]
  df <- (nl.x - 1) * (nl.y - 1) * prod(nl.S)
  if (n < n.min) {
    warning(gettextf("n=%d is too small (n < n.min = %d ) for G^2 test (=> treated as independence)", 
                     n, n.min), domain = NA)
    return(1)
  }
  i.ny <- seq_len(nl.y)
  lenS <- length(S)
  d.x1 <- dm[, x] + 1L
  d.y1 <- dm[, y] + 1L
  if (lenS <= 4) {
    switch(lenS + 1L, {
      nij <- array(0L, c(nl.x, nl.y))
      for (i in 1:nl.x) {
        d.x.i <- d.x1 == i
        for (j in i.ny) nij[i, j] <- sum(d.x.i & d.y1 == 
                                           j)
      }
      t.X <- rowSums(nij)
      t.Y <- colSums(nij)
      t.log <- n * (nij/tcrossprod(t.X, t.Y))
      t.G2 <- 2 * nij * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }, {
      in.S <- seq_len(nl.S)
      dmS.1 <- dm[, S] + 1L
      nijk <- array(0L, c(nl.x, nl.y, nl.S))
      for (i in 1:nl.x) {
        d.x.i <- d.x1 == i
        for (j in i.ny) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for (k in in.S) nijk[i, j, k] <- sum(d.x.i.y.j & 
                                                 dmS.1 == k)
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(nl.x, nl.y, prod(nl.S)))
      for (k in 1:prod(nl.S)) t.log[, , k] <- nijk[, , 
                                                   k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }, {
      in.S1 <- seq_len(nl.S1 <- nl.S[1])
      in.S2 <- seq_len(nl.S2 <- nl.S[2])
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      nijk <- array(0L, c(nl.x, nl.y, nl.S1 * nl.S2))
      for (i in 1:nl.x) {
        d.x.i <- d.x1 == i
        for (j in i.ny) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for (k in in.S1) {
            d.x.y.S1 <- d.x.i.y.j & dmS1.1 == k
            for (l in in.S2) nijk[i, j, nl.S2 * (k - 
                                                   1) + l] <- sum(d.x.y.S1 & dmS2.1 == l)
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(nl.x, nl.y, prod(nl.S)))
      for (k in 1:prod(nl.S)) t.log[, , k] <- nijk[, , 
                                                   k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }, {
      in.S1 <- seq_len(nl.S1 <- nl.S[1])
      in.S2 <- seq_len(nl.S2 <- nl.S[2])
      in.S3 <- seq_len(nl.S3 <- nl.S[3])
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      nijk <- array(0L, c(nl.x, nl.y, prod(nl.S)))
      for (i1 in 1:nl.x) {
        d.x.i1 <- d.x1 == i1
        for (i2 in i.ny) {
          d.x.i1.y.i2 <- d.x.i1 & d.y1 == i2
          for (i3 in in.S1) {
            d.x.y.S1 <- d.x.i1.y.i2 & dmS1.1 == i3
            for (i4 in in.S2) {
              d.x.y.S1.2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in in.S3) nijk[i1, i2, nl.S3 * 
                                       nl.S2 * (i3 - 1) + nl.S3 * (i4 - 1) + 
                                       i5] <- sum(d.x.y.S1.2 & dmS3.1 == i5)
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(nl.x, nl.y, prod(nl.S)))
      for (k in 1:prod(nl.S)) t.log[, , k] <- nijk[, , 
                                                   k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
      t.G2 <- 2 * nijk * log(t.log)
      t.G2[which(is.nan(t.G2), arr.ind = TRUE)] <- 0
      G2 <- sum(t.G2)
    }, {
      in.S1 <- seq_len(nl.S1 <- nl.S[1])
      in.S2 <- seq_len(nl.S2 <- nl.S[2])
      in.S3 <- seq_len(nl.S3 <- nl.S[3])
      in.S4 <- seq_len(nl.S4 <- nl.S[4])
      dmS1.1 <- dm[, S[1]] + 1L
      dmS2.1 <- dm[, S[2]] + 1L
      dmS3.1 <- dm[, S[3]] + 1L
      dmS4.1 <- dm[, S[4]] + 1L
      nijk <- array(0L, c(nl.x, nl.y, prod(nl.S)))
      for (i1 in 1:nl.x) {
        d.x.i1 <- d.x1 == i1
        for (i2 in i.ny) {
          d.x.i1.y.i2 <- d.x.i1 & d.y1 == i2
          for (i3 in in.S1) {
            d.x.y.S1 <- d.x.i1.y.i2 & dmS1.1 == i3
            for (i4 in in.S2) {
              d.x.y.S1.2 <- d.x.y.S1 & dmS2.1 == i4
              for (i5 in in.S3) {
                d.x.y.S1.2.3 <- d.x.y.S1.2 & dmS3.1 == 
                  i5
                for (i6 in in.S4) nijk[i1, i2, nl.S4 * 
                                         nl.S3 * nl.S2 * (i3 - 1) + nl.S4 * 
                                         nl.S3 * (i4 - 1) + nl.S4 * (i5 - 1) + 
                                         i6] <- sum(d.x.y.S1.2.3 & dmS4.1 == 
                                                      i6)
              }
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)
      t.log <- array(0, c(nl.x, nl.y, prod(nl.S)))
      for (k in 1:prod(nl.S)) t.log[, , k] <- nijk[, , 
                                                   k] * (nk[k]/tcrossprod(nik[, k], njk[, k]))
      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    })
  }
  else {
    nijk <- array(0L, c(nl.x, nl.y, 1L))
    i <- d.x1[1]
    j <- d.y1[1]
    k <- NULL
    lapply(as.list(S), function(x) {
      k <<- cbind(k, d.x1)
      NULL
    })
    parents.count <- 1L
    parents.val <- t(k[1, ])
    nijk[i, j, parents.count] <- 1L
    for (it.sample in 2:n) {
      flag <- 0
      i <- d.x1[it.sample]
      j <- d.y1[it.sample]
      t.comp <- t(parents.val[1:parents.count, ]) == k[it.sample, 
                                                       ]
      dim(t.comp) <- c(lenS, parents.count)
      for (it.parents in 1:parents.count) {
        if (all(t.comp[, it.parents])) {
          nijk[i, j, it.parents] <- nijk[i, j, it.parents] + 
            1L
          flag <- 1
          break
        }
      }
      if (flag == 0) {
        parents.count <- parents.count + 1L
        if (verbose >= 2) 
          cat(sprintf(" adding new parents (count = %d) at sample %d\n", 
                      parents.count, it.sample))
        parents.val <- rbind(parents.val, k[it.sample, 
                                            ])
        nijk <- abind(nijk, array(0L, c(nl.x, nl.y, 1)))
        nijk[i, j, parents.count] <- 1L
      }
    }
    if (verbose && verbose < 2) 
      cat(sprintf(" added a total of %d new parents\n", 
                  parents.count))
    nik <- apply(nijk, 3, rowSums)
    njk <- apply(nijk, 3, colSums)
    nk <- colSums(njk)
    t.log <- array(0, c(nl.x, nl.y, parents.count))
    for (k in 1:parents.count) t.log[, , k] <- nijk[, , k] * 
      (nk[k]/tcrossprod(nik[, k], njk[, k]))
    t.G2 <- 2 * nijk * log(t.log)
    t.G2[which(is.nan(t.G2), arr.ind = TRUE)] <- 0
    G2 <- sum(t.G2)
  }
  if (adaptDF && lenS > 0) {
    zero.counts <- if (lenS == 0) 
      length(which(nij == 0))
    else length(which(nijk == 0)) + 4 * (2^lenS - dim(nijk)[3])
    df <- max((df - zero.counts), 1)
  }
  pchisq(G2, df, lower.tail = FALSE)
}




