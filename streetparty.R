### add neighbours of neighbours to prediction matrix

streetparty <- function(m) {
  ## m    adjacency matrix

  m2 <- m

  for (i in 1:nrow(m2)) {
    m2[i, ] <- m[i, ] + apply( m[ ,m[i, ]==1, drop=FALSE], 1, sum ) 
  }

  m3 <- 1*(m2 > 0)
  diag(m3) <- 0

  return(m3)
}
