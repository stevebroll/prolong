

get_delta_X <- function(X, t, n, p) {
  DXarray <- array(0, dim = c(n, p, t - 1))
  for (i in 1:(t - 1)) {
    DXarray[, , i] <- scale(X[, , i + 1] - X[, , i], center = F)/sqrt(n-1)
  }

  DX = matrix(0, nrow = n*(t-1), ncol = sum(1:(t-1)) * p)
  rowvec = 1:n
  colvec = 1:p
  for(i in (t-1):1){
    for(j in i:1){
      DX[rowvec, colvec] = DXarray[,,j]
      colvec = colvec + p
    }
    rowvec = rowvec + n
  }
  return(DX)
}


get_delta_Y <- function(Y, t) {
  DY <- rep(0, n*(t-1))
  rowvec = 1:n
  for (i in 1:(t - 1)) {
    DY[rowvec] <- Y[, i + 1] - Y[, i]
    rowvec = rowvec + n
  }
  return(DY)
}



