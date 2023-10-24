get_delta_X <- function(X, n, p, t) {
  DXarray <- array(0, dim = c(n, p, t - 1))
  for (k in 1:(t - 1)) {
    DXarray[, , k] <-
      scale(t(X[, , k + 1]) - t(X[, , k]), center = F) / sqrt(n - 1)
  }

  DX <- matrix(0, nrow = n * (t - 1), ncol = sum(1:(t - 1)) * p)
  rowvec <- 1:n
  colvec <- 1:p
  for (i in 1:(t - 1)) {
    for (j in 1:i) {
      DX[rowvec, colvec] <- DXarray[, , j]
      colvec <- colvec + p
    }
    rowvec <- rowvec + n
  }

  return(list("DX" = DX, "DXarray" = DXarray))
}



get_delta_Y <- function(Y, n, t) {
  DY <- rep(0, n * (t - 1)) # DY will be n(t-1) length vector
  rowvec <- 1:n
  # loop through time points, assumes n_t same for all t
  for (i in 1:(t - 1)) {
    DY[rowvec] <- Y[, i + 1] - Y[, i] # First difference
    rowvec <- rowvec + n
  }
  return(DY)
}
