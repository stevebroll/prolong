get_delta_X <- function(X, t, n, p) {
  NULL

}



get_delta_Y <- function(Y, n, t) {
  DY <- rep(0, n * (t - 1)) # DY will be n(t-1) length vector
  rowvec <- 1:n
  # loop through time points, assumes n_t same for all t
  for (t in 1:(t - 1)) {
    DY[rowvec] <- Y[, i + 1] - Y[, i] # First difference
    rowvec <- rowvec + n
  }
  return(DY)
}
