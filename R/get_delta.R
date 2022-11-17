

get_delta_X <- function(X, t) {
  DX <- array(NULL, dim = c(nrow(X), ncol(X), t - 1))
  for (i in 1:(t - 1)) {
    DX[, , i] <- X[, , i + 1] - X[, , i]
  }
  return(DX)
}


get_delta_Y <- function(Y, t) {
  DY <- matrix(NULL, nrow = nrow(Y), ncol = t - 1)
  for (i in 1:(t - 1)) {
    DY[, i] <- Y[, i + 1] - Y[, i]
  }
  return(DY)
}
