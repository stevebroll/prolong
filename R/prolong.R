
prolong <- function(Y, X, wald.filter = F) {
  if (nrow(Y) != nrow(X)) stop("incompatible dimensions for X and Y")
  if (ncol(Y) != nrow(X)) stop("incompatible dimensions or differing # of time points for X and Y")

  n <- nrow(X)
  p <- ncol(X)
  t <- ncol(Y)

  if (wald.filter) print(0)
  DX <- get_delta_X(X, t)
  DY <- get_delta_Y(Y, t)
}
