

get_cor_matrix <- function(DX, n, p, t) {
  cordim <- sum(1:(t - 1)) * p
  cormat <- matrix(0, nrow = cordim, ncol = cordim)
  rowvec <- 1:n
  colstart <- 1
  colend <- 0
  for (i in (t - 1):1) {
    colend <- colend + i * p

    colvec <- seq(from = colstart, to = colend, by = 1)
    cormat[colvec, colvec] <- cor(DX[rowvec, colvec])

    rowvec <- rowvec + n
    colstart <- colend + 1
  }


  return(cormat)
}
