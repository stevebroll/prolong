get_cor_matrix <- function(DXarray, n, p, t) {
  cormat <- matrix(0, nrow = sum(1:(t - 1)) * p, ncol = sum(1:(t - 1)) * p)
  colvec = 0
  for (i in 1:(t - 1)) {
    DXtemp <- NULL
    colvec = 1:(p * i) + max(colvec)
    for (j in 1:i) {
      DXtemp = cbind(DXtemp, DXarray[, , j])
    }
    cormat[colvec, colvec] = stats::cor(DXtemp)
  }
  return(abs(cormat))
}
