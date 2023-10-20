prolong <- function(Y, X, lambda1 = NULL, lambda2 = NULL, groups = NULL) {
  if (nrow(Y) != nrow(X)) stop("Incompatible dimensions for X and Y, X and Y should have same # rows, 1 for each sample")
  if (ncol(Y) != dim(X)[3]) stop("Incompatible dimensions for X and Y, Y should have t columns the third component of dim(X) should also be t")
  if (is.null(groups)) message("No groups supplied, ordinary lasso will be used instead of group lasso")
  if(!is.null(groups) & groups != sort(groups)) stop("Groups must consist of consecutive columns, with group numbers counting from 1 to # groups")
  n <- nrow(X)
  p <- ncol(X)
  t <- ncol(Y)
  g <- length(unique(groups))

}
