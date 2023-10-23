prolong <-
  function(Y,
           X,
           lambda1 = NULL,
           lambda2 = NULL,
           lambdar = NULL,
           groups = NULL) {
    if (nrow(Y) != nrow(X))
      stop("Incompatible dimensions for X and Y, X and Y should have same # rows, 1 for each sample")
    if (ncol(Y) != dim(X)[3])
      stop(
        "Incompatible dimensions for X and Y, Y should have t columns the third component of dim(X) should also be t"
      )
    if (is.null(groups))
      message("No groups supplied, ordinary lasso will be used instead of group lasso")
    if (!is.null(groups) &
        groups != sort(groups))
      stop(
        "Groups must consist of consecutive columns, with group numbers counting from 1 to # groups"
      )
    n <- nrow(X)
    p <- ncol(X)
    t <- ncol(Y)
    g <- length(unique(groups))

    DXout <- get_delta_X(X, n, p, t)
    DY <- get_delta_Y(Y, n, t)
    cormat <- get_cor_matrix(DXout$DXarray)
    graph <-
      igraph::graph_from_adjacency_matrix(cormat,
                                          mode = 'undirected',
                                          weighted = T,
                                          diag = F)
    lap <-
      as.matrix(igraph::laplacian_matrix(graph, normalized = T))

    # Optimization for l2
    if (is.null(lambda2) | is.null(lambdar)) {
      message('lambda2 and/or lambdar missing, optimizing over both')
      ZTZ = crossprod(DXout$DX)
      ZTY = crossprod(DXout$DX, DY)
      YZT = crossprod(DY, DXout$DX)

      minfun =  function(l) {
        n * log(crossprod(Y, Y) - YTZ %*% solve(abs(l[1]) * (Q + diag(l[2], nrow(
          Q
        )))
        + ZTZ) %*% ZTY) + log(abs(det(abs(l[1]) *
                                        (
                                          Q + diag(l[2], nrow(Q))
                                        )
                                      + crossprod(Z, Z)))) - log(abs(det(abs(l[1]) *
                                                                           (
                                                                             Q + diag(l[2], nrow(Q))
                                                                           ))))
      }
      opt = optim(c(1, 0.1), minfun)
      print(opt)
      lambda2 = opt$par[1]
      lambdar = opt$par[2]
    }



  }
