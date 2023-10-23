prolong <-
  function(Y,
           X,
           lambda1 = NULL,
           lambda2 = NULL,
           lambdar = NULL,
           groups = TRUE) {
    if (nrow(Y) != nrow(X))
      stop("Incompatible dimensions for X and Y, X and Y should have same # rows, 1 for each sample")
    if (ncol(Y) != dim(X)[3])
      stop(
        "Incompatible dimensions for X and Y, Y should have t columns the third component of dim(X) should also be t"
      )
    if (is.null(groups) | isFALSE(groups))
      message("No groups supplied r suggested, ordinary lasso will be used instead of group lasso")
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
      opt <- optim(c(1, 1), minfun)
      lambda2 = opt$par[1]
      lambdar = opt$par[2]
    }

    # get incidence matrix
    LDL <- fastmatrix::ldl(lap + diag(lambdar, nrow = nrow(Q)))
    incidence <- LDL$lower %*% diag(sqrt(abs(LDL$d)))

    tri <- t*(t-1)/2
    Xaug <- (1/sqrt(1 + lambda2))*rbind(Xcomb, sqrt(lambda2) *t(incidence))
    Xaug <- Xaug[,rep((1:p),each = tri) + rep(seq(0,p*(tri-1),p),p)]
    Yaug <- c(Ycomb, rep(0, nrow(Xaug) - (t-1)*n))

    foldids <- c(rep(caret::createFolds(1:n,5, list = F),(t-1)))

    if(is.null(groups) | isFALSE(groups)){
      NULL

    } else {
      if(isTRUE(groups)){
        groups = rep(1:p, each = tri)
      }
      cv = gglasso::cv.gglasso(Xaug, Yaug, intercept = F, group = groups, foldid = foldids)
      gllmod = gglasso(Xaug, Yaug, intercept = F, group = groups, lambda = cv$lambda.1se)
      coefs = coef(gllmod)[-1,]
      coefs = coefs/(sqrt(1 + lambda2))
      names(coefs) = rep(colnames(DXout$DXarray), each = tri)


    }




    return(coefs)

  }
