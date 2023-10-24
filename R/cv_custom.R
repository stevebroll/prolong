getmin <- function(lambda, cvm, cvsd) # same as unexported gglasso:::getmin
{
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

cv.gglasso_prolong <-
  function(x,
           y,
           group,
           lambda = NULL,
           pred.loss = c("misclass", "loss", "L1", "L2"),
           nfolds = 5,
           foldid,
           delta) {
    if (missing(pred.loss)) {
      pred.loss <- "default"
    } else {
      pred.loss <- match.arg(pred.loss)
    }
    N <- nrow(x)
    y <- drop(y)
    if (missing(delta)) {
      delta <- 1
    }
    if (delta < 0) {
      stop("delta must be non-negtive")
    }
    gglasso.object <- gglasso::gglasso(x, y, group,
      lambda = lambda,
      delta = delta, intercept = F
    )
    lambda <- gglasso.object$lambda
    nfolds <- max(foldid)
    if (nfolds < 3) {
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
    }
    outlist <- as.list(seq(nfolds))
    for (i in seq(nfolds)) {
      test <- which(foldid == i)
      y_sub <- y[-test]
      outlist[[i]] <- gglasso::gglasso(
        x = x[-test, , drop = FALSE],
        y = y_sub, group = group, lambda = lambda, delta = delta, intercept = F
      )
    }
    fun <- paste("gglasso::cv", class(gglasso.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(
      outlist, lambda, x, y, foldid,
      pred.loss, delta
    ))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm +
      cvsd, cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gglasso"
    obj
  }
