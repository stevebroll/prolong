# unexported gglasso functions

cv.ls <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(L2 = "Least-Squared loss", L1 = "Absolute loss")
    if (pred.loss == "default")
        pred.loss <- "L2"
    if (!match(pred.loss, c("L2", "L1"), FALSE)) {
        warning("Only 'L2' and 'L1'  available for least squares models; 'L2' used")
        pred.loss <- "L2"
    }
    predmat <- matrix(NA, length(y), length(lambda))
    nfolds <- max(foldid)
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- stats::predict(fitobj, x[which, , drop = FALSE])
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, L2 = (y - predmat)^2, L1 = abs(y - predmat))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}


getmin <- function(lambda, cvm, cvsd) {
    cvmin <- min(cvm)
    idmin <- cvm <= cvmin
    lambda.min <- max(lambda[idmin])
    idmin <- match(lambda.min, lambda)
    semin <- (cvm + cvsd)[idmin]
    idmin <- cvm <= semin
    lambda.1se <- max(lambda[idmin])
    list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

cv.gglasso_prolong <- function(x, y, group, lambda = NULL, pred.loss = c("L2"), nfolds = 5, foldid, delta) {
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
    gglasso.object <- gglasso::gglasso(x, y, group, lambda = lambda, delta = delta, intercept = F)
    lambda <- gglasso.object$lambda
    nfolds <- max(foldid)
    if (nfolds < 3) {
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    }
    outlist <- as.list(seq(nfolds))
    for (i in seq(nfolds)) {
        test <- which(foldid == i)
        y_sub <- y[-test]
        outlist[[i]] <- gglasso::gglasso(x = x[-test, , drop = FALSE], y = y_sub, group = group, lambda = lambda, delta = delta, intercept = F)
    }
    # fun <- paste('gglasso::cv', class(gglasso.object)[[2]], sep = '.')
    cvstuff <- do.call(cv.ls, list(outlist, lambda, x, y, foldid, pred.loss, delta))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gglasso"
    obj
}

# Unexported glmnet functions

cvtype <- function(type.measure = "mse", subclass = "elnet") {
    type.measures <- c("mse", "deviance", "class", "auc", "mae", "C")
    devname <- switch(subclass, elnet = "Mean-squared Error", lognet = "Binomial Deviance", fishnet = "Poisson Deviance", coxnet = "Partial Likelihood Deviance",
        multnet = "Multinomial Deviance", mrelnet = "Mean-squared Error", glmnetfit = "GLM Deviance")
    typenames <- c(deviance = devname, mse = "Mean-Squared Error", mae = "Mean Absolute Error", auc = "AUC", class = "Misclassification Error", C = "C-index")
    subclass.ch <- switch(subclass, elnet = c(1, 2, 5), lognet = c(2, 3, 4, 1, 5), fishnet = c(2, 1, 5), coxnet = c(2, 6), multnet = c(2, 3, 1, 5), mrelnet = c(1,
        2, 5), glmnetfit = c(2, 1, 5))
    subclass.type <- type.measures[subclass.ch]
    if (type.measure == "default") {
        type.measure <- subclass.type[1]
    }
    model.name <- switch(subclass, elnet = "Gaussian", lognet = "Binomial", fishnet = "Poisson", coxnet = "Cox", multnet = "Multinomial", mrelnet = "Multi-response Gaussian",
        glmnetfit = "GLM")
    if (!match(type.measure, subclass.type, FALSE)) {
        type.measure <- subclass.type[1]
        warning(paste("Only ", paste(subclass.type, collapse = ", "), " available as type.measure for ", model.name, " models; ", type.measure, " used instead",
            sep = ""), call. = FALSE)
    }
    names(type.measure) <- typenames[type.measure]
    type.measure
}

cvcompute <- function(cvstuff, foldid, nlams) {
    weights <- cvstuff$weights
    mat <- cvstuff$cvraw
    wisum <- tapply(weights, foldid, sum)
    nfolds <- max(foldid)
    outmat <- matrix(NA, nfolds, ncol(mat))
    good <- matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] <- NA
    for (i in seq(nfolds)) {
        mati <- mat[foldid == i, , drop = FALSE]
        wi <- weights[foldid == i]
        outmat[i, ] <- apply(mati, 2, stats::weighted.mean, w = wi, na.rm = TRUE)
        good[i, seq(nlams[i])] <- 1
    }
    N <- apply(good, 2, sum)
    list(cvraw = outmat, weights = wisum, N = N, type.measure = cvstuff$type.measure)
}


cvstats <- function(cvstuff, foldid, nfolds, lambda, nz, grouped, ...) {
    if (grouped) {
        nlams <- rep(dim(cvstuff$cvraw)[2], nfolds)
        cvstuff <- cvcompute(cvstuff, foldid, nlams)
    }
    cvm <- with(cvstuff, apply(cvraw, 2, stats::weighted.mean, w = weights, na.rm = TRUE))
    cvsd <- with(cvstuff, sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, stats::weighted.mean, w = weights, na.rm = TRUE)/(N - 1)))
    nas <- is.na(cvsd)
    if (any(nas)) {
        lambda <- lambda[!nas]
        cvm <- cvm[!nas]
        cvsd <- cvsd[!nas]
        nz <- nz[!nas]
    }
    list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd, cvlo = cvm - cvsd, nzero = nz)
}

cv.elnet <- function(predmat, y, type.measure, weights, foldid, grouped) {
    N = length(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure, mse = (y - predmat)^2, deviance = (y - predmat)^2, mae = abs(y - predmat))
    list(cvraw = cvraw, weights = weights, N = N, type.measure = type.measure, grouped = grouped)
}

getOptcv.glmnet <- function(lambda, cvm, cvsd, cvname) {
    if (match(cvname, c("AUC", "C-index"), 0))
        cvm = -cvm
    cvmin = min(cvm, na.rm = TRUE)
    idmin = cvm <= cvmin
    lambda.min = max(lambda[idmin], na.rm = TRUE)
    idmin = match(lambda.min, lambda)
    semin = (cvm + cvsd)[idmin]
    id1se = cvm <= semin
    lambda.1se = max(lambda[id1se], na.rm = TRUE)
    id1se = match(lambda.1se, lambda)
    index = matrix(c(idmin, id1se), 2, 1, dimnames = list(c("min", "1se"), "Lambda"))
    list(lambda.min = lambda.min, lambda.1se = lambda.1se, index = index)
}


cv.glmnet.raw <- function(x, y, weights, offset, lambda, type.measure, nfolds, foldid, alignment, grouped, keep, parallel, trace.it, glmnet.call, cv.call, ...) {
    if (trace.it) {
        cat("Training\n")
    }
    glmnet.object <- glmnet::glmnet(x, y, weights = weights, offset = offset, lambda = lambda, trace.it = trace.it, intercept = F, ...)
    glmnet.object$call <- glmnet.call
    subclass <- class(glmnet.object)[[1]]
    type.measure <- cvtype(type.measure, subclass)
    is.offset <- glmnet.object$offset
    if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
        nz <- stats::predict(glmnet.object, type = "nonzero")
        nz <- sapply(nz, function(x) {
            sapply(x, length)
        })
        nz <- ceiling(apply(nz, 1, stats::median))
    } else {
        nz <- sapply(stats::predict(glmnet.object, type = "nonzero"), length)
    }
    outlist <- as.list(seq(nfolds))
    N <- nrow(x)
    if (parallel) {
        `%dopar%` <- foreach::`%dopar%`
        outlist <- foreach::foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar% {
            fold <- which(foldid == i)
            if (length(dim(y)) > 1) {
                y_sub <- y[-fold, ]
            } else {
                y_sub <- y[-fold]
            }
            if (is.offset) {
                offset_sub <- as.matrix(offset)[-fold, ]
            } else {
                offset_sub <- NULL
            }
            glmnet::glmnet(x[-fold, , drop = FALSE], y_sub, lambda = lambda, offset = offset_sub, weights = weights[-fold], intercept = F, ...)
        }
    } else {
        for (i in seq(nfolds)) {
            if (trace.it) {
                cat(sprintf("Fold: %d/%d\n", i, nfolds))
            }
            fold <- which(foldid == i)
            if (length(dim(y)) > 1) {
                y_sub <- y[-fold, ]
            } else {
                y_sub <- y[-fold]
            }
            if (is.offset) {
                offset_sub <- as.matrix(offset)[-fold, ]
            } else {
                offset_sub <- NULL
            }
            outlist[[i]] <- glmnet::glmnet(x[-fold, , drop = FALSE], y_sub, lambda = lambda, offset = offset_sub, weights = weights[-fold], trace.it = trace.it,
                intercept = F, ...)
        }
    }
    lambda <- glmnet.object$lambda
    class(outlist) <- paste0(subclass, "list")
    predmat <- glmnet::buildPredmat(outlist, lambda, x, offset, foldid, alignment, y = y, weights = weights, grouped = grouped, type.measure = type.measure, family = stats::family(glmnet.object))
    fun <- paste("cv", subclass, sep = ".")
    cvstuff <- do.call(fun, list(predmat, y, type.measure, weights, foldid, grouped))
    grouped <- cvstuff$grouped
    if ((N/nfolds < 3) && grouped) {
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold", call. = FALSE)
        grouped <- FALSE
    }
    out <- cvstats(cvstuff, foldid, nfolds, lambda, nz, grouped)
    cvname <- names(cvstuff$type.measure)
    names(cvname) <- cvstuff$type.measure
    out <- c(out, list(call = cv.call, name = cvname, glmnet.fit = glmnet.object))
    if (keep) {
        out <- c(out, list(fit.preval = predmat, foldid = foldid))
    }
    lamin <- with(out, getOptcv.glmnet(lambda, cvm, cvsd, cvname))
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.glmnet"
    obj
}



cv.glmnet_prolong <- function(x, y, weights = NULL, offset = NULL, lambda = NULL, type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
    nfolds = NULL, foldid = NULL, alignment = c("lambda", "fraction"), keep = FALSE, parallel = FALSE, gamma = c(0, 0.25, 0.5, 0.75, 1), relax = FALSE, trace.it = 0,
    grouped = FALSE, ...) {
    type.measure <- match.arg(type.measure)
    alignment <- match.arg(alignment)
    if (!is.null(lambda) && length(lambda) < 2) {
        stop("Need more than one value of lambda for cv.glmnet")
    }
    if (!is.null(lambda) && alignment == "fraction") {
        warning("fraction of path alignment not available if lambda given as argument; switched to alignment=`lambda`")
        alignment <- "lambda"
    }
    N <- nrow(x)
    if (is.null(weights)) {
        weights <- rep(1, N)
    } else {
        weights <- as.double(weights)
    }
    y <- drop(y)
    cv.call <- glmnet.call <- match.call(expand.dots = TRUE)
    which <- match(c("type.measure", "nfolds", "foldid", "grouped", "keep"), names(glmnet.call), FALSE)
    if (any(which)) {
        glmnet.call <- glmnet.call[-which]
    }
    glmnet.call[[1]] <- as.name("glmnet")
    if (glmnet::glmnet.control()$itrace) {
        trace.it <- 1
    } else {
        if (trace.it) {
            glmnet::glmnet.control(itrace = 1)
            on.exit(glmnet::glmnet.control(itrace = 0))
        }
    }
    if (is.null(foldid)) {
        foldid <- sample(rep(seq(nfolds), length = N))
    } else {
        nfolds <- max(foldid)
    }
    if (nfolds < 3) {
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    }
    cv.glmnet.raw(x, y, weights, offset, lambda, type.measure, nfolds, foldid, alignment, grouped, keep, parallel, trace.it, glmnet.call, cv.call, ...)
}
