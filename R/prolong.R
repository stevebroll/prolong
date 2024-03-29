#' Fit prolong Model
#'
#' @description
#'
#' `prolong()` is designed to take an \eqn{n \times t} outcome matrix and an
#' \eqn{n \times p \times t} array of covariates and automatically run the
#' entire processing and model fitting process, returning a list of coefficients
#' with corresponding variable names.
#'
#' @details
#' First, the data is reshaped into the first-differenced vectorized
#' \eqn{n(t-1)} length Y and first-differenced matricized \eqn{n(t-1) \times
#' pt(t-1)/2} X with a corresponding dependence matrix.
#' Next, hyperparameters for the graph laplacian-based network
#' penalty are found via MLE and the parameter for lasso/group lasso is found
#' via a careful implementation of cross-validation.
#' Lastly, a group lasso + laplacian or lasso + laplacian model is implemented,
#' and its bias-adjusted coefficients are returned.
#'
#' @param x Input covariate array, with n rows, p columns, and t slices
#' @param y Input response matrix, with n rows and t columns
#' @param lambda1 Lasso/group lasso parameter, if left as `NULL` this parameter
#' will be chosen via a cross-validation that keeps each subject's observations
#' across time points together. It is recommended to save the lambda2 and
#' lambdar values from the first run of prolong for future runs, since the
#' optimization step takes the longest
#' @param lambda2 Laplacian penalty parameter, if left as `NULL` will be chosen along
#'  with lambdar via MLE, see supplementary material in
#'  \insertCite{spatial-connectivity}{prolong}
#' @param lambdar Nuisance parameter, if left as `NULL` will be chosen along with
#'  lambdar via MLE. Added to the the diagonal elements of the laplacian matrix
#'   to get the invertibility required for the MLE
#' @param groups Optional pre-specified groups. If `NULL` or `FALSE`, lasso will be
#'  used. If left as `TRUE`, then there will be p groups each containing the
#'   observations across time points
#' @param foldids Optional pre-specified foldids for the cv. Should be of length
#'  n. If left as `NULL`, subjects will be automatically split into 5 folds
#' @param optimvals Initial values of lambda2 and lambdar to be used by
#' `optim()` for the MLE of lambda2 and lambdar
#'
#'
#'
#' @return A list object with S3 class `prolong`
#' \tabular{ll}{
#' \code{beta} \tab `p*t*t-1/2 *length(lambda1)` matrix of coefficients. Currently, `lambda1` is chosen via cross-validation and is just a single value \cr
#' \code{selected} \tab the names of the variables with at least one non-zero coefficient \cr
#' \code{df} \tab number of non-zero coefficients (out of `p*t(t-1)/2`, not `p`)\cr
#' \code{dim} \tab dimension of full coefficient matrix over `lambda1` values. Currently, `lambda1` is chosen via cross-validation and is just a single value \cr
#' \code{lambda1} \tab sequence of `lambda1` values used in final `gglasso`/`glmnet` call \cr
#' \code{lambda2} \tab `lambda2` value either passed by user or chosen via MLE, parameter of interest for the network penalty \cr
#' \code{lambdar} \tab `lambdar` value either passed by user or chosen via MLE, nuisance parameter needed to estimate `lambda2` via MLE \cr
#' \code{npasses} \tab total number of iterations summed over `lambda1` values for final `gglasso`/`glmnet` call \cr
#' \code{jerr} \tab error flag for final `gglasso`/`glmnet` call, 0 if no error \cr
#' \code{group} \tab vector of consecutive integers describing the grouping of coefficients \cr
#' \code{call} \tab the `gglasso`/`glmnet` call that produced this object \cr
#' }
#'
#' @export
#'
#' @importFrom Rdpack reprompt
#' @importFrom foreach %dopar%
#'
#' @examples
#' \dontrun{
#' promod <- prolong(Ymatrix, Xarray)
#' promod$beta
#' promod$selected
#'
#' promod <- prolong(Ymatrix, Xarray, lambda2 = .001, lambdar = 10, groups = FALSE)
#' promod$beta
#' promod$selected
#' }
#'
#' @references
#' \insertRef{prolong}{prolong}
#' \insertAllCited{}
#' \insertRef{gglasso}{prolong}
#' \insertRef{group-lasso}{prolong}
#' \insertRef{glmnet}{prolong}
#' \insertRef{lasso}{prolong}
#'
prolong <- function(x, y, lambda1 = NULL, lambda2 = NULL, lambdar = NULL, groups = TRUE,
    foldids = NULL, optimvals = c(1, 0.01)) {
    if (nrow(y) != nrow(x)) {
        stop("Incompatible dimensions for X and Y, X and Y should have same #
           rows, 1 for each sample")
    }
    if (ncol(y) != dim(x)[3]) {
        stop("Incompatible dimensions for X and Y, Y should have t columns the third
        component of dim(X) should also be t")
    }
    if (is.null(groups) | isFALSE(groups)) {
        message("No groups supplied or suggested, ordinary lasso will be used instead of group lasso")
    }
    if (!is.null(groups) & !identical(groups, sort(groups))) {
        stop("Groups must consist of consecutive columns, with group numbers counting from 1 to # groups")
    }
    n <- nrow(x)
    p <- ncol(x)
    t <- ncol(y)

    DXout <- get_delta_x(x, n, p, t)
    DY <- get_delta_y(y, n, t)
    cormat <- get_cor_matrix(DXout$DXarray, n, p, t)
    graph <- igraph::graph_from_adjacency_matrix(cormat, mode = "undirected", weighted = T,
        diag = F)
    lap <- as.matrix(igraph::laplacian_matrix(graph, normalized = T))

    # Optimization for l2
    if (is.null(lambda2) | is.null(lambdar)) {
        message("lambda2 and/or lambdar missing, optimizing over both")
        ZTZ <- crossprod(DXout$DX)
        ZTY <- crossprod(DXout$DX, DY)
        YTZ <- crossprod(DY, DXout$DX)
        YTY <- crossprod(DY)
        dn <- nrow(DXout$DX)
        opt <- stats::optim(par = optimvals, fn = minfun, lap = lap, dn = dn, YTY = YTY,
            YTZ = YTZ, ZTZ = ZTZ, ZTY = ZTY)
        lambda2 <- abs(opt$par[1])
        lambdar <- abs(opt$par[2])
        cat(paste("lambda2 = ", lambda2, "\nlambdar = ", lambdar, sep = ""))
    }

    # get incidence matrix
    LDL <- fastmatrix::ldl(lap + diag(lambdar, nrow = nrow(lap)))
    incidence <- LDL$lower %*% diag(sqrt(abs(LDL$d)))

    tri <- t * (t - 1)/2
    Xaug <- (1/sqrt(1 + lambda2)) * rbind(DXout$DX, sqrt(lambda2) * t(incidence))
    colnames(Xaug) <- rep(colnames(DXout$DXarray), tri)
    Xaug <- Xaug[, rep((1:p), each = tri) + rep(seq(0, p * (tri - 1), p), p)]
    Yaug <- c(DY, rep(0, nrow(Xaug) - (t - 1) * n))

    if (!is.null(foldids)) {
        foldids <- rep(foldids, (t - 1))
    } else {
        foldids <- c(rep(caret::createFolds(1:n, 5, list = F), (t - 1)))
    }

    if (is.null(groups) | isFALSE(groups)) {
        if (is.null(lambda1)) {
            cv <- cv.glmnet_prolong(Xaug, Yaug, foldid = foldids)
            lambda1 <- cv$lambda.1se
        } else {
            lambda1 <- lambda1
        }
        llmod <- glmnet::glmnet(Xaug, Yaug, intercept = F, lambda = lambda1)
        coefs <- stats::coef(llmod)[-1, ]
        coefs <- coefs/(sqrt(1 + lambda2))
        npasses <- llmod$npasses
        jerr <- llmod$jerr
        call <- llmod$call
    } else {
        if (isTRUE(groups)) {
            groups <- rep(1:p, each = tri)
        }
        if (is.null(lambda1)) {
            cv <- cv.gglasso_prolong(Xaug, Yaug, group = groups, foldid = foldids)
            lambda1 <- cv$lambda.1se
        } else {
            lambda1 <- lambda1
        }
        gllmod <- gglasso::gglasso(Xaug, Yaug, intercept = F, group = groups, lambda = lambda1)
        coefs <- stats::coef(gllmod)[-1, ]
        coefs <- coefs/(sqrt(1 + lambda2))
        npasses <- gllmod$npasses
        jerr <- gllmod$jerr
        call <- gllmod$call
    }

    names(coefs) <- colnames(Xaug)
    df <- length(which(coefs != 0))
    selected <- unique(names(coefs)[which(coefs != 0)])
    output <- list(beta = as.matrix(coefs), selected = selected, df = df, dim = c(p,
        length(lambda1)), lambda1 = lambda1, lambda2 = lambda2, lambdar = lambdar,
        npasses = npasses, jerr = jerr, group = groups, call = call)
    class(output) <- c("prolong", class(output))

    return(output)
}
