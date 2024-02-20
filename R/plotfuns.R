#' Interactive Heatmaps for Delta-Scale Pairwise Correlations
#'
#' @description
#' `delta_heatmap()` is designed to take an \eqn{n \times p \times t} array
#' and display a heatmap to visualize the observed delta-scalecorrelation
#' matrix. This delta-scale correlation matrix is used in the `prolong` model as
#' the adjacency matrix for the graph whose laplacian is used in the (group)
#' lasso + laplacian penalty used.
#'
#' @inheritParams prolong
#' @param timediff Pair of time points to use in heatmap. Should be in format `'t2-t1'`
#' @param object A `prolong` model object. Can be left `NULL` if `selected` is provided or if all variables are to be used
#' @param selected A character list of variable names or a numeric index of variables of interest to be included in the heatmap Can be left `NULL` if object is provided or if all variables are to be used
#' @param interactive If `TRUE`, a shiny app will open in browser that will display an interactive version of the heatmap where sub-heatmaps can be selected for display
#' @param grayscale If `TRUE` the viridis colors will be desaturated to grayscale
#'
#' @return A heatmap, possibly interactive in a shiny app
#'
#' @export
#'
#' @examples
#' \dontrun{
#' delta_heatmap(Xarray)
#' delta_heatmap(Xarray, timediff = '4-2', interactive = F, grayscale = T)
#' }
#'
#' @references
#' \insertRef{CH2}{prolong}
#' \insertRef{CH3}{prolong}
delta_heatmap <- function(x, timediff = "2-1", object = NULL, selected = NULL, interactive = TRUE,
    grayscale = FALSE) {
    t1 <- as.numeric(unlist(strsplit(timediff, "-"))[1])
    t2 <- as.numeric(unlist(strsplit(timediff, "-"))[2])
    x1 <- x[, , t2] - x[, , t1]

    if (!is.null(object) & !("prolong" %in% class(object))) {
        stop("Model object must be of class `prolong`")
    }
    if (!is.null(object)) {
        varnames <- object$selected
        varlist <- which(colnames(x) %in% varnames)
        x1 <- x1[, varlist]
    } else if (!is.null(selected)) {
        if (is.character(selected)) {
            varnames <- selected
            varlist <- which(colnames(x) %in% varnames)
            x1 <- x1[, varlist]
        } else if (is.numeric(selected)) {
            x1 <- x1[, varlist]
        }
    }

    rr1 <- abs(stats::cor(x1))
    cols <- grDevices::hcl.colors(100, palette = "viridis")
    if (grayscale) {
        cols <- colorspace::desaturate(cols)
    }
    ht <- complexheatmap.2(rr1, distfun = function(v) {
        stats::as.dist(1 - v)
    }, col = cols, trace = "none", symm = T, keysize = 0.5, offsetRow = 0, offsetCol = 0)
    if (interactive) {
        InteractiveComplexHeatmap::htShiny(ht, title = paste(timediff, "Heatmap"),
            description = "")
    } else {
        ht
    }
}

#' Automatically Plot Trajectories of Variables Selected by `prolong()`
#'
#' @inheritParams delta_scatter
#' @param object A `prolong` model object. Can be left `NULL` if `selected` is provided
#' @param selected A character list of variable names or a numeric index of variables of interest whose trajectories are to be plotted. Can be left `NULL` if object is provided
#' @param colors Either a single color to plot all trajectories or an n length vector with a color for each subject
#' @param timelabs Optional numeric or character vector of labels for the t time points. If left `NULL`, 1:t will be used
#'
#' @return A 2d sequence of trajectory plots from facet_wrap()
#' @export
#'
#' @examples
#' \dontrun{
#' promod <- prolong(Ymatrix, Xarray)
#' plot_trajectories(Xarray, promod)
#' plot_trajectories(Xarray, selected = promod$selected)
#' }
plot_trajectories <- function(x, object = NULL, selected = NULL, timelabs = NULL,
    colors = "#00BFC4") {
    if (is.null(object) & is.null(selected)) {
        stop("Either a `prolong` model object or list of variable names must be supplied")
    }
    if (!is.null(object) & !("prolong" %in% class(object))) {
        stop("Model object must be of class `prolong`")
    }
    if (!is.null(object)) {
        varnames <- object$selected
        varlist <- which(colnames(x) %in% varnames)
    } else if (!is.null(selected)) {
        if (is.character(selected)) {
            varnames <- selected
            varlist <- which(colnames(x) %in% varnames)
        } else if (is.numeric(selected)) {
            varlist <- selected
        } else {
            stop("selected must be a numeric or character vector")
        }
    }


    n <- nrow(x)
    t <- dim(x)[3]
    plotdat <- matrix(NA, n * t, length(varlist))
    for (i in 1:n) {
        plotdat[(i - 1) * t + seq(t), ] <- t(x[i, varlist, ])
    }
    colnames(plotdat) <- colnames(x)[varlist]
    if (is.null(timelabs)) {
        timelabs <- 1:t
    }
    time <- factor(rep(timelabs, n), levels = timelabs)
    id <- rep(1:n, each = t)
    if (is.null(colors)) {
        colors <- "#00BFC4"
    } else if (length(colors == 1)) {
        colors <- colors
    } else if (length(colors == n)) {
        colors <- rep(colors, each = t)
    } else {
        stop("colors should either a single color for all lines or an n length
         vector with a color for each subject")
    }
    meltdf <- cbind(reshape2::melt(plotdat)[, 3:2], time, id)
    colnames(meltdf)[2] <- "varname"
    value <- meltdf$value  # to avoid check notes
    p <- ggplot2::ggplot(meltdf, ggplot2::aes(x = time, y = value, group = id)) +
        ggplot2::geom_line(color = colors)
    p + ggplot2::facet_wrap("varname", scales = "free")
}




#' Experimental 2D and 3D Scatter Plots for Delta-Scale Pairwise Correlations
#'
#' @inheritParams prolong
#' @param timediff1 First pair of time points for the x-axis of the scatter plot. Should be in format `'t2-t1'`
#' @param timediff2 Second pair of time points for the y-axis of the scatter plot. Should be in format `'t2-t1'`
#' @param timediff3 Optional third pair of time points for the z-axis of the scatter plot. Should be in format `'t2-t1'`
#' @param fisherz If `TRUE`, the Fisher z-transformation (atanh) will be applied to the correlations
#' @param interactive If `TRUE` and timediff3 is `NULL`, an interactive `plotly` plot with hover text containing correlations and variable names will be generated. If `FALSE` and timediff3 is `NULL`, a regular `ggplot2` scatterplot will be shown. If timediff3 is provided, this parameter will only determine whether hovertext is shown
#' @param digits Number of digits for rounding in the hovertext
#'
#' @return Either a 2D or 3D scatterplot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' delta_scatter(Xarray)
#' delta_scatter(Xarray, timediff3 = '4-3')
#' delta_scatter(Xarray, timediff1 = '3-1', timediff2 = '5-3', timediff3 = '7-5')
#' }
#'
#' @references
#' \insertRef{ggplot2}{prolong}
#' \insertRef{plotly}{prolong}
delta_scatter <- function(x, timediff1 = "2-1", timediff2 = "3-2", timediff3 = NULL,
    fisherz = TRUE, interactive = TRUE, digits = 3) {
    p <- ncol(x)

    t1 <- as.numeric(unlist(strsplit(timediff1, "-"))[1])
    t2 <- as.numeric(unlist(strsplit(timediff1, "-"))[2])
    x1 <- x[, , t2] - x[, , t1]

    t1 <- as.numeric(unlist(strsplit(timediff2, "-"))[1])
    t2 <- as.numeric(unlist(strsplit(timediff2, "-"))[2])
    x2 <- x[, , t2] - x[, , t1]

    corr1 <- stats::cor(x1)
    corr2 <- stats::cor(x2)

    if (fisherz) {
        corr1 <- atanh(corr1)
        corr2 <- atanh(corr2)
    }

    vec1 <- corr1[upper.tri(corr1)]
    vec2 <- corr2[upper.tri(corr2)]

    labs <- expand.grid(colnames(x), colnames(x))[as.vector(upper.tri(matrix(0, p,
        p))), ]
    if (is.null(timediff3)) {
        # get long mat
        rmat <- data.frame(round(vec1, digits), round(vec2, digits), labs)
        if (fisherz) {
            colnames(rmat) <- c(paste(timediff1, "Fisher Z-Transformed Correlations"),
                paste(timediff2, "Fisher Z-Transformed Correlations"), "Variable 1",
                "Variable 2")
        } else {
            colnames(rmat) <- c(paste(timediff1, "Correlations"), paste(timediff2,
                "Correlations"), "Variable 1", "Variable 2")
        }
        name1 <- colnames(rmat)[1]
        name2 <- colnames(rmat)[2]
        p <- ggplot2::ggplot(data = rmat, ggplot2::aes(x = get(name1), y = get(name2),
            text = c(paste("Variable 1:", rmat[, 3], "\nVariable 2:", rmat[, 4])))) +
            ggplot2::geom_point(size = 0.5) + ggplot2::xlab(name1) + ggplot2::ylab(name2)
        if (interactive) {
            # p <- p %>% plotly::toWebGL()
            plotly::ggplotly(p) %>%
                plotly::layout(hoverlabel = list(align = "left")) %>%
                plotly::partial_bundle()
        } else {
            p
        }
    } else {
        t1 <- as.numeric(unlist(strsplit(timediff3, "-"))[1])
        t2 <- as.numeric(unlist(strsplit(timediff3, "-"))[2])
        x3 <- x[, , t2] - x[, , t1]

        corr3 <- stats::cor(x3)

        if (fisherz) {
            corr3 <- atanh(corr3)
        }

        vec3 <- corr3[upper.tri(corr3)]

        # get long mat
        rmat <- data.frame(round(vec1, digits), round(vec2, digits), round(vec3,
            digits), labs)
        if (fisherz) {
            colnames(rmat) <- c(paste(timediff1, "Fisher Z-Transformed Correlations"),
                paste(timediff2, "Fisher Z-Transformed Correlations"), paste(timediff3,
                  "Fisher Z-Transformed Correlations"), "Variable 1", "Variable 2")
        } else {
            c(paste(timediff1, "Correlations"), paste(timediff2, "Correlations"),
                paste(timediff3, "Correlations"), "Variable 1", "Variable 2")
        }
        name <- colnames(rmat)[1:3]
        p <- ggplot2::ggplot(data = rmat, ggplot2::aes(x = rlang::.data[[name[1]]],
            y = rlang::.data[[name[2]]], text = c(paste("Variable 1:", rmat[, 3],
                "\nVariable 2:", rmat[, 4])))) + ggplot2::geom_point(size = 0.5)
        if (interactive) {
            p <- plotly::plot_ly(rmat, x = ~get(name[1]), y = ~get(name[2]), z = ~get(name[3]),
                hovertemplate = paste(timediff1, " Corr: %{x}<br>", timediff2, " Corr: %{y}<br>",
                  timediff3, " Corr: %{z}<br>", "%{text}", "<extra></extra>", sep = ""),
                text = c(paste("Variable 1:", rmat[, 4], "\nVariable 2:", rmat[,
                  5])))
            p <- p %>%
                plotly::add_markers(size = 1)
            p <- p %>%
                plotly::partial_bundle()
            p
        } else {
            p <- plotly::plot_ly(rmat, x = ~get(name[1]), y = ~get(name[2]), z = ~get(name[3]))
            p <- p %>%
                plotly::add_markers(size = 1)
            p <- p %>%
                plotly::partial_bundle()
            p
        }
    }
}




#' Experimental Network Plots for Delta-Scale Pairwise (Partial) Correlations
#'
#' @inheritParams prolong
#' @inheritParams delta_heatmap
#' @param partial If `TRUE`, partial correlations using the `ppcor` package will be used. If `FALSE`, correlations will be used
#' @param corr_thresh Correlation value threshold for edge inclusion in the graph
#' @param interactive If `TRUE`, an interactive `visNetwork` plot will open in the Viewer panel. If `FALSE`, a static `igraph` plot will be produced
#' @param method A character string indicating which correlation coefficient is to be computed, with 'pearson' as default
#'
#' @return Either an interactive `visNetwork` or static `igraph` plot
#' @export
#'
#' @examples
#' \dontrun{
#' delta_network(Xarray)
#' delta_network(Xarray, timediff = '4-2', corr_thresh = .9, method = 'spearman')
#' }
delta_network <- function(x, timediff = "2-1", partial = TRUE, corr_thresh = 0.75,
    interactive = TRUE, method = c("pearson", "kendall", "spearman")) {
    t1 <- as.numeric(unlist(strsplit(timediff, "-"))[1])
    t2 <- as.numeric(unlist(strsplit(timediff, "-"))[2])
    x1 <- x[, , t2] - x[, , t1]
    if (partial) {
        suppressWarnings(rr1 <- abs(ppcor::pcor(x1, method = method)$estimate))
    } else {
        rr1 <- abs(stats::cor(x1, method = method))
    }
    colnames(rr1) <- rownames(rr1) <- colnames(x)
    diag(rr1) <- 0
    rr1[which(rr1 < corr_thresh)] <- 0
    g1 <- igraph::graph_from_adjacency_matrix(rr1, mode = "undirected", weighted = TRUE)
    igraph::V(g1)$name <- colnames(x)
    g1 <- igraph::delete.vertices(g1, which(igraph::degree(g1) == 0))
    if (interactive) {
        visNetwork::visIgraph(g1)
    } else {
        igraph::plot.igraph(g1)
    }
}
