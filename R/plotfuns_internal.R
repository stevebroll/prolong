warning_wrap <- function(...) {
    # From ComplexHeatmap
    x <- paste0(...)
    x <- paste(strwrap(x), collapse = "\n")
    warning(x, call. = FALSE)
}

stop_wrap <- function(...) {
    # From ComplexHeatmap
    x <- paste0(...)
    x <- paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

message_wrap <- function(...) {
    # from ComplexHeatmap
    x <- paste0(...)
    x <- paste(strwrap(x), collapse = "\n")
    message(x)
}


# From ComplexHeatmap
complexheatmap.2 <- function(x, Rowv = TRUE, Colv = TRUE, distfun = stats::dist,
    hclustfun = stats::hclust, dendrogram = c("both", "row", "column", "none"), reorderfun = function(d,
        w) {
        stats::reorder(d, w)
    }, symm = FALSE, scale = c("none", "row", "column"), na.rm = TRUE, revC = identical(Colv,
        "Rowv"), add.expr, breaks, symbreaks = any(x < 0, na.rm = TRUE) || scale !=
        "none", col = "heat.colors", colsep, rowsep, sepcolor = "white", sepwidth = c(0.05,
        0.05), cellnote, notecex = 0.6, notecol = "cyan", na.color = graphics::par("bg"),
    trace = c("column", "row", "both", "none"), tracecol = "cyan", hline = stats::median(breaks),
    vline = stats::median(breaks), linecol = tracecol, margins = c(5, 5), ColSideColors,
    RowSideColors, cexRow = 0.6, cexCol = 0.6, labRow = NULL, labCol = NULL, srtRow = NULL,
    srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
    colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5, density.info = c("histogram",
        "density", "none"), denscol = tracecol, symkey = any(x < 0, na.rm = TRUE) ||
        symbreaks, densadj = 0.25, key.title = NULL, key.xlab = NULL, key.ylab = NULL,
    key.xtickfun = NULL, key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL,
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, ..., run_draw = FALSE) {
    if (is.data.frame(x)) {
        warning_wrap("The input is a data frame, convert it to the matrix.")
        mat <- as.matrix(x)
    } else {
        mat <- x
    }
    nr <- nrow(mat)
    nc <- ncol(mat)
    ht_param <- list()
    if (identical(Rowv, FALSE) && identical(symm, TRUE)) {
        Colv <- FALSE
    }
    if (identical(Rowv, NA) || identical(Rowv, NULL) || identical(Rowv, FALSE)) {
        ht_param$cluster_rows <- FALSE
    } else {
        if (inherits(Rowv, c("dendrogram", "hclust"))) {
            row_dend <- stats::as.dendrogram(Rowv)
        } else if (is.vector(Rowv)) {
            if (length(Rowv) == nr) {
                row_dend <- stats::as.dendrogram(hclustfun(distfun(mat)))
                row_dend <- reorderfun(row_dend, Rowv)
            } else if (length(Rowv) == 1 && is.logical(Rowv)) {
                if (Rowv) {
                  row_dend <- stats::as.dendrogram(hclustfun(distfun(mat)))
                  row_dend <- reorderfun(row_dend, -rowMeans(mat))
                }
            } else {
                stop_wrap("Wrong value for 'Rowv'.")
            }
        }
        ht_param$cluster_rows <- row_dend
    }
    if (identical(Colv, NA) || identical(Colv, NULL) || identical(Colv, FALSE)) {
        ht_param$cluster_columns <- FALSE
    } else {
        if (inherits(Colv, c("dendrogram", "hclust"))) {
            column_dend <- stats::as.dendrogram(Colv)
        } else if (is.vector(Colv)) {
            if (length(Colv) == nc) {
                column_dend <- stats::as.dendrogram(hclustfun(distfun(t(mat))))
                column_dend <- reorderfun(column_dend, Colv)
            } else if (length(Colv) == 1 && is.logical(Colv)) {
                if (Colv) {
                  column_dend <- stats::as.dendrogram(hclustfun(distfun(t(mat))))
                  column_dend <- reorderfun(column_dend, colMeans(mat))
                }
            } else {
                stop_wrap("Wrong value for 'Colv'.")
            }
        }
        ht_param$cluster_columns <- column_dend
    }
    dendrogram <- match.arg(dendrogram)[1]
    if ("both" %in% dendrogram) {
        ht_param$show_row_dend <- TRUE
        ht_param$show_column_dend <- TRUE
    } else if ("row" %in% dendrogram) {
        ht_param$show_row_dend <- TRUE
        ht_param$show_column_dend <- FALSE
    } else if ("column" %in% dendrogram) {
        ht_param$show_row_dend <- FALSE
        ht_param$show_column_dend <- TRUE
    } else if ("none" %in% dendrogram) {
        ht_param$show_row_dend <- FALSE
        ht_param$show_column_dend <- FALSE
    }
    ht_param$row_dend_width <- grid::unit(4, "cm")
    ht_param$column_dend_height <- grid::unit(3, "cm")
    scale <- match.arg(scale)[1]
    if ("row" %in% scale) {
        if (any(is.na(mat))) {
            mat <- (mat - rowMeans(mat, na.rm = TRUE))/matrixStats::rowSds(mat, na.rm = TRUE)
        } else {
            mat <- t(scale(t(mat)))
        }
        message_wrap("Note, in 'heatmap.2()', when rows are scaled, the row dendrogram is still calculated from the original matrix (not from the scaled matrix).")
    } else if ("column" %in% scale) {
        if (any(is.na(mat))) {
            mat <- t((t(mat) - colMeans(mat, na.rm = TRUE))/matrixStats::colSds(mat,
                na.rm = TRUE))
        } else {
            mat <- scale(mat)
        }
        message_wrap("Note, in 'heatmap.2()', when columns are scaled, the column dendrogram is still calculated from the original matrix (not from the scaled matrix).")
    }
    ht_param$matrix <- mat
    if (is.character(col) && length(col) == 1) {
        col <- get(col, mode = "function")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col)) {
            breaks <- 16
        } else {
            breaks <- length(col) + 1
        }
    }
    if (length(breaks) == 1) {
        if (!symbreaks) {
            breaks <- seq(min(mat, na.rm = na.rm), max(mat, na.rm = na.rm), length.out = breaks)
        } else {
            extreme <- max(abs(mat), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length.out = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (inherits(col, "function")) {
        col <- col(ncol)
    }
    n_col <- ncol
    if (exists("extreme")) {
        lim <- max(abs(mat), na.rm = TRUE)
        ht_param$col <- circlize::colorRamp2(seq(-lim, lim, length.out = n_col),
            col)
    } else {
        ht_param$col <- circlize::colorRamp2(seq(min(mat, na.rm = TRUE), max(mat,
            na.rm = TRUE), length.out = n_col), col)
    }
    if (!missing(colsep)) {
        warning_wrap("argument `colsep` is not supported in heatmap.2 -> Heatmap translation, skip it. Suggest to use `column_split` argument in Heatmap() which can be directly used here.")
    }
    if (!missing(rowsep)) {
        warning_wrap("argument `rowsep` is not supported in heatmap.2 -> Heatmap translation, skip it. Suggest to use `row_split` argument in Heatmap() which can be directly used here.")
    }
    if (!missing(cellnote)) {
        ht_param$layer_fun <- function(j, i, x, y, w, h, fill) {
            grid::grid.text(ComplexHeatmap::pindex(cellnote, i, j), x, y, gp = grid::gpar(cex = notecex,
                col = notecol))
        }
    }
    ht_param$na_col <- na.color
    if (!missing(ColSideColors)) {
        ht_param$top_annotation <- ComplexHeatmap::HeatmapAnnotation(column = ColSideColors,
            col = list(column = structure(unique(ColSideColors), names = unique(ColSideColors))),
            show_legend = FALSE, show_annotation_name = FALSE)
    }
    if (!missing(RowSideColors)) {
        ht_param$left_annotation <- function(...) {
            x <- paste0(...)
            x <- paste(strwrap(x), collapse = "\n")
            stop(x, call. = FALSE)
        }
        ComplexHeatmap::rowAnnotation(row = RowSideColors, col = list(row = structure(unique(RowSideColors),
            names = unique(RowSideColors))), show_legend = FALSE, show_annotation_name = FALSE)
    }
    if (identical(ht_param$cluster_rows, FALSE) || dendrogram %in% c("none", "column")) {
        if (is.null(ht_param$left_annotation)) {
            ht_param$left_annotation <- ComplexHeatmap::rowAnnotation(foo1 = ComplexHeatmap::anno_empty(width = grid::unit(4,
                "cm"), border = FALSE))
        } else {
            ht_param$left_annotation <- c(ht_param$left_annotation, ComplexHeatmap::rowAnnotation(foo1 = ComplexHeatmap::anno_empty(width = grid::unit(3.5,
                "cm"))))
        }
    }
    if (identical(ht_param$cluster_columns, FALSE) || dendrogram %in% c("none", "row")) {
        if (is.null(ht_param$top_annotation)) {
            ht_param$top_annotation <- ComplexHeatmap::HeatmapAnnotation(foo2 = ComplexHeatmap::anno_empty(height = grid::unit(3,
                "cm"), border = FALSE))
        } else {
            ht_param$top_annotation <- c(ht_param$top_annotation, ComplexHeatmap::HeatmapAnnotation(foo2 = ComplexHeatmap::anno_empty(height = grid::unit(3.5,
                "cm"))))
        }
    }
    if (!is.null(labRow)) {
        ht_param$row_labels <- labRow
    } else if (is.null(rownames(mat))) {
        ht_param$row_labels <- seq_len(nr)
    } else {
        ht_param$row_labels <- rownames(mat)
    }
    if (!is.null(labCol)) {
        ht_param$column_labels <- labCol
    } else if (is.null(colnames(mat))) {
        ht_param$column_labels <- seq_len(nc)
    } else {
        ht_param$column_labels <- colnames(mat)
    }
    if (is.null(colRow)) {
        colRow <- "black"
    }
    ht_param$row_names_gp <- grid::gpar(fontsize = 12 * cexRow, col = colRow)
    if (is.null(colCol)) {
        colCol <- "black"
    }
    ht_param$column_names_gp <- grid::gpar(fontsize = 12 * cexCol, col = colCol)
    if (!is.null(srtRow)) {
        ht_param$row_names_rot <- srtRow
    } else {
        ht_param$row_names_rot <- 0
    }
    if (!is.null(srtCol)) {
        ht_param$column_names_rot <- srtCol
    } else {
        ht_param$column_names_rot <- 90
    }
    if (!is.null(main)) {
        ht_param$column_title <- main
    }
    if (!is.null(xlab)) {
        if (is.null(ht_param$column_labels)) {
            ht_param$bottom_annotation <- ComplexHeatmap::HeatmapAnnotation(xlab = ComplexHeatmap::anno_block(labels = xlab,
                gp = grid::gpar(col = NA)))
        } else {
            ht_param$bottom_annotation <- ComplexHeatmap::HeatmapAnnotation(colnames = ComplexHeatmap::anno_text(ht_param$column_labels,
                gp = ht_param$column_names_gp, rot = ht_param$column_names_rot),
                xlab = ComplexHeatmap::anno_block(labels = xlab, gp = grid::gpar(col = NA)))
            ht_param$show_column_names <- FALSE
        }
    }
    if (!is.null(ylab)) {
        if (is.null(ht_param$row_labels)) {
            ht_param$right_annotation <- ComplexHeatmap::rowAnnotation(ylab = ComplexHeatmap::anno_block(labels = ylab,
                gp = grid::gpar(col = NA)))
        } else {
            ht_param$right_annotation <- ComplexHeatmap::rowAnnotation(rownames = ComplexHeatmap::anno_text(ht_param$row_labels,
                gp = ht_param$row_names_gp, rot = ht_param$row_names_rot), ylab = ComplexHeatmap::anno_block(labels = ylab,
                gp = grid::gpar(col = NA)))
            ht_param$show_row_names <- FALSE
        }
    }
    ht_param$show_heatmap_legend <- FALSE
    trace <- match.arg(trace)[1]
    min <- breaks[1]
    max <- breaks[length(breaks)]
    rg <- max - min
    layer_fun <- NULL
    if (trace == "both") {
        warning_wrap("trace = 'both' is not supported, change to trace = 'column'.")
        trace <- "column"
    }
    if (trace == "column") {
        layer_fun <- function(j, i, x, y, w, h, fill) {
            ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
            ind <- ind_mat[1, ]
            grid::grid.segments(x[ind], grid::unit(0, "npc"), x[ind], grid::unit(1,
                "npc"), gp = grid::gpar(col = linecol, lty = 2))
            for (ki in seq_len(ncol(ind_mat))) {
                ind <- ind_mat[, ki]
                offset <- (mat[i[ind], j[ind[1]]] - min)/(max - min) * w[ind]
                pos_x <- rep(x[ind] - w[ind] * 0.5 + offset, each = 2)
                pos_y <- rep(y[ind] + h[ind] * 0.5, each = 2)
                pos_y[seq_along(pos_y)%%2 == 0] <- y[ind] - h[ind] * 0.5
                grid::grid.lines(pos_x, pos_y, gp = grid::gpar(col = tracecol))
            }
        }
    } else if (trace == "row") {
        layer_fun <- function(j, i, x, y, w, h, fill) {
            ind_mat <- ComplexHeatmap::restore_matrix(j, i, x, y)
            ind <- ind_mat[, 1]
            grid::grid.segments(grid::unit(0, "npc"), y[ind], grid::unit(1, "npc"),
                y[ind], gp = grid::gpar(col = linecol, lty = 2))
            for (ki in seq_len(nrow(ind_mat))) {
                ind <- ind_mat[ki, ]
                offset <- (mat[i[ind[1]], j[ind]] - min)/(max - min) * h[ind]
                pos_x <- rep(x[ind] - w[ind] * 0.5, each = 2)
                pos_x[seq_along(pos_x)%%2 == 0] <- x[ind] + w[ind] * 0.5
                pos_y <- rep(y[ind] - h[ind] * 0.5 + offset, each = 2)
                grid::grid.lines(pos_x, pos_y, gp = grid::gpar(col = tracecol))
            }
        }
    }
    if (!is.null(ht_param$layer_fun) && trace %in% c("row", "column")) {
        fun1 <- ht_param$layer_fun
        fun2 <- layer_fun
        fun3 <- function(j, i, x, y, w, h, fill) {
            fun1(j, i, x, y, w, h, fill)
            fun2(j, i, x, y, w, h, fill)
        }
        layer_fun <- fun3
        ht_param$layer_fun <- layer_fun
    } else if (is.null(ht_param$layer_fun) && trace %in% c("row", "column")) {
        ht_param$layer_fun <- layer_fun
    }
    random_str <- paste(sample(c(letters, LETTERS, 0:9), 8), collapse = "")
    ht_param$name <- paste0("heatmap.2_", random_str)
    density.info <- match.arg(density.info)[1]
    if (is.null(key.xlab)) {
        key.xlab <- "Value"
    }
    if (is.null(key.ylab)) {
        if (density.info == "histogram") {
            key.ylab <- "Count"
        }
        if (density.info == "density") {
            key.ylab <- "Density"
        }
    }
    if (density.info == "none") {
        post_fun <- NULL
    } else {
        post_fun <- function(ht) {
            ComplexHeatmap::decorate_heatmap_body(paste0("heatmap.2_", random_str),
                {
                  grid::pushViewport(grid::viewport(grid::unit(0, "npc"), grid::unit(1,
                    "npc"), width = grid::unit(4, "cm"), height = grid::unit(3, "cm"),
                    just = c("right", "bottom")))
                  left_width <- grid::unit(1, "cm")
                  bottom_height <- grid::unit(1, "cm")
                  top_height <- grid::unit(0.5, "cm")
                  if (density.info == "histogram") {
                    tb <- graphics::hist(mat, breaks = breaks, plot = FALSE)
                    x_at <- pretty(range(tb$breaks), n = 3)
                    y_at <- pretty(range(tb$counts), n = 3)
                    x_range <- range(tb$breaks)
                    y_range <- range(tb$counts)
                    y_range[2] <- y_range[2] + (y_range[2] - y_range[1]) * 0.05
                  } else if (density.info == "density") {
                    den <- stats::density(mat, na.rm = TRUE, from = min(breaks),
                      to = max(breaks), adjust = densadj)
                    den_x <- den$x
                    den_y <- den$y
                    x_range <- range(breaks)
                    l <- den_x >= x_range[1] & den_x <= x_range[2]
                    den_x <- den_x[l]
                    den_y <- den_y[l]
                    x_at <- pretty(range(den_x), n = 3)
                    y_at <- pretty(range(den_y), n = 3)
                    y_range <- range(den_y)
                    y_range[2] <- y_range[2] + (y_range[2] - y_range[1]) * 0.05
                  }
                  x_at <- x_at[x_at >= x_range[1] & x_at <= x_range[2]]
                  y_at <- y_at[y_at >= y_range[1] & y_at <= y_range[2]]
                  grid::pushViewport(grid::viewport(x = left_width, y = bottom_height,
                    width = grid::unit(1, "npc") - left_width - grid::unit(2, "mm"),
                    height = grid::unit(1, "npc") - bottom_height - top_height, just = c("left",
                      "bottom"), xscale = x_range, yscale = y_range))
                  x <- seq(min(breaks), max(breaks), length.out = 101)
                  grid::grid.rect(x = x[1:100], width = (x_range[2] - x_range[1])/100,
                    default.units = "native", just = "left", gp = grid::gpar(fill = ht_param$col(x +
                      (x_range[2] - x_range[1])/100 * 0.5), col = NA))
                  if (density.info == "histogram") {
                    x <- rep(tb$breaks, each = 2)
                    x <- x[-c(1, length(x))]
                    y <- rep(tb$counts, each = 2)
                  } else if (density.info == "density") {
                    x <- den_x
                    y <- den_y
                  }
                  grid::grid.lines(x, y, default.units = "native", gp = grid::gpar(col = denscol))
                  if (trace == "column") {
                    grid::grid.lines(c(vline, vline), grid::unit(c(0, 1), "npc"),
                      default.units = "native", gp = grid::gpar(col = linecol, lty = 2))
                  } else if (trace == "row") {
                    grid::grid.lines(grid::unit(c(0, 1), "npc"), c(hline, hline),
                      default.units = "native", gp = grid::gpar(col = linecol, lty = 2))
                  }
                  grid::grid.rect(gp = grid::gpar(fill = "transparent"))
                  grid::grid.xaxis(at = x_at, gp = grid::gpar(fontsize = 8))
                  grid::grid.text(key.xlab, y = grid::unit(0, "npc") - grid::unit(8,
                    "mm"), gp = grid::gpar(fontsize = 8))
                  grid::grid.yaxis(at = y_at, gp = grid::gpar(fontsize = 8))
                  grid::grid.text(key.ylab, x = grid::unit(0, "npc") - grid::unit(8,
                    "mm"), gp = grid::gpar(fontsize = 8), rot = 90)
                  grid::grid.text("Color Key", y = grid::unit(1, "npc") + top_height *
                    0.5, gp = grid::gpar(fontface = "bold", fontsize = 10))
                  grid::popViewport()
                  grid::popViewport()
                })
        }
    }
    ht_param$post_fun <- post_fun
    if (!missing(add.expr)) {
        warning_wrap("argument `add.expr` is not supported in heatmap.2 -> Heatmap translation, skip it.")
    }
    if (!is.null(lmat)) {
        warning_wrap("argument `lmat` is not supported in heatmap.2 -> Heatmap translation, skip it.")
    }
    if (!is.null(lhei)) {
        warning_wrap("argument `lhei` is not supported in heatmap.2 -> Heatmap translation, skip it.")
    }
    if (!is.null(lwid)) {
        warning_wrap("argument `lwid` is not supported in heatmap.2 -> Heatmap translation, skip it.")
    }
    if (!is.null(extrafun)) {
        warning_wrap("argument `extrafun` is not supported in heatmap.2 -> Heatmap translation, skip it.")
    }
    ht_param <- c(ht_param, list(...))
    ht <- do.call(ComplexHeatmap::Heatmap, ht_param)
    attr(ht, "translate_from") <- "heatmap"
    if (run_draw) {
        ComplexHeatmap::draw(ht)
    } else {
        ht
    }
}
