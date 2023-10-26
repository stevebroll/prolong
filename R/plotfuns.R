delta_scatter <- function(x,
                          timediff1 = "2-1",
                          timediff2 = "3-2",
                          timediff3 = NULL,
                          fisherz = T,
                          interactive = T,
                          digits = 3) {
  p <- ncol(x)

  t1 <- unlist(strsplit(timediff1, "-"))[1]
  t2 <- unlist(strsplit(timediff1, "-"))[2]
  x1 <- x[, , t2] - x[, , t1]

  t1 <- unlist(strsplit(timediff2, "-"))[1]
  t2 <- unlist(strsplit(timediff2, "-"))[2]
  x2 <- x[, , t2] - x[, , t1]

  corr1 <- stats::cor(x1)
  corr2 <- stats::cor(x2)

  if (fisherz) {
    corr1 <- atanh(corr1)
    corr2 <- atanh(corr2)
  }

  vec1 <- corr1[upper.tri(corr1)]
  vec2 <- corr2[upper.tri(corr2)]

  labs <-
    expand.grid(colnames(x), colnames(x))[as.vector(upper.tri(matrix(0, p, p))), ]
  if (is.null(timediff3)) {
    # get long mat
    rmat <-
      data.frame(round(vec1, digits), round(vec2, digits), labs)
    if (fisherz) {
      colnames(rmat) <-
        c(
          paste(timediff1, "Fisher Z-Transformed Correlations"),
          paste(timediff2, "Fisher Z-Transformed Correlations"),
          "Variable 1",
          "Variable 2"
        )
    } else {
      c(
        paste(timediff1, "Correlations"),
        paste(timediff2, "Correlations"),
        "Variable 1",
        "Variable 2"
      )
    }
    name <- colnames(rmat)[1:2]
    p <- ggplot2::ggplot(data = rmat, ggplot2::aes(
      x = rlang::.data[[name[1]]],
      y = rlang::.data[[name[2]]],
      text = c(paste(
        "Variable 1:", rmat[, 3], "\nVariable 2:", rmat[, 4]
      ))
    )) +
      ggplot2::geom_point(size = .5)
    if (interactive) {
      p <- p %>% plotly::toWebGL()
      plotly::ggplotly(p) %>% plotly::layout(hoverlabel = list(align = "left"))
    } else {
      p
    }
  } else {
    t1 <- unlist(strsplit(timediff3, "-"))[1]
    t2 <- unlist(strsplit(timediff3, "-"))[2]
    x3 <- x[, , t2] - x[, , t1]

    corr3 <- stats::cor(x3)

    if (fisherz) {
      corr3 <- atanh(corr3)
    }

    vec3 <- corr3[upper.tri(corr3)]

    # get long mat
    rmat <-
      data.frame(
        round(vec1, digits),
        round(vec2, digits),
        round(vec3, digits),
        labs
      )
    if (fisherz) {
      colnames(rmat) <-
        c(
          paste(timediff1, "Fisher Z-Transformed Correlations"),
          paste(timediff2, "Fisher Z-Transformed Correlations"),
          paste(timediff3, "Fisher Z-Transformed Correlations"),
          "Variable 1",
          "Variable 2"
        )
    } else {
      c(
        paste(timediff1, "Correlations"),
        paste(timediff2, "Correlations"),
        paste(timediff3, "Correlations"),
        "Variable 1",
        "Variable 2"
      )
    }
    name <- colnames(rmat)[1:3]
    p <- ggplot2::ggplot(data = rmat, ggplot2::aes(
      x = rlang::.data[[name[1]]],
      y = rlang::.data[[name[2]]],
      text = c(paste(
        "Variable 1:", rmat[, 3], "\nVariable 2:", rmat[, 4]
      ))
    )) +
      ggplot2::geom_point(size = .5)
    if (interactive) {
      p <- plotly::plot_ly(
        rmat,
        x = ~ get(name[1]),
        y = ~ get(name[2]),
        z = ~ get(name[3]),
        hovertemplate = paste(
          timediff1,
          " Corr: %{x}<br>",
          timediff2,
          " Corr: %{y}<br>",
          timediff3,
          " Corr: %{z}<br>",
          "%{text}",
          "<extra></extra>",
          sep = ""
        ),
        text = c(paste(
          "Variable 1:", rmat[, 4], "\nVariable 2:", rmat[, 5]
        ))
      )
      p <- p %>% plotly::add_markers(size = 1)
      p <- p %>% plotly::partial_bundle()
      p
    } else {
      p <- plotly::plot_ly(
        rmat,
        x = ~ get(name[1]),
        y = ~ get(name[2]),
        z = ~ get(name[3])
      )
      p <- p %>% plotly::add_markers(size = 1)
      p <- p %>% plotly::partial_bundle()
      p
    }
  }
}
