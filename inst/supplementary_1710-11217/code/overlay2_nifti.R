overlay2.nifti <- function (x, y, z = 1, w = 1, col.x = gray(0:64/64), col.y = hotmetal(),
    zlim.x = NULL, zlim.y = NULL, plane = c("axial", "coronal",
        "sagittal"), plot.type = c("multiple", "single"), xlab = "",
    ylab = "", axes = FALSE, oma = rep(0, 4), mar = rep(0, 4),
    bg = "black", NA.x = FALSE, NA.y = FALSE, title = "", ...)
{
    switch(plane[1], axial = {
        aspect <- x@pixdim[3]/x@pixdim[2]
    }, coronal = {
        if (length(dim(x)) == 3) {
            x@.Data <- aperm(x, c(1, 3, 2))
        } else {
            x@.Data <- aperm(x, c(1, 3, 2, 4))
        }
        y@.Data <- aperm(y, c(1, 3, 2))
        aspect <- x@pixdim[4]/x@pixdim[2]
    }, sagittal = {
        if (length(dim(x)) == 3) {
            x@.Data <- aperm(x, c(2, 3, 1))
        } else {
            x@.Data <- aperm(x, c(2, 3, 1, 4))
        }
        y@.Data <- aperm(y, c(2, 3, 1))
        aspect <- x@pixdim[4]/x@pixdim[3]
    }, stop(paste("Orthogonal plane", plane[1], "is not valid.")))
    if (!all(dim(x)[1:3] == dim(y)[1:3])) {
        stop("dimensions of \"x\" and \"y\" must be equal")
    }
    if (NA.x) {
        x[x == 0] = NA
        if (all(is.na(x))) {
            stop(paste0("x has no non-zero values and NA.x = TRUE.  ",
                "Likely set NA.x = FALSE."))
        }
    }
    if (NA.y) {
        y[y == 0] = NA
        if (all(is.na(y))) {
            stop(paste0("y has no non-zero values and NA.y = TRUE.  ",
                "Either remove the overlay, or set NA.y = FALSE"))
        }
    }
    X <- nrow(x)
    Y <- ncol(x)
    Z <- nsli(x)
    W <- ntim(x)
    if (X == 0 || Y == 0 || Z == 0) {
        stop("size of NIfTI volume is zero, nothing to plot")
    }
    if (is.null(zlim.x)) {
        zlim.x <- c(x@cal_min, x@cal_max)
        if (any(!is.finite(zlim.x)) || diff(zlim.x) == 0) {
            zlim.x <- c(x@glmin, x@glmax)
        }
        if (any(!is.finite(zlim.x)) || diff(zlim.x) == 0) {
            zlim.x <- range(x, na.rm = TRUE)
        }
    }
    breaks.x <- c(min(x, zlim.x, na.rm = TRUE), seq(min(zlim.x,
        na.rm = TRUE), max(zlim.x, na.rm = TRUE), length = length(col.x) -
        1), max(x, zlim.x, na.rm = TRUE))
    if (is.null(zlim.y)) {
        zlim.y <- c(y@cal_min, y@cal_max)
        if (any(!is.finite(zlim.y)) || diff(zlim.y) == 0) {
            zlim.y <- c(y@glmin, y@glmax)
        }
        if (any(!is.finite(zlim.y)) || diff(zlim.y) == 0) {
            zlim.y <- range(y, na.rm = TRUE)
        }
    }
    if (plot.type[1] == "multiple") {
        index <- 1:Z
    }
    else {
        index <- z
    }
    lz <- length(index)
    if (z < 1 || z > Z) {
        stop("slice \"z\" out of range")
    }
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = ceiling(rep(sqrt(lz), 2)), oma = oma, mar = mar,
        bg = bg)
    if (is.na(W)) {
        for (z in index) {
            graphics::image(1:X, 1:Y, x[, , z], col = col.x,
                breaks = breaks.x, zlim = zlim.x, asp = aspect,
                axes = axes, xlab = xlab, ylab = ylab, ...)
            graphics::image(1:X, 1:Y, y[, , z], col = col.y,
                            zlim = zlim.y, add = TRUE)
        }
    }
    else {
        if (w < 1 || w > W) {
            stop("volume \"w\" out of range")
        }
        for (z in index) {
            graphics::image(1:X, 1:Y, x[, , z, w], col = col.x,
                breaks = breaks.x, zlim = zlim.x, asp = aspect,
                axes = axes, xlab = xlab, ylab = ylab, ...)
            graphics::image(1:X, 1:Y, y[, , z], col = col.y,
                            zlim = zlim.y, add = TRUE)
        }
    }
    xu <- X/100
    yu <- Y/100
    ## Assumes min and max values of y are -u and u
    plotrix::color.legend(xl = X - 6*xu, xr = X - 2*xu, yb = Y - 40*yu, yt = Y - 10*yu,
                          legend = c(min(y, na.rm = TRUE), max(y, na.rm = TRUE)),
                          rect.col = col.y, gradient = "y", col = "white")
    mtext(title, col = "white", line = -2)
    par(oldpar)
    invisible()
}
