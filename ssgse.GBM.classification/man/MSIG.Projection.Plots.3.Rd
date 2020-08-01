\name{MSIG.Projection.Plots.3}
\alias{MSIG.Projection.Plots.3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Projection.Plots.3
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Projection.Plots.3(input.ds, input.cls = "", model.set = "ALL", output.2D.proj.file, output.2D.proj.plot, output.3D.proj.file, output.3D.1.proj.plot, output.3D.2.proj.plot, output.3D.3.proj.plot, output.heatmap.plot, output.heatmap.sorted.plot, output.heatmap.sorted.2.plot, output.hclust.plot, use.feature.names = FALSE, use.biplot = TRUE, title = "", seed = 1234, non.interactive.run = F, heatmap.row.norm = T, heatmap.cmap.type = 1, symbol.scaling = 1, col = c("greny3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"), symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{model.set}{
%%     ~~Describe \code{model.set} here~~
}
  \item{output.2D.proj.file}{
%%     ~~Describe \code{output.2D.proj.file} here~~
}
  \item{output.2D.proj.plot}{
%%     ~~Describe \code{output.2D.proj.plot} here~~
}
  \item{output.3D.proj.file}{
%%     ~~Describe \code{output.3D.proj.file} here~~
}
  \item{output.3D.1.proj.plot}{
%%     ~~Describe \code{output.3D.1.proj.plot} here~~
}
  \item{output.3D.2.proj.plot}{
%%     ~~Describe \code{output.3D.2.proj.plot} here~~
}
  \item{output.3D.3.proj.plot}{
%%     ~~Describe \code{output.3D.3.proj.plot} here~~
}
  \item{output.heatmap.plot}{
%%     ~~Describe \code{output.heatmap.plot} here~~
}
  \item{output.heatmap.sorted.plot}{
%%     ~~Describe \code{output.heatmap.sorted.plot} here~~
}
  \item{output.heatmap.sorted.2.plot}{
%%     ~~Describe \code{output.heatmap.sorted.2.plot} here~~
}
  \item{output.hclust.plot}{
%%     ~~Describe \code{output.hclust.plot} here~~
}
  \item{use.feature.names}{
%%     ~~Describe \code{use.feature.names} here~~
}
  \item{use.biplot}{
%%     ~~Describe \code{use.biplot} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{heatmap.row.norm}{
%%     ~~Describe \code{heatmap.row.norm} here~~
}
  \item{heatmap.cmap.type}{
%%     ~~Describe \code{heatmap.cmap.type} here~~
}
  \item{symbol.scaling}{
%%     ~~Describe \code{symbol.scaling} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{symbs}{
%%     ~~Describe \code{symbs} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (input.ds, input.cls = "", model.set = "ALL", output.2D.proj.file, 
    output.2D.proj.plot, output.3D.proj.file, output.3D.1.proj.plot, 
    output.3D.2.proj.plot, output.3D.3.proj.plot, output.heatmap.plot, 
    output.heatmap.sorted.plot, output.heatmap.sorted.2.plot, 
    output.hclust.plot, use.feature.names = FALSE, use.biplot = TRUE, 
    title = "", seed = 1234, non.interactive.run = F, heatmap.row.norm = T, 
    heatmap.cmap.type = 1, symbol.scaling = 1, col = c("greny3", 
        "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", 
        "orange", "blue2", "lightblue", "pink", "coral"), symbs = c(22, 
        21, 20, 23, 24, 25, 21, 20, 23, 24, 25)) 
{
    print(c("Running MSIG.Projection.Plots... on:", input.ds))
    library("scatterplot3d")
    library(MASS)
    set.seed(seed = seed, kind = NULL)
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dim(m)
    Ns <- length(m[1, ])
    k.proj <- length(m[, 1])
    if (input.cls != "") {
        CLS <- ReadClsFile(file = input.cls)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
    }
    else {
        class.labels <- rep(1, Ns)
        class.phen <- "Samples"
    }
    if (model.set == "ALL") {
        model.set <- seq(1, Ns)
    }
    m.train <- as.matrix(m[, model.set])
    num.samples.train <- length(model.set)
    sample.names.train <- sample.names[model.set]
    if (input.cls != "") {
        class.labels.train <- class.labels[model.set]
    }
    m.test <- as.matrix(m[, -model.set])
    sample.names.test <- sample.names[-model.set]
    if (input.cls != "") {
        class.labels.test <- class.labels[-model.set]
    }
    pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)
    S1 <- pca$x[, 1]
    S2 <- pca$x[, 2]
    S3 <- pca$x[, 3]
    X1 <- pca$rotation[, 1]
    X2 <- pca$rotation[, 2]
    X3 <- pca$rotation[, 3]
    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1/max.X
    X2 <- max.S * X2/max.X
    max.A <- max(max.S, max.X)
    c0 <- col
    c1 <- col
    color <- c1[class.labels]
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.2D.proj.plot
            x11(height = 20, width = 30)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.2D.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 30)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.2D.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 30)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.2D.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 30)
        }
    }
    nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(3, 
        3, 1), heights = 1, respect = FALSE)
    plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), 
        type = "n", main = paste(title, " -- Model Samples Biplot", 
            sep = ""), sub = input.ds)
    for (j in 1:num.samples.train) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
            symb <- symbs[class.labels[j] + 1]
        }
        else {
            color.code <- c1[class.labels[j]]
            symb <- symbs[class.labels[j]]
        }
        points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 
            3, bg = color.code, col = "black")
    }
    if (use.biplot == TRUE) {
        for (j in 1:k.proj) {
            x.coor <- X1[j] * 0.925
            y.coor <- X2[j] * 0.925
            arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, 
                angle = 20, col = "black")
            if (use.feature.names == FALSE) {
                leg.txt <- paste("F", j, sep = "")
            }
            else {
                leg.txt <- gs.names[j]
            }
            text(X1[j], X2[j], labels = leg.txt, adj = NULL, 
                pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 
                  2, col = "black")
        }
        ang <- vector(length = k.proj, mode = "numeric")
        for (j in 1:k.proj) {
            ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], 
                X1[j]), 2 * pi + atan2(X2[j], X1[j]))
        }
        ang.index <- order(ang, decreasing = F)
        ang2 <- ang[ang.index]
        for (j in 1:k.proj) {
            if (j == k.proj) {
                angle.in.between <- (ang2[1] - ang2[j] - 2 * 
                  pi)/2 + ang2[j] - 2 * pi
            }
            else {
                angle.in.between <- (ang2[j + 1] - ang2[j])/2 + 
                  ang2[j]
            }
            x <- max.S * cos(angle.in.between)
            y <- max.S * sin(angle.in.between)
            arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, 
                col = "grey50")
        }
    }
    test.scores <- predict(pca, t(m.test))
    S1 <- c(pca$x[, 1], test.scores[, 1])
    S2 <- c(pca$x[, 2], test.scores[, 2])
    S3 <- c(pca$x[, 3], test.scores[, 3])
    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1/max.X
    X2 <- max.S * X2/max.X
    num.samples <- length(S1)
    plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), 
        type = "n", main = paste(title, " -- All Samples Biplot", 
            sep = ""), sub = input.ds)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            symb <- symbs[class.labels[j] + 1]
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            symb <- symbs[class.labels[j]]
            color.code <- c1[class.labels[j]]
        }
        points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 
            3, bg = color.code, col = "black")
    }
    if (use.biplot == TRUE) {
        for (j in 1:k.proj) {
            x.coor <- X1[j] * 0.925
            y.coor <- X2[j] * 0.925
            arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, 
                angle = 20, col = "black")
            if (use.feature.names == FALSE) {
                leg.txt <- paste("F", j, sep = "")
            }
            else {
                leg.txt <- gs.names[j]
            }
            text(X1[j], X2[j], labels = leg.txt, adj = NULL, 
                pos = NULL, offset = 0.5, vfont = NULL, cex = symbol.scaling * 
                  2, col = "black")
        }
        ang <- vector(length = k.proj, mode = "numeric")
        for (j in 1:k.proj) {
            ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], 
                X1[j]), 2 * pi + atan2(X2[j], X1[j]))
        }
        ang.index <- order(ang, decreasing = F)
        ang2 <- ang[ang.index]
        for (j in 1:k.proj) {
            if (j == k.proj) {
                angle.in.between <- (ang2[1] - ang2[j] - 2 * 
                  pi)/2 + ang2[j] - 2 * pi
            }
            else {
                angle.in.between <- (ang2[j + 1] - ang2[j])/2 + 
                  ang2[j]
            }
            x <- max.S * cos(angle.in.between)
            y <- max.S * sin(angle.in.between)
            arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, 
                col = "grey50")
        }
    }
    class.phen.train <- unique(class.labels.train)
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = symbol.scaling * 1.35, pt.cex = symbol.scaling * 
            3)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    max.S <- max(sqrt(S1 * S1 + S2 * S2 + S3 * S3))
    max.X <- max(sqrt(X1 * X1 + X2 * X2 + X3 * X3))
    X1 <- max.S * X1/max.X
    X2 <- max.S * X2/max.X
    X3 <- max.S * X3/max.X
    max.A <- max(max.S, max.X)
    color <- c1[class.labels]
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.3D.1.proj.plot
            x11(height = 20, width = 25)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.1.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.1.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.3D.1.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 
        1), heights = 1, respect = FALSE)
    x <- scatterplot3d(S1, S2, S3, xlab = "F1", ylab = "F2", 
        zlab = "F3", type = "n", angle = 45, pch = 20, main = paste(title, 
            " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = symbol.scaling * 
            1)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
            symb <- symbs[class.labels[j] + 1]
        }
        else {
            symb <- symbs[class.labels[j]]
            color.code <- c1[class.labels[j]]
        }
        cex <- 2.5 * (max(S2) - S2[j])/(max(S2) - min(S2)) + 
            1.5
        x$points3d(S1[j], S2[j], S3[j], col = "black", pch = symb, 
            bg = color.code, cex = 1 * symbol.scaling * cex)
    }
    if (use.biplot == TRUE) {
        origin.3D <- x$xyz.convert(0, 0, 0)
        for (j in 1:k.proj) {
            x.coor <- X1[j] * 0.925
            y.coor <- X2[j] * 0.925
            z.coor <- X3[j] * 0.925
            end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
            arrows(origin.3D$x, origin.3D$y, end.point.3D$x, 
                end.point.3D$y, lwd = 2, length = 0.15, angle = 20, 
                col = "grey")
            if (use.feature.names == FALSE) {
                leg.txt <- paste("F", j, sep = "")
            }
            else {
                leg.txt <- gs.names[j]
            }
            x.coor <- X1[j]
            y.coor <- X2[j]
            z.coor <- X3[j]
            end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
            text(end.point.3D$x, end.point.3D$y, labels = leg.txt, 
                adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, 
                cex = symbol.scaling * 2, col = "grey")
        }
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = symbol.scaling * 1.2, pt.cex = symbol.scaling * 
            1.5)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    S3 <- -S3
    X3 <- -X3
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.3D.2.proj.plot
            x11(height = 20, width = 25)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.2.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.2.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.3D.2.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 
        1), heights = 1, respect = FALSE)
    x <- scatterplot3d(S1, S3, S2, xlab = "F1", ylab = "F3", 
        zlab = "F2", type = "n", angle = 45, pch = 20, main = paste(title, 
            " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = symbol.scaling * 
            1)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            symb <- symbs[class.labels[j] + 1]
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            symb <- symbs[class.labels[j]]
            color.code <- c1[class.labels[j]]
        }
        cex <- 2.5 * (max(S3) - S3[j])/(max(S3) - min(S3)) + 
            1.5
        x$points3d(S1[j], S3[j], S2[j], col = "black", pch = symb, 
            bg = color.code, cex = 1 * symbol.scaling * cex)
    }
    if (use.biplot == TRUE) {
        origin.3D <- x$xyz.convert(0, 0, 0)
        for (j in 1:k.proj) {
            x.coor <- X1[j] * 0.925
            y.coor <- X3[j] * 0.925
            z.coor <- X2[j] * 0.925
            end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
            arrows(origin.3D$x, origin.3D$y, end.point.3D$x, 
                end.point.3D$y, lwd = 2, length = 0.15, angle = 20, 
                col = "grey")
            if (use.feature.names == FALSE) {
                leg.txt <- paste("F", j, sep = "")
            }
            else {
                leg.txt <- gs.names[j]
            }
            x.coor <- X1[j]
            y.coor <- X3[j]
            z.coor <- X2[j]
            end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
            text(end.point.3D$x, end.point.3D$y, labels = leg.txt, 
                adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, 
                cex = symbol.scaling * 2, col = "grey")
        }
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = symbol.scaling * 1.2, pt.cex = symbol.scaling * 
            1.5)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    S1 <- -S1
    X1 <- -X1
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.3D.3.proj.plot
            x11(height = 20, width = 25)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.3.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.3.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.3D.3.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 
        1), heights = 1, respect = FALSE)
    x <- scatterplot3d(S2, S1, S3, xlab = "F2", ylab = "F1", 
        zlab = "F3", type = "n", angle = 45, pch = 20, main = paste(title, 
            " -- 3D Biplot", sep = ""), sub = " ", cex.symbols = symbol.scaling * 
            1)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            symb <- symbs[class.labels[j] + 1]
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            symb <- symbs[class.labels[j]]
            color.code <- c1[class.labels[j]]
        }
        cex <- 2.5 * (max(S1) - S1[j])/(max(S1) - min(S1)) + 
            1.5
        x$points3d(S2[j], S1[j], S3[j], col = "black", pch = symb, 
            bg = color.code, cex = 1 * symbol.scaling * cex)
    }
    if (use.biplot == TRUE) {
        origin.3D <- x$xyz.convert(0, 0, 0)
        for (j in 1:k.proj) {
            x.coor <- X2[j] * 0.925
            y.coor <- X1[j] * 0.925
            z.coor <- X3[j] * 0.925
            end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
            arrows(origin.3D$x, origin.3D$y, end.point.3D$x, 
                end.point.3D$y, lwd = 2, length = 0.15, angle = 20, 
                col = "grey")
            if (use.feature.names == FALSE) {
                leg.txt <- paste("F", j, sep = "")
            }
            else {
                leg.txt <- gs.names[j]
            }
            x.coor <- X2[j]
            y.coor <- X1[j]
            z.coor <- X3[j]
            end.point.3D <- x$xyz.convert(x.coor, y.coor, z.coor)
            text(end.point.3D$x, end.point.3D$y, labels = leg.txt, 
                adj = NULL, pos = NULL, offset = 0.5, vfont = NULL, 
                cex = symbol.scaling * 2, col = "grey")
        }
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = symbol.scaling * 1.2, pt.cex = 1.5)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    height <- ifelse(k.proj > 50, 20, 0.5 * k.proj + 7)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.heatmap.plot
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
    }
    MSIG.HeatMapPlot.5(V = m, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, col.names = sample.names, main = paste(title, 
            " -- Heat Map ", sep = ""), sub = " ", xlab = " ", 
        ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, 
        char.rescale = 1)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    dist.matrix <- dist(t(m))
    HC <- hclust(dist.matrix, method = "complete")
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.hclust.plot
            x11(height = 20, width = 25)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.hclust.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.hclust.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.hclust.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 20, width = 25)
        }
    }
    HC$labels <- class.phen[class.labels]
    dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2)
    local({
        colLab <<- function(n) {
            if (is.leaf(n)) {
                a <- attributes(n)
                i <<- i + 1
                attr(n, "nodePar") <- c(a$nodePar, list(lab.col = mycols[i], 
                  pch = c(0, 0), col = c(0, 0), bg = c(0, 0), 
                  cex = c(0.8, 0.8), lab.font = i\%\%1))
            }
            n
        }
        mycols <- col[class.labels[HC$order]]
        i <- 0
    })
    dL <- dendrapply(dhc, colLab)
    plot(dL, cex = 1, edge.root = T)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    dist.matrix <- dist(t(m))
    HC <- hclust(dist.matrix, method = "complete")
    m2 <- m[, HC$order]
    sample.names2 <- sample.names[HC$order]
    class.labels2 <- class.labels[HC$order]
    dist.matrix <- dist(m)
    HC <- hclust(dist.matrix, method = "complete")
    m2 <- m2[HC$order, ]
    gs.names2 <- gs.names[HC$order]
    height <- ifelse(k.proj > 50, 20, 0.5 * k.proj + 7)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.heatmap.sorted.plot
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.sorted.plot, 
                ".pdf", sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.sorted.plot, 
                ".pdf", sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.heatmap.sorted.plot, 
                ".pdf", sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
    }
    MSIG.HeatMapPlot.5(V = m2, row.names = gs.names2, col.labels = class.labels2, 
        col.classes = class.phen, col.names = sample.names2, 
        main = paste(title, " -- Heat Map (sorted)", sep = ""), 
        sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, 
        cmap.type = heatmap.cmap.type, char.rescale = 1)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    m2 <- m
    gs.names2 <- gs.names
    sample.names2 <- sample.names
    max.classes <- max(class.labels)
    for (k in 1:max.classes) {
        if (sum(class.labels == k) > 1) {
            m3 <- m2[, class.labels == k]
            sn <- sample.names2[class.labels == k]
            dist.matrix <- dist(t(m3))
            HC <- hclust(dist.matrix, method = "complete")
            m3 <- m3[, HC$order]
            sn <- sn[HC$order]
            m2[, class.labels == k] <- m3
            sample.names2[class.labels == k] <- sn
        }
    }
    height <- ifelse(k.proj > 50, 20, 0.5 * k.proj + 7)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.heatmap.sorted.2.plot
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.sorted.2.plot, 
                ".pdf", sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.sorted.2.plot, 
                ".pdf", sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.heatmap.sorted.2.plot, 
                ".pdf", sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 25)
        }
    }
    MSIG.HeatMapPlot.5(V = m2, row.names = gs.names2, col.labels = class.labels, 
        col.classes = class.phen, col.names = sample.names2, 
        main = paste(title, " -- Heat Map (sorted inside phenotype class)", 
            sep = ""), sub = " ", xlab = " ", ylab = " ", row.norm = heatmap.row.norm, 
        cmap.type = heatmap.cmap.type, char.rescale = 1)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
