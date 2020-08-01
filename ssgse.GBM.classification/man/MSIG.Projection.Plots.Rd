\name{MSIG.Projection.Plots}
\alias{MSIG.Projection.Plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Projection.Plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Projection.Plots(input.ds, input.cls = "", output.2D.sammon.file, output.2D.sammon.plot, output.3D.sammon.file, output.3D.sammon.plot, output.heatmap.plot, output.hclust.plot, title = "", niter = 1000, seed = 1234, non.interactive.run = F, col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.2D.sammon.file}{
%%     ~~Describe \code{output.2D.sammon.file} here~~
}
  \item{output.2D.sammon.plot}{
%%     ~~Describe \code{output.2D.sammon.plot} here~~
}
  \item{output.3D.sammon.file}{
%%     ~~Describe \code{output.3D.sammon.file} here~~
}
  \item{output.3D.sammon.plot}{
%%     ~~Describe \code{output.3D.sammon.plot} here~~
}
  \item{output.heatmap.plot}{
%%     ~~Describe \code{output.heatmap.plot} here~~
}
  \item{output.hclust.plot}{
%%     ~~Describe \code{output.hclust.plot} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
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
function (input.ds, input.cls = "", output.2D.sammon.file, output.2D.sammon.plot, 
    output.3D.sammon.file, output.3D.sammon.plot, output.heatmap.plot, 
    output.hclust.plot, title = "", niter = 1000, seed = 1234, 
    non.interactive.run = F, col = c("grey3", "cadetblue", "darkgreen", 
        "chartreuse2", "red1", "darkred", "orange", "blue2", 
        "lightblue", "pink", "coral")) 
{
    print("Running MSIG.Projection.Plots...")
    library("scatterplot3d")
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
    sammon.map <- sammon(d = dist(t(m)), k = 2, niter = niter, 
        trace = TRUE, magic = 0.2, tol = 1e-06)
    S1 <- sammon.map$points[, 1]
    S2 <- sammon.map$points[, 2]
    S1 <- (S1 - min(S1))/(max(S1) - min(S1))
    S2 <- (S2 - min(S2))/(max(S2) - min(S2))
    V <- data.frame(cbind(S1, S2))
    names(V) <- c("S1", "S2")
    row.names(V) <- sample.names
    write.gct(gct.data.frame = V, filename = output.2D.sammon.file)
    c0 <- col
    c1 <- colors()[match(c0, colors())]
    color <- c1[class.labels]
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.2D.sammon.plot
            x11(height = 8, width = 14)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.2D.sammon.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 8, width = 14)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.2D.sammon.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 8, width = 14)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.2D.sammon.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 8, width = 14)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 
        1), heights = 1, respect = FALSE)
    num.samples <- length(S1)
    plot(S1, S2, xlim = c(0, 1), ylim = c(0, 1), type = "n", 
        main = paste(title, " -- 2D Sammon Map", sep = ""), sub = input.ds)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            color.code <- c1[class.labels[j]]
        }
        points(S1[j], S2[j], pch = 21, type = "p", cex = 2.5, 
            bg = color.code, col = "black")
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(21, n.phen)
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.5, pt.cex = 2)
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
    sammon.map <- sammon(d = dist(t(m)), k = 3, niter = niter, 
        trace = TRUE, magic = 0.2, tol = 1e-06)
    S1 <- sammon.map$points[, 3]
    S2 <- sammon.map$points[, 2]
    S3 <- sammon.map$points[, 1]
    S1 <- (S1 - min(S1))/(max(S1) - min(S1))
    S2 <- (S2 - min(S2))/(max(S2) - min(S2))
    S3 <- (S3 - min(S3))/(max(S3) - min(S3))
    V <- data.frame(cbind(S1, S2, S3))
    names(V) <- c("S1", "S2", "S3")
    row.names(V) <- sample.names
    write.gct(gct.data.frame = V, filename = output.3D.sammon.file)
    color <- c1[class.labels]
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.3D.sammon.plot
            x11(height = 14, width = 22)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.sammon.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 22)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.3D.sammon.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 22)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.2D.sammon.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 22)
        }
    }
    nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(3, 
        3, 1), heights = 1, respect = FALSE)
    x <- scatterplot3d(S1, S2, S3, xlab = "F1", ylab = "F2", 
        zlab = "F3", type = "p", color, angle = 45, pch = 20, 
        main = paste(title, " -- 3D Sammon Map", sep = ""), sub = " ", 
        cex.symbols = 1)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            color.code <- c1[class.labels[j]]
        }
        cex <- 3 * (max(S2) - S2[j])/(max(S2) - min(S2)) + 1
        x$points3d(S1[j], S2[j], S3[j], col = "black", pch = 21, 
            bg = color.code, cex = cex)
    }
    x <- scatterplot3d(S2, S1, S3, xlab = "F2", ylab = "F1", 
        zlab = "F3", type = "p", color, angle = 45, pch = 20, 
        main = paste(title, " -- 3D Sammon Map", sep = ""), sub = input.ds, 
        cex.symbols = 1)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            color.code <- c1[class.labels[j]]
        }
        cex <- 3 * (max(S1) - S1[j])/(max(S1) - min(S1)) + 1
        x$points3d(S2[j], S1[j], S3[j], col = "black", pch = 21, 
            bg = color.code, cex = cex)
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(21, n.phen)
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 2, pt.cex = 2.5)
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
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.heatmap.plot
            x11(height = 14, width = 15)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 15)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 15)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 15)
        }
    }
    MSIG.HeatMapPlot(V = m, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, col.names = sample.names, main = paste(title, 
            " -- Heat Map", sep = ""), sub = " ", xlab = " ", 
        ylab = " ")
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
            x11(height = 15, width = 15)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.hclust.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 15, width = 15)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.hclust.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 15, width = 15)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.hclust.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 15, width = 15)
        }
    }
    plot(HC, xlab = "samples", cex = 0.75, labels = class.phen[class.labels], 
        col = "blue", main = paste(title, " -- Hierarchical Clustering", 
            sep = ""))
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
