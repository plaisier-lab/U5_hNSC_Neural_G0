\name{MSIG.Select.Top.Markers}
\alias{MSIG.Select.Top.Markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Select.Top.Markers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Select.Top.Markers(input.ds, input.cls, output.marker.report, output.marker.file, output.marker.gene.set.file = "", output.marker.plot, up.and.down.markers = F, topgs = 10, seed = 1234, non.interactive.run = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.marker.report}{
%%     ~~Describe \code{output.marker.report} here~~
}
  \item{output.marker.file}{
%%     ~~Describe \code{output.marker.file} here~~
}
  \item{output.marker.gene.set.file}{
%%     ~~Describe \code{output.marker.gene.set.file} here~~
}
  \item{output.marker.plot}{
%%     ~~Describe \code{output.marker.plot} here~~
}
  \item{up.and.down.markers}{
%%     ~~Describe \code{up.and.down.markers} here~~
}
  \item{topgs}{
%%     ~~Describe \code{topgs} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
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
function (input.ds, input.cls, output.marker.report, output.marker.file, 
    output.marker.gene.set.file = "", output.marker.plot, up.and.down.markers = F, 
    topgs = 10, seed = 1234, non.interactive.run = F) 
{
    print("Running MSIG.Select.Top.Markers...")
    if (output.marker.gene.set.file != "") {
        gs.file <- file(output.marker.gene.set.file, "w")
    }
    if (regexpr(pattern = ".gct", input.ds) == -1) {
        dataset <- GSEA.Res2Frame(filename = input.ds)
        gs.names <- row.names(dataset)
        sample.names <- names(dataset)
        m <- data.matrix(dataset)
    }
    else {
        dataset <- MSIG.Gct2Frame(filename = input.ds)
        gs.names <- dataset$row.names
        gs.descs <- dataset$descs
        sample.names <- dataset$names
        m <- data.matrix(dataset$ds)
    }
    m1 <- m
    dim(m)
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
    if (up.and.down.markers == F) {
        topgs <- ifelse(topgs > floor(Ng/length(class.phen)), 
            floor(Ng/length(class.phen)), topgs)
    }
    else {
        topgs <- ifelse(topgs > floor(Ng/(2 * length(class.phen))), 
            floor(Ng/(2 * length(class.phen))), topgs)
    }
    if (up.and.down.markers == F) {
        sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * 
            topgs, ncol = Ns)
        sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * 
            topgs, mode = "character")
        sample.molsig.sorted.s2n <- vector(length = length(class.phen) * 
            topgs, mode = "character")
        sample.molsig.sorted.class <- vector(length = length(class.phen) * 
            topgs, mode = "character")
    }
    else {
        sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * 
            2 * topgs, ncol = Ns)
        sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * 
            2 * topgs, mode = "character")
        sample.molsig.sorted.s2n <- vector(length = length(class.phen) * 
            2 * topgs, mode = "character")
        sample.molsig.sorted.class <- vector(length = length(class.phen) * 
            2 * topgs, mode = "character")
    }
    for (k in 1:length(class.phen)) {
        class.k.labels <- ifelse(class.labels == k, 0, 1)
        col.index <- order(class.k.labels, decreasing = F)
        class.k.labels <- class.k.labels[col.index]
        for (j in 1:Ng) {
            m1[j, ] <- m[j, col.index]
        }
        names(m1) <- sample.names
        set.seed(seed)
        O <- GSEA.GeneRanking(m1, class.k.labels, gene.labels, 
            1, permutation.type = 0, sigma.correction = "GeneCluster", 
            fraction = 1, replace = F, reverse.sign = F)
        order.matrix <- O$order.matrix
        obs.order.matrix <- O$obs.order.matrix
        correl.matrix <- O$s2n.matrix
        obs.correl.matrix <- O$obs.s2n.matrix
        rm(O)
        obs.s2n.orig <- apply(obs.correl.matrix, 1, median)
        obs.index <- order(obs.s2n.orig, decreasing = T)
        obs.s2n <- sort(obs.s2n.orig, decreasing = T)
        sample.molsig.sorted <- m[obs.index, ]
        gs.names.sorted <- gs.names[obs.index]
        if (up.and.down.markers == F) {
            start <- (k - 1) * topgs + 1
            end <- k * topgs
            sample.molsig.sorted.subset[start:end, ] <- sample.molsig.sorted[1:topgs, 
                ]
            sample.molsig.sorted.subset.gs[start:end] <- gs.names.sorted[1:topgs]
            sample.molsig.sorted.s2n[start:end] <- signif(obs.s2n[1:topgs], 
                digits = 3)
            sample.molsig.sorted.class[start:end] <- class.phen[k]
        }
        else {
            start <- (k - 1) * 2 * topgs + 1
            sample.molsig.sorted.subset[start:(start + topgs - 
                1), ] <- sample.molsig.sorted[1:topgs, ]
            sample.molsig.sorted.subset[(start + topgs):(start + 
                2 * topgs - 1), ] <- sample.molsig.sorted[seq(Ng, 
                Ng - topgs + 1, -1), ]
            sample.molsig.sorted.subset.gs[start:(start + topgs - 
                1)] <- gs.names.sorted[1:topgs]
            sample.molsig.sorted.subset.gs[(start + topgs):(start + 
                2 * topgs - 1)] <- gs.names.sorted[seq(Ng, Ng - 
                topgs + 1, -1)]
            sample.molsig.sorted.s2n[start:(start + topgs - 1)] <- signif(obs.s2n[1:topgs], 
                digits = 3)
            sample.molsig.sorted.s2n[(start + topgs):(start + 
                2 * topgs - 1)] <- signif(obs.s2n[seq(Ng, Ng - 
                topgs + 1, -1)], digits = 3)
            sample.molsig.sorted.class[start:(start + 2 * topgs - 
                1)] <- class.phen[k]
        }
        if (output.marker.gene.set.file != "") {
            if (up.and.down.markers == F) {
                gene.set <- paste(gs.names.sorted[1:topgs], sep = "\t")
                gene.set.name <- paste("Markers_of_", class.phen[k], 
                  sep = "")
                gene.set.desc <- paste("Top markers of phenotype: ", 
                  class.phen[k], sep = "")
                cat(gene.set.name, gene.set.desc, gene.set, "\n", 
                  file = gs.file, append = TRUE, sep = "\t")
            }
            else {
                gene.set <- paste(gs.names.sorted[1:topgs], sep = "\t")
                gene.set.name <- paste("UP_Markers_of_", class.phen[k], 
                  sep = "")
                gene.set.desc <- paste("Top UP markers of phenotype: ", 
                  class.phen[k], sep = "")
                cat(gene.set.name, gene.set.desc, gene.set, "\n", 
                  file = gs.file, append = TRUE, sep = "\t")
                gene.set <- paste(gs.names.sorted[seq(Ng, Ng - 
                  topgs + 1, -1)], sep = "\t")
                gene.set.name <- paste("DOWN_Markers_of_", class.phen[k], 
                  sep = "")
                gene.set.desc <- paste("Top DOWN markers of phenotype: ", 
                  class.phen[k], sep = "")
                cat(gene.set.name, gene.set.desc, gene.set, "\n", 
                  file = gs.file, append = TRUE, sep = "\t")
            }
        }
    }
    if (output.marker.gene.set.file != "") {
        close(gs.file)
    }
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            glob.filename <- output.marker.plot
            x11(height = 22, width = 40)
        }
        else if (.Platform$OS.type == "unix") {
            glob.filename <- paste(output.marker.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = glob.filename, height = 22, width = 40)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            glob.filename <- paste(output.marker.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = glob.filename, height = 22, width = 40)
        }
        else if (.Platform$OS.type == "windows") {
            glob.filename <- paste(output.marker.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = glob.filename, height = 22, width = 40)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 
        1), heights = 1, respect = FALSE)
    c1 <- c("grey3", "cadetblue", "darkgreen", "chartreuse2", 
        "red1", "darkred", "orange", "blue2", "lightblue", "pink", 
        "coral")
    MSIG.HeatMapPlot.3(V = sample.molsig.sorted.subset, row.names = sample.molsig.sorted.subset.gs, 
        col.labels = class.labels, col.classes = class.phen, 
        phen.cmap = c1[1:length(class.phen)], col.names = sample.names, 
        main = "Top Markers -- Heat Map", xlab = " ", ylab = " ", 
        sub = " ", row.norm = T, cmap.type = 4)
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(22, n.phen)
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.1, pt.cex = 1.1)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = glob.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    print("Saving markers dataset...")
    n.markers <- length(sample.molsig.sorted.subset.gs)
    report <- data.frame(cbind(sample.molsig.sorted.subset.gs, 
        sample.molsig.sorted.s2n, sample.molsig.sorted.class))
    names(report) <- c("Name", "S2N", "class")
    row.names(report) <- seq(1, n.markers)
    write.table(report, file = output.marker.report, quote = F, 
        sep = "\t")
    V <- data.frame(sample.molsig.sorted.subset)
    names(V) <- sample.names
    row.names(V) <- sample.molsig.sorted.subset.gs
    write.gct(gct.data.frame = V, descs = seq(1, n.markers), 
        filename = output.marker.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
