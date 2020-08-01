\name{MSIG.Select.UP.DN.Markers}
\alias{MSIG.Select.UP.DN.Markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Select.UP.DN.Markers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Select.UP.DN.Markers(input.ds, input.cls, output.marker.file, output.marker.plot, output.marker.gct, num.of.markers = 10, disc.metric = "S2N.DIFF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.marker.file}{
%%     ~~Describe \code{output.marker.file} here~~
}
  \item{output.marker.plot}{
%%     ~~Describe \code{output.marker.plot} here~~
}
  \item{output.marker.gct}{
%%     ~~Describe \code{output.marker.gct} here~~
}
  \item{num.of.markers}{
%%     ~~Describe \code{num.of.markers} here~~
}
  \item{disc.metric}{
%%     ~~Describe \code{disc.metric} here~~
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
function (input.ds, input.cls, output.marker.file, output.marker.plot, 
    output.marker.gct, num.of.markers = 10, disc.metric = "S2N.DIFF") 
{
    print("Running MSIG.Select.UP.DN.Markers...")
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
    num.of.markers <- ifelse(num.of.markers > floor(Ng/(2 * length(class.phen))), 
        floor(Ng/(2 * length(class.phen))), num.of.markers)
    sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * 
        2 * num.of.markers, ncol = Ns)
    sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * 
        2 * num.of.markers, mode = "character")
    sample.molsig.sorted.s2n <- vector(length = length(class.phen) * 
        2 * num.of.markers, mode = "character")
    sample.molsig.sorted.class <- vector(length = length(class.phen) * 
        2 * num.of.markers, mode = "character")
    class.k.labels <- class.labels
    col.index <- order(class.k.labels, decreasing = F)
    class.k.labels <- class.k.labels[col.index]
    for (j in 1:Ng) {
        m1[j, ] <- m[j, col.index]
    }
    names(m1) <- sample.names
    print("Executing marker selection...")
    if (disc.metric == "S2N") {
        obs.s2n <- Gene.ranking(m1, class.k.labels, method = "S2N")
        fold.changes <- fold.changes(m1, class.k.labels, method = "MEAN.DIFF", 
            thres = 1)
    }
    else if (disc.metric == "RS2N") {
        obs.s2n <- Gene.ranking(m1, class.k.labels, method = "RS2N")
        fold.changes <- fold.changes(m1, class.k.labels, method = "MEDIAN.DIFF", 
            thres = 1)
    }
    else if (disc.metric == "MEAN.DIFF") {
        obs.s2n <- Gene.ranking(m1, class.k.labels, method = "MEAN.DIFF")
        fold.changes <- fold.changes(m1, class.k.labels, method = "MEAN.DIFF", 
            thres = 1)
    }
    else if (disc.metric == "MEDIAN.DIFF") {
        obs.s2n <- Gene.ranking(m1, class.k.labels, method = "MEDIAN.DIFF")
        fold.changes <- fold.changes(m1, class.k.labels, method = "MEDIAN.DIFF", 
            thres = 1)
    }
    obs.index <- order(obs.s2n, decreasing = T)
    obs.s2n <- sort(obs.s2n, decreasing = T)
    sample.molsig.sorted <- m[obs.index, ]
    gs.names.sorted <- gs.names[obs.index]
    msig.up <- m1[obs.index[1:num.of.markers], ]
    msig.up.genes <- gs.names[obs.index[1:num.of.markers]]
    msig.up.genes <- paste(msig.up.genes, signif(fold.changes[obs.index[1:num.of.markers]], 
        digits = 3), sep = "_")
    msig.dn <- m1[obs.index[seq(Ng, Ng - num.of.markers + 1, 
        -1)], ]
    msig.dn.genes <- gs.names[obs.index[seq(Ng, Ng - num.of.markers + 
        1, -1)]]
    msig.dn.genes <- paste(msig.dn.genes, signif(1/fold.changes[obs.index[seq(Ng, 
        Ng - num.of.markers + 1, -1)]], digits = 3), sep = "_")
    descs.up <- paste("UP", seq(1, num.of.markers), sep = "_")
    descs.dn <- paste("DN", seq(1, num.of.markers), sep = "_")
    descs <- c(descs.up, descs.dn)
    s2n.up <- signif(obs.s2n[1:num.of.markers], digits = 4)
    s2n.dn <- signif(obs.s2n[seq(Ng, Ng - num.of.markers + 1, 
        -1)], digits = 4)
    s2n <- c(s2n.up, s2n.dn)
    markers <- data.frame(cbind(descs, c(msig.up.genes, msig.dn.genes), 
        s2n))
    names(markers) <- c("UP/DN Rank", "Gene", disc.metric)
    descs.up <- paste(descs.up, s2n.up, sep = "_")
    descs.dn <- paste(descs.dn, s2n.dn, sep = "_")
    print(markers)
    write.table(markers, file = output.marker.file, quote = F, 
        row.names = F, sep = "\t")
    c1 <- c("red", "blue")
    V <- data.frame(rbind(msig.up, msig.dn))
    names(V) <- sample.names
    row.names(V) <- c(msig.up.genes, msig.dn.genes)
    write.gct(gct.data.frame = V, descs = descs, filename = output.marker.gct)
    msig.up.genes <- paste(msig.up.genes, descs.up, sep = "     ")
    glob.filename <- paste(output.marker.plot, ".UP", sep = "")
    x11(height = 9, width = 12)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 
        1), heights = 1, respect = FALSE)
    MSIG.HeatMapPlot.3(V = msig.up, row.names = msig.up.genes, 
        col.labels = class.labels, col.classes = class.phen, 
        phen.cmap = c1[1:length(class.phen)], col.names = sample.names, 
        main = "Top UP Markers", xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1.25)
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(22, n.phen)
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.5, pt.cex = 1.5)
    savePlot(filename = glob.filename, type = "jpeg", device = dev.cur())
    msig.dn.genes <- paste(msig.dn.genes, descs.dn, sep = "     ")
    glob.filename <- paste(output.marker.plot, ".DN", sep = "")
    x11(height = 9, width = 12)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 
        1), heights = 1, respect = FALSE)
    MSIG.HeatMapPlot.3(V = msig.dn, row.names = msig.dn.genes, 
        col.labels = class.labels, col.classes = class.phen, 
        phen.cmap = c1[1:length(class.phen)], col.names = sample.names, 
        main = "Top DN Markers", xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1.25)
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(22, n.phen)
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.5, pt.cex = 1.5)
    savePlot(filename = glob.filename, type = "jpeg", device = dev.cur())
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
