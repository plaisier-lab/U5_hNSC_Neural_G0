\name{MSIG.File.to.HeatMap}
\alias{MSIG.File.to.HeatMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.File.to.HeatMap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.File.to.HeatMap(input.ds, input.cls = "", output.heatmap.plot, output.heatmap.sorted.plot, title = "", non.interactive.run = F, heatmap.row.norm = F, phen.cmap = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"), heatmap.cmap.type = 1, char.rescale = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.heatmap.plot}{
%%     ~~Describe \code{output.heatmap.plot} here~~
}
  \item{output.heatmap.sorted.plot}{
%%     ~~Describe \code{output.heatmap.sorted.plot} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{heatmap.row.norm}{
%%     ~~Describe \code{heatmap.row.norm} here~~
}
  \item{phen.cmap}{
%%     ~~Describe \code{phen.cmap} here~~
}
  \item{heatmap.cmap.type}{
%%     ~~Describe \code{heatmap.cmap.type} here~~
}
  \item{char.rescale}{
%%     ~~Describe \code{char.rescale} here~~
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
function (input.ds, input.cls = "", output.heatmap.plot, output.heatmap.sorted.plot, 
    title = "", non.interactive.run = F, heatmap.row.norm = F, 
    phen.cmap = c("grey3", "cadetblue", "darkgreen", "chartreuse2", 
        "red1", "darkred", "orange", "blue2", "lightblue", "pink", 
        "coral"), heatmap.cmap.type = 1, char.rescale = 0.8) 
{
    print("Running MSIG.File.to.HeatMap...")
    print(heatmap.row.norm)
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
    height <- ifelse(k.proj >= 25, 25, k.proj * 0.8 + 5)
    x11(width = 30, height = height)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(4, 
        1), respect = FALSE)
    MSIG.HeatMapPlot.3(V = m, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, 
        main = paste(title, " -- Heat Map", sep = ""), xlab = " ", 
        ylab = " ", sub = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, 
        char.rescale = char.rescale)
    leg.txt <- class.phen
    p.vec <- rep(21, 21)
    c.vec <- phen.cmap
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1, pt.cex = 1.5)
    savePlot(filename = output.heatmap.plot, type = "jpeg", device = dev.cur())
    height <- ifelse(k.proj >= 25, 25, k.proj * 0.8 + 5)
    x11(width = 30, height = height)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(4, 
        1), respect = FALSE)
    dist.matrix <- dist(t(m))
    HC <- hclust(dist.matrix, method = "complete")
    m <- m[, HC$order]
    sample.names <- sample.names[HC$order]
    class.labels <- class.labels[HC$order]
    dist.matrix <- dist(m)
    HC <- hclust(dist.matrix, method = "complete")
    m <- m[HC$order, ]
    gs.names <- gs.names[HC$order]
    MSIG.HeatMapPlot.3(V = m, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, 
        main = paste(title, " -- Heat Map (sorted)", sep = ""), 
        xlab = " ", ylab = " ", sub = " ", row.norm = heatmap.row.norm, 
        cmap.type = heatmap.cmap.type, char.rescale = char.rescale)
    leg.txt <- class.phen
    p.vec <- rep(21, 21)
    c.vec <- phen.cmap
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1, pt.cex = 1.5)
    savePlot(filename = output.heatmap.sorted.plot, type = "jpeg", 
        device = dev.cur())
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
