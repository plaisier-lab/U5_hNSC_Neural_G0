\name{MSIG.File.to.HeatMap.2}
\alias{MSIG.File.to.HeatMap.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.File.to.HeatMap.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.File.to.HeatMap.2(input.ds, input.cls = "", output.heatmap.plot, output.heatmap.sorted.plot, output.heatmap.sorted.2.plot, title = "", non.interactive.run = F, heatmap.row.norm = F, phen.cmap = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"), heatmap.cmap.type = 1, char.rescale = 1)
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
  \item{output.heatmap.sorted.2.plot}{
%%     ~~Describe \code{output.heatmap.sorted.2.plot} here~~
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
    output.heatmap.sorted.2.plot, title = "", non.interactive.run = F, 
    heatmap.row.norm = F, phen.cmap = c("grey3", "cadetblue", 
        "darkgreen", "chartreuse2", "red1", "darkred", "orange", 
        "blue2", "lightblue", "pink", "coral"), heatmap.cmap.type = 1, 
    char.rescale = 1) 
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
    MSIG.HeatMapPlot.4(V = m, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, 
        main = title, xlab = " ", ylab = " ", sub = "heat map ", 
        row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, 
        char.rescale = char.rescale)
    savePlot(filename = output.heatmap.plot, type = "jpeg", device = dev.cur())
    dist.matrix <- dist(t(m))
    HC <- hclust(dist.matrix, method = "complete")
    m1 <- m[, HC$order]
    sample.names1 <- sample.names[HC$order]
    class.labels1 <- class.labels[HC$order]
    dist.matrix <- dist(m1)
    HC <- hclust(dist.matrix, method = "complete")
    m1 <- m1[HC$order, ]
    gs.names1 <- gs.names[HC$order]
    MSIG.HeatMapPlot.4(V = m1, row.names = gs.names1, col.labels = class.labels1, 
        col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names1, 
        main = title, xlab = " ", ylab = " ", sub = "sorted heat map ", 
        row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, 
        char.rescale = char.rescale)
    savePlot(filename = output.heatmap.sorted.plot, type = "jpeg", 
        device = dev.cur())
    dist.matrix <- dist(m)
    HC <- hclust(dist.matrix, method = "complete")
    m2 <- m[HC$order, ]
    gs.names2 <- gs.names[HC$order]
    sample.names2 <- sample.names
    max.classes <- max(class.labels)
    for (k in 1:max.classes) {
        m3 <- m2[, class.labels == k]
        sn <- sample.names2[class.labels == k]
        dist.matrix <- dist(t(m3))
        HC <- hclust(dist.matrix, method = "complete")
        m3 <- m3[, HC$order]
        sn <- sn[HC$order]
        m2[, class.labels == k] <- m3
        sample.names2[class.labels == k] <- sn
    }
    MSIG.HeatMapPlot.4(V = m2, row.names = gs.names2, col.labels = class.labels, 
        col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names2, 
        main = title, xlab = " ", ylab = " ", sub = "sorted heat map (inside phenotype)", 
        row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type, 
        char.rescale = char.rescale)
    savePlot(filename = output.heatmap.sorted.2.plot, type = "jpeg", 
        device = dev.cur())
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
