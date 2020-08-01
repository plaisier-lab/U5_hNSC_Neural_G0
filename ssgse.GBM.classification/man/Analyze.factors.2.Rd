\name{Analyze.factors.2}
\alias{Analyze.factors.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyze.factors.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Analyze.factors.2(input.ds, input.cls, output.dir, identifier, model.size, seed = 123, nchar.phen = 2, use.biplot = TRUE, non.interactive.run = FALSE, heatmap.row.norm = FALSE, heatmap.cmap.type = 3, use.feature.names = FALSE, high.conf.thres = 0.5, col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"), symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25), symbol.scaling = 1, levels = NULL, nlevels = 20, kernel = "radial", cost = 1, gamma = 0.05, theta = 0, model.set.refinement = T, produce.contours = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
}
  \item{identifier}{
%%     ~~Describe \code{identifier} here~~
}
  \item{model.size}{
%%     ~~Describe \code{model.size} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{nchar.phen}{
%%     ~~Describe \code{nchar.phen} here~~
}
  \item{use.biplot}{
%%     ~~Describe \code{use.biplot} here~~
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
  \item{use.feature.names}{
%%     ~~Describe \code{use.feature.names} here~~
}
  \item{high.conf.thres}{
%%     ~~Describe \code{high.conf.thres} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{symbs}{
%%     ~~Describe \code{symbs} here~~
}
  \item{symbol.scaling}{
%%     ~~Describe \code{symbol.scaling} here~~
}
  \item{levels}{
%%     ~~Describe \code{levels} here~~
}
  \item{nlevels}{
%%     ~~Describe \code{nlevels} here~~
}
  \item{kernel}{
%%     ~~Describe \code{kernel} here~~
}
  \item{cost}{
%%     ~~Describe \code{cost} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{model.set.refinement}{
%%     ~~Describe \code{model.set.refinement} here~~
}
  \item{produce.contours}{
%%     ~~Describe \code{produce.contours} here~~
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
function (input.ds, input.cls, output.dir, identifier, model.size, 
    seed = 123, nchar.phen = 2, use.biplot = TRUE, non.interactive.run = FALSE, 
    heatmap.row.norm = FALSE, heatmap.cmap.type = 3, use.feature.names = FALSE, 
    high.conf.thres = 0.5, col = c("green", "blue", "pink", "red", 
        "orange", "red4", "steelblue2", "violet"), symbs = c(22, 
        21, 20, 23, 24, 25, 21, 20, 23, 24, 25), symbol.scaling = 1, 
    levels = NULL, nlevels = 20, kernel = "radial", cost = 1, 
    gamma = 0.05, theta = 0, model.set.refinement = T, produce.contours = T) 
{
    print(c("Running MSIG.Analyze.factors.2... on: ", input.ds, 
        " ", input.cls))
    O <- MSIG.Evaluate.Projection.2(input.ds = input.ds, input.cls = input.cls, 
        model.set = seq(1, model.size), prediction.results.file = paste(output.dir, 
            identifier, ".pred.txt", sep = ""), prediction.matrix.file = paste(output.dir, 
            identifier, ".pred.gct", sep = ""), train.pred.plot = paste(output.dir, 
            identifier, ".train.pred", sep = ""), test.pred.plot = paste(output.dir, 
            identifier, ".test.pred", sep = ""), pred.2D.plot = paste(output.dir, 
            identifier, ".2D.pred", sep = ""), col = col, use.feature.names = use.feature.names, 
        nchar.phen = nchar.phen, high.conf.thres = high.conf.thres, 
        symbol.scaling = symbol.scaling, symbs = symbs, levels = levels, 
        nlevels = nlevels, kernel = kernel, cost = cost, gamma = gamma, 
        produce.contours = produce.contours)
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m.ds <- data.matrix(dataset$ds)
    N.ds <- length(m.ds[, 1])
    M.ds <- length(m.ds[1, ])
    ds.names <- dataset$row.names
    ds.descs <- dataset$descs
    ds.sample.names <- dataset$names
    CLS <- MSIG.ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
    class.labels <- match(class.list, class.phen)
    class.phen <- unique(class.list)
    dist.matrix <- dist(t(m.ds))
    HC <- hclust(dist.matrix, method = "complete")
    x11(height = 24, width = 30)
    HC$labels <- class.list
    dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2, 
        edgePar = list(col = c("blue", "green"), lty = c(1, 1), 
            lwd = c(2, 2), t.col = c(1, 1)))
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
    plot(dL, cex = 1, edge.root = T, main = " Hierarchical Clustering (original data)", 
        xlab = "samples")
    plot.filename <- paste(output.dir, identifier, ".htree", 
        sep = "")
    savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m.ds <- data.matrix(dataset$ds)
    N.ds <- length(m.ds[, 1])
    M.ds <- length(m.ds[1, ])
    ds.names <- dataset$row.names
    ds.descs <- dataset$descs
    ds.sample.names <- dataset$names
    dist.matrix <- dist(t(m.ds))
    HC <- hclust(dist.matrix, method = "complete")
    x11(height = 24, width = 30)
    HC$labels <- class.list
    dhc <- as.dendrogram(HC, hang = 0.01, edge.root = T, dLeaf = 2, 
        edgePar = list(col = c("blue", "green"), lty = c(1, 1), 
            lwd = c(2, 2), t.col = c(1, 1)))
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
    plot(dL, cex = 1, edge.root = T, main = " Hierarchical Clustering (projected data)", 
        xlab = "samples")
    plot.filename <- paste(output.dir, identifier, ".H.htree", 
        sep = "")
    savePlot(filename = plot.filename, type = "jpeg", device = dev.cur())
    membership <- vector(length = M.ds, mode = "numeric")
    for (j in 1:M.ds) {
        membership[j] <- order(m.ds[, j], decreasing = T)
    }
    mem.order <- order(membership, decreasing = F)
    membership.sorted <- membership[mem.order]
    ds.sample.names <- paste(class.list, ds.sample.names, sep = "_")
    ds.sample.names.sorted <- ds.sample.names[mem.order]
    class.list.sorted <- class.list[mem.order]
    mem.table <- data.frame(cbind(class.list, ds.sample.names, 
        membership, rep(" ", M.ds), class.list.sorted, ds.sample.names.sorted, 
        membership.sorted))
    row.names(mem.table) <- seq(1, M.ds)
    names(mem.table) <- c("Phen", "Sample Names", "Membership", 
        " ", "Phen Sorted", "Sample Names Sorted", "Membership Sorted")
    mem.filename <- paste(output.dir, identifier, ".H.mem.txt", 
        sep = "")
    write.table(file = mem.filename, mem.table, quote = F, sep = "\t")
    table(class.list.sorted, membership.sorted)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
