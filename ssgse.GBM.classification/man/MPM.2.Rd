\name{MPM.2}
\alias{MPM.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MPM.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MPM.2(model.dataset.table, test.datasets.table, identifier, k.proj, alg = "NMF.div", niter = 1500, seed = 123, nchar.phen = 2, postprojnorm = TRUE, use.biplot = TRUE, non.interactive.run = FALSE, heatmap.row.norm = FALSE, heatmap.cmap.type = 3, use.feature.names = FALSE, high.conf.thres = 0.5, output.dir, col = c("green", "blue", "pink", "red", "orange", "red4", "steelblue2", "violet"), symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25), symbol.scaling = 1, levels = NULL, nlevels = 20, kernel = "radial", cost = 1, gamma = 0.05, theta = 0, model.set.refinement = T, produce.contours = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model.dataset.table}{
%%     ~~Describe \code{model.dataset.table} here~~
}
  \item{test.datasets.table}{
%%     ~~Describe \code{test.datasets.table} here~~
}
  \item{identifier}{
%%     ~~Describe \code{identifier} here~~
}
  \item{k.proj}{
%%     ~~Describe \code{k.proj} here~~
}
  \item{alg}{
%%     ~~Describe \code{alg} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{nchar.phen}{
%%     ~~Describe \code{nchar.phen} here~~
}
  \item{postprojnorm}{
%%     ~~Describe \code{postprojnorm} here~~
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
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
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
function (model.dataset.table, test.datasets.table, identifier, 
    k.proj, alg = "NMF.div", niter = 1500, seed = 123, nchar.phen = 2, 
    postprojnorm = TRUE, use.biplot = TRUE, non.interactive.run = FALSE, 
    heatmap.row.norm = FALSE, heatmap.cmap.type = 3, use.feature.names = FALSE, 
    high.conf.thres = 0.5, output.dir, col = c("green", "blue", 
        "pink", "red", "orange", "red4", "steelblue2", "violet"), 
    symbs = c(22, 21, 20, 23, 24, 25, 21, 20, 23, 24, 25), symbol.scaling = 1, 
    levels = NULL, nlevels = 20, kernel = "radial", cost = 1, 
    gamma = 0.05, theta = 0, model.set.refinement = T, produce.contours = T) 
{
    print(c(model.dataset.table))
    print(c(test.datasets.table))
    set.seed(seed = seed, kind = NULL)
    O <- MSIG.Subset.Dataset(input.ds = model.dataset.table$gct.file, 
        input.cls = model.dataset.table$cls.file, column.subset = model.dataset.table$column.subset, 
        column.sel.type = model.dataset.table$column.sel.type, 
        row.subset = "ALL", output.ds = paste(output.dir, "temp1.gct", 
            sep = ""), output.cls = paste(output.dir, identifier, 
            ".model_set.1.cls", sep = ""))
    O <- MSIG.Preprocess.Dataset(input.ds = paste(output.dir, 
        "temp1.gct", sep = ""), output.ds = paste(output.dir, 
        identifier, ".model_set.1.gct", sep = ""), thres = model.dataset.table$thres, 
        ceil = model.dataset.table$ceil, fold = model.dataset.table$fold, 
        delta = model.dataset.table$delta, normalization = model.dataset.table$norm)
    O <- MSIG.Extract.Factors(input.ds = paste(output.dir, identifier, 
        ".model_set.1.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".model_set.1.cls", sep = ""), output.W.file = paste(output.dir, 
        identifier, ".model.W.gct", sep = ""), output.H.file = paste(output.dir, 
        identifier, ".model.H.gct", sep = ""), k.proj = k.proj, 
        alg = alg, niter = niter, seed = seed, theta = theta, 
        sort.factors = T)
    O <- MSIG.Factors.Project(input.ds = paste(output.dir, identifier, 
        ".model_set.1.gct", sep = ""), factors.ds = paste(output.dir, 
        identifier, ".model.W.gct", sep = ""), postprojnorm = postprojnorm, 
        output.file = paste(output.dir, identifier, ".model_dataset.H.gct", 
            sep = ""))
    CLS <- MSIG.ReadClsFile(file = paste(output.dir, identifier, 
        ".model_set.1.cls", sep = ""))
    model.size <- length(CLS$class.v)
    O <- MSIG.Evaluate.Projection.2(input.ds = paste(output.dir, 
        identifier, ".model_dataset.H.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".model_set.1.cls", sep = ""), model.set = seq(1, 
        model.size), prediction.results.file = paste(output.dir, 
        identifier, ".prelim.pred.txt", sep = ""), prediction.matrix.file = paste(output.dir, 
        identifier, ".prelim.pred.gct", sep = ""), train.pred.plot = paste(output.dir, 
        identifier, ".prelim.train.pred", sep = ""), test.pred.plot = paste(output.dir, 
        identifier, ".prelim.test.pred", sep = ""), pred.2D.plot = paste(output.dir, 
        identifier, ".prelim.2D.pred", sep = ""), col = col, 
        use.feature.names = use.feature.names, nchar.phen = nchar.phen, 
        high.conf.thres = high.conf.thres, symbs = symbs, symbol.scaling = symbol.scaling, 
        levels = levels, nlevels = nlevels, kernel = kernel, 
        cost = cost, gamma = gamma, produce.contours = F)
    input.txt <- paste(output.dir, identifier, ".prelim.pred.txt", 
        sep = "")
    pred.table <- read.table(file = input.txt, skip = 2, nrow = model.size, 
        sep = "\t", header = T, comment.char = "", as.is = T)
    conf.list <- pred.table["Conf..H.L."]
    actual.list <- pred.table["Actual"]
    predicted <- pred.table["Predicted"]
    sample.select <- ((conf.list == " H ") & (actual.list == 
        predicted))
    if (model.set.refinement == T) {
        high.conf.set <- seq(1, model.size)[sample.select]
    }
    else {
        high.conf.set <- seq(1, model.size)
    }
    print(c("pred table from file=", pred.table))
    print(c("model sizet=", model.size))
    print(c("high.conf.set=", high.conf.set))
    print(c("sample.select", sample.select))
    print(paste("Original: ", model.size, " (samples); New: ", 
        length(high.conf.set), " (samples); diff: ", model.size - 
            length(high.conf.set), sep = " "))
    O <- MSIG.Subset.Dataset(input.ds = paste(output.dir, identifier, 
        ".model_set.1.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".model_set.1.cls", sep = ""), column.subset = high.conf.set, 
        row.subset = "ALL", output.ds = paste(output.dir, identifier, 
            ".model_set.2.gct", sep = ""), output.cls = paste(output.dir, 
            identifier, ".model_set.2.cls", sep = ""))
    O <- MSIG.Subset.Dataset(input.ds = model.dataset.table$gct.file, 
        input.cls = model.dataset.table$cls.file, column.subset = high.conf.set, 
        row.subset = "ALL", output.ds = paste(output.dir, identifier, 
            ".model_set.0.gct", sep = ""), output.cls = paste(output.dir, 
            identifier, ".model_set.0.cls", sep = ""))
    O <- MSIG.Extract.Factors(input.ds = paste(output.dir, identifier, 
        ".model_set.2.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".model_set.2.cls", sep = ""), output.W.file = paste(output.dir, 
        identifier, ".model.W.gct", sep = ""), output.H.file = paste(output.dir, 
        identifier, ".model.H.gct", sep = ""), k.proj = k.proj, 
        alg = alg, niter = niter, seed = seed, theta = theta, 
        sort.factors = T)
    O <- MSIG.Factors.Project(input.ds = paste(output.dir, identifier, 
        ".model_set.2.gct", sep = ""), factors.ds = paste(output.dir, 
        identifier, ".model.W.gct", sep = ""), postprojnorm = postprojnorm, 
        output.file = paste(output.dir, identifier, ".model_dataset.H.gct", 
            sep = ""))
    O <- MSIG.Subset.Dataset(input.ds = paste(output.dir, identifier, 
        ".model_dataset.H.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".model_set.2.cls", sep = ""), column.subset = "ALL", 
        column.sel.type = "samples", row.subset = "ALL", output.ds = paste(output.dir, 
            identifier, ".all.H.gct", sep = ""), output.cls = paste(output.dir, 
            identifier, ".all.H.cls", sep = ""))
    O <- MSIG.Subset.Dataset(input.ds = paste(output.dir, identifier, 
        ".model_set.2.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".model_set.2.cls", sep = ""), column.subset = "ALL", 
        column.sel.type = "samples", row.subset = "ALL", output.ds = paste(output.dir, 
            identifier, ".all.gct", sep = ""), output.cls = paste(output.dir, 
            identifier, ".all.cls", sep = ""))
    if (!is.null(test.datasets.table)) {
        max.files <- length(test.datasets.table)
        for (ds in 1:max.files) {
            print(c("Processing test file: ", test.datasets.table[[ds]]$gct.file))
            O <- MSIG.Subset.Dataset(input.ds = test.datasets.table[[ds]]$gct.file, 
                input.cls = test.datasets.table[[ds]]$cls.file, 
                column.subset = test.datasets.table[[ds]]$column.subset, 
                column.sel.type = test.datasets.table[[ds]]$column.sel.type, 
                row.subset = "ALL", output.ds = paste(output.dir, 
                  "temp1.gct", sep = ""), output.cls = paste(output.dir, 
                  "temp1.cls", sep = ""))
            O <- MSIG.Preprocess.Dataset(input.ds = paste(output.dir, 
                "temp1.gct", sep = ""), output.ds = paste(output.dir, 
                "temp2.gct", sep = ""), thres = test.datasets.table[[ds]]$thres, 
                ceil = test.datasets.table[[ds]]$ceil, normalization = "NULL")
            O <- MSIG.Match.and.Select(input1.ds = paste(output.dir, 
                identifier, ".model_set.2.gct", sep = ""), input2.ds = paste(output.dir, 
                "temp2.gct", sep = ""), output.ds = paste(output.dir, 
                "temp3.gct", sep = ""))
            O <- MSIG.Preprocess.Dataset(input.ds = paste(output.dir, 
                "temp3.gct", sep = ""), output.ds = paste(output.dir, 
                "temp4.gct", sep = ""), normalization = test.datasets.table[[ds]]$norm)
            O <- MSIG.Factors.Project(input.ds = paste(output.dir, 
                "temp4.gct", sep = ""), factors.ds = paste(output.dir, 
                identifier, ".model.W.gct", sep = ""), postprojnorm = postprojnorm, 
                output.file = paste(output.dir, "temp5.gct", 
                  sep = ""))
            O <- MSIG.Match.and.Merge(input1.ds = paste(output.dir, 
                identifier, ".all.H.gct", sep = ""), input1.cls = paste(output.dir, 
                identifier, ".all.H.cls", sep = ""), input2.ds = paste(output.dir, 
                "temp5.gct", sep = ""), input2.cls = paste(output.dir, 
                "temp1.cls", sep = ""), output.ds = paste(output.dir, 
                identifier, ".all.H.gct", sep = ""), output.cls = paste(output.dir, 
                identifier, ".all.H.cls", sep = ""))
            O <- MSIG.Match.and.Merge(input1.ds = paste(output.dir, 
                identifier, ".all.gct", sep = ""), input1.cls = paste(output.dir, 
                identifier, ".all.cls", sep = ""), input2.ds = paste(output.dir, 
                "temp4.gct", sep = ""), input2.cls = paste(output.dir, 
                "temp1.cls", sep = ""), output.ds = paste(output.dir, 
                identifier, ".all.gct", sep = ""), output.cls = paste(output.dir, 
                identifier, ".all.cls", sep = ""))
            print(c("Done processing test file: ", test.datasets.table[[ds]]$gct.file))
        }
    }
    CLS <- MSIG.ReadClsFile(file = paste(output.dir, identifier, 
        ".model_set.2.cls", sep = ""))
    model.size <- length(CLS$class.v)
    O <- MSIG.Projection.Plots.3(input.ds = paste(output.dir, 
        identifier, ".all.H.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".all.H.cls", sep = ""), model.set = seq(1, 
        model.size), output.2D.proj.file = paste(output.dir, 
        identifier, ".2D.proj.gct", sep = ""), output.2D.proj.plot = paste(output.dir, 
        identifier, ".2D.proj", sep = ""), output.3D.proj.file = paste(output.dir, 
        identifier, ".3D.proj.gct", sep = ""), output.3D.1.proj.plot = paste(output.dir, 
        identifier, ".3D.1.proj", sep = ""), output.3D.2.proj.plot = paste(output.dir, 
        identifier, ".3D.2.proj", sep = ""), output.3D.3.proj.plot = paste(output.dir, 
        identifier, ".3D.3.proj", sep = ""), output.heatmap.plot = paste(output.dir, 
        identifier, ".heatmap", sep = ""), output.heatmap.sorted.plot = paste(output.dir, 
        identifier, ".heatmap.sorted", sep = ""), output.heatmap.sorted.2.plot = paste(output.dir, 
        identifier, ".heatmap.sorted.2", sep = ""), output.hclust.plot = paste(output.dir, 
        identifier, ".hclust", sep = ""), use.biplot = use.biplot, 
        title = identifier, seed = seed, non.interactive.run = non.interactive.run, 
        heatmap.row.norm = heatmap.row.norm, heatmap.cmap.type = heatmap.cmap.type, 
        symbol.scaling = symbol.scaling, col = col, symbs = symbs)
    O <- MSIG.Evaluate.Projection.2(input.ds = paste(output.dir, 
        identifier, ".all.H.gct", sep = ""), input.cls = paste(output.dir, 
        identifier, ".all.H.cls", sep = ""), model.set = seq(1, 
        model.size), prediction.results.file = paste(output.dir, 
        identifier, ".pred.txt", sep = ""), prediction.matrix.file = paste(output.dir, 
        identifier, ".pred.gct", sep = ""), train.pred.plot = paste(output.dir, 
        identifier, ".train.pred", sep = ""), test.pred.plot = paste(output.dir, 
        identifier, ".test.pred", sep = ""), pred.2D.plot = paste(output.dir, 
        identifier, ".2D.pred", sep = ""), col = col, use.feature.names = use.feature.names, 
        nchar.phen = nchar.phen, high.conf.thres = high.conf.thres, 
        symbol.scaling = symbol.scaling, symbs = symbs, levels = levels, 
        nlevels = nlevels, kernel = kernel, cost = cost, gamma = gamma, 
        produce.contours = produce.contours)
    input.ds <- paste(output.dir, identifier, ".all.gct", sep = "")
    input.cls <- paste(output.dir, identifier, ".all.cls", sep = "")
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
        }
        mycols <- col[class.labels[HC$order]]
        i <- 0
    })
    dL <- dendrapply(dhc, colLab)
    plot.filename <- paste(output.dir, identifier, ".htree", 
        sep = "")
    input.ds <- paste(output.dir, identifier, ".all.H.gct", sep = "")
    input.cls <- paste(output.dir, identifier, ".all.H.cls", 
        sep = "")
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
        }
        mycols <- col[class.labels[HC$order]]
        i <- 0
    })
    dL <- dendrapply(dhc, colLab)
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
