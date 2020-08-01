\name{MSIG.GSEA.Project}
\alias{MSIG.GSEA.Project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.GSEA.Project
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.GSEA.Project(input.ds, input.cls = "", gs.db, output.gct.file, output.plot.file, non.interactive.run = F, nperm = 25, projection.type = 0, weighted.score.type = 1, topgs = 10, preproc.type = 0, gs.size.threshold.min = 5, gs.size.threshold.max = 1e+05, reverse.sign = F, seed = 3338)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{gs.db}{
%%     ~~Describe \code{gs.db} here~~
}
  \item{output.gct.file}{
%%     ~~Describe \code{output.gct.file} here~~
}
  \item{output.plot.file}{
%%     ~~Describe \code{output.plot.file} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{projection.type}{
%%     ~~Describe \code{projection.type} here~~
}
  \item{weighted.score.type}{
%%     ~~Describe \code{weighted.score.type} here~~
}
  \item{topgs}{
%%     ~~Describe \code{topgs} here~~
}
  \item{preproc.type}{
%%     ~~Describe \code{preproc.type} here~~
}
  \item{gs.size.threshold.min}{
%%     ~~Describe \code{gs.size.threshold.min} here~~
}
  \item{gs.size.threshold.max}{
%%     ~~Describe \code{gs.size.threshold.max} here~~
}
  \item{reverse.sign}{
%%     ~~Describe \code{reverse.sign} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
function (input.ds, input.cls = "", gs.db, output.gct.file, output.plot.file, 
    non.interactive.run = F, nperm = 25, projection.type = 0, 
    weighted.score.type = 1, topgs = 10, preproc.type = 0, gs.size.threshold.min = 5, 
    gs.size.threshold.max = 1e+05, reverse.sign = F, seed = 3338) 
{
    print(" *** Running GSEA projection...")
    if (.Platform$OS.type == "windows") {
    }
    set.seed(seed = seed, kind = NULL)
    if (is.data.frame(input.ds)) {
        dataset <- input.ds
    }
    else {
        if (regexpr(pattern = ".gct", input.ds) == -1) {
            dataset <- GSEA.Res2Frame(filename = input.ds)
        }
        else {
            dataset <- GSEA.Gct2Frame(filename = input.ds)
        }
    }
    gene.labels <- row.names(dataset)
    sample.names <- names(dataset)
    A <- data.matrix(dataset)
    dim(A)
    cols <- length(A[1, ])
    rows <- length(A[, 1])
    if (preproc.type == 1) {
        A <- GSEA.Threshold(A, 0, 1e+05)
        A <- GSEA.NormalizeCols(A)
    }
    else if (preproc.type == 2) {
        A <- GSEA.Threshold(A, 0, 1e+05)
        for (j in 1:cols) {
            A[, j] <- rank(A[, j])
        }
    }
    else if (preproc.type == 3) {
    }
    if (input.cls == "") {
        class.labels <- rep(0, cols)
        class.phen <- "   "
    }
    else {
        CLS <- ReadClsFile(file = input.cls)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
    }
    if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
        temp <- gs.db
    }
    else {
        temp <- readLines(gs.db)
    }
    max.Ng <- length(temp)
    temp.size.G <- vector(length = max.Ng, mode = "numeric")
    for (i in 1:max.Ng) {
        temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 
            2
    }
    max.size.G <- max(temp.size.G)
    gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, 
        ncol = max.size.G)
    temp.names <- vector(length = max.Ng, mode = "character")
    temp.desc <- vector(length = max.Ng, mode = "character")
    gs.count <- 1
    for (i in 1:max.Ng) {
        gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 
            2
        gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
        gene.set.name <- gs.line[1]
        gene.set.desc <- gs.line[1]
        gene.set.tags <- vector(length = gene.set.size, mode = "character")
        for (j in 1:gene.set.size) {
            gene.set.tags[j] <- gs.line[j + 2]
        }
        existing.set <- is.element(gene.set.tags, gene.labels)
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < gs.size.threshold.min) || (set.size > 
            gs.size.threshold.max)) 
            next
        temp.size.G[gs.count] <- set.size
        gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", 
            max.size.G - temp.size.G[gs.count]))
        temp.names[gs.count] <- gene.set.name
        temp.desc[gs.count] <- gene.set.desc
        gs.count <- gs.count + 1
    }
    Ng <- gs.count - 1
    gs.names <- vector(length = Ng, mode = "character")
    gs.desc <- vector(length = Ng, mode = "character")
    size.G <- vector(length = Ng, mode = "numeric")
    gs.names <- temp.names[1:Ng]
    gs.desc <- temp.desc[1:Ng]
    size.G <- temp.size.G[1:Ng]
    N <- length(A[, 1])
    Ns <- length(A[1, ])
    print(c("Number of genes:", N))
    print(c("Number of Gene Sets:", Ng))
    print(c("Number of samples:", Ns))
    print(c("Original number of Gene Sets:", max.Ng))
    print(c("Maximum gene set size:", max.size.G))
    sample.ES <- matrix(0, nrow = Ng, ncol = Ns)
    sample.NES <- matrix(0, nrow = Ng, ncol = Ns)
    sample.p.val <- matrix(0, nrow = Ng, ncol = Ns)
    correl.vector <- vector(length = N, mode = "numeric")
    order.matrix <- matrix(nrow = N, ncol = Ns)
    all.gs.descs <- vector(length = Ng, mode = "character")
    phi <- array(0, c(Ns, Ng, nperm))
    for (sample.index in 1:Ns) {
        order.matrix[, sample.index] <- order(A[, sample.index], 
            decreasing = T)
        gene.list2 <- order.matrix[, sample.index]
        for (i in 1:Ng) {
            print(paste("Computing observed enrichment for gene set:", 
                i, "in sample:", sample.index, sep = " "))
            gene.set <- gs[i, gs[i, ] != "null"]
            gene.set2 <- vector(length = length(gene.set), mode = "numeric")
            gene.set2 <- match(gene.set, gene.labels)
            if (weighted.score.type == 0) {
                correl.vector <- rep(1, N)
            }
            else if (weighted.score.type > 0) {
                correl.vector <- A[gene.list2, sample.index]
            }
            if (projection.type == 0) {
                GSEA.results <- GSEA.EnrichmentScore2(gene.list = gene.list2, 
                  gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                  correl.vector = correl.vector)
            }
            else if (projection.type == 1) {
                GSEA.results <- GSEA.EnrichmentScore3(gene.list = gene.list2, 
                  gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                  correl.vector = correl.vector)
            }
            else if (projection.type == 2) {
                GSEA.results <- Wilcox.Score(gene.list = gene.list2, 
                  gene.set = gene.set2)
            }
            sample.ES[i, sample.index] <- GSEA.results$ES
            for (r in 1:nperm) {
                reshuffled.gene.labels <- sample(1:rows)
                if (weighted.score.type == 0) {
                  correl.vector <- rep(1, N)
                }
                else if (weighted.score.type == 1) {
                  correl.vector <- A[reshuffled.gene.labels, 
                    sample.index]
                }
                if (projection.type == 0) {
                  GSEA.results <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = correl.vector)
                }
                else if (projection.type == 1) {
                  GSEA.results <- GSEA.EnrichmentScore3(gene.list = reshuffled.gene.labels, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = correl.vector)
                }
                else if (projection.type == 2) {
                  GSEA.results <- Wilcox.Score(gene.list = reshuffled.gene.labels, 
                    gene.set = gene.set2)
                }
                phi[sample.index, i, r] <- GSEA.results$ES
            }
            if (sample.ES[i, sample.index] >= 0) {
                pos.phi <- phi[sample.index, i, phi[sample.index, 
                  i, ] >= 0]
                if (length(pos.phi) == 0) {
                  pos.phi <- 0.5
                }
                pos.m <- mean(pos.phi)
                sample.NES[i, sample.index] <- sample.ES[i, sample.index]/pos.m
                s <- sum(pos.phi >= sample.ES[i, sample.index])/length(pos.phi)
                s <- ifelse(s == 0, 1/nperm, s)
                sample.p.val[i, sample.index] <- 1/s
            }
            else {
                neg.phi <- phi[sample.index, i, phi[sample.index, 
                  i, ] < 0]
                if (length(neg.phi) == 0) {
                  neg.phi <- 0.5
                }
                neg.m <- mean(neg.phi)
                sample.NES[i, sample.index] <- sample.ES[i, sample.index]/abs(neg.m)
                s <- sum(neg.phi <= sample.ES[i, sample.index])/length(neg.phi)
                s <- ifelse(s == 0, 1/nperm, s)
                sample.p.val[i, sample.index] <- -1/s
            }
        }
    }
    phen.cmap <- c("grey3", "cadetblue", "darkgreen", "chartreuse2", 
        "red1", "darkred", "orange", "blue2", "lightblue", "pink", 
        "coral")
    MSIG.HeatMapPlot.4(V = sample.ES, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, phen.cmap = phen.cmap, col.names = sample.names, 
        main = "Single Sample GSEA Projection", xlab = " ", ylab = " ", 
        sub = gct.file, row.norm = T, cmap.type = 4, char.rescale = 0.75)
    savePlot(filename = output.plot.file, type = "jpeg", device = dev.cur())
    print("Saving projected dataset...")
    V <- data.frame(sample.NES)
    names(V) <- sample.names
    row.names(V) <- gs.names
    write.gct(gct.data.frame = V, descs = all.gs.descs, filename = output.gct.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
