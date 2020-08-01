\name{OPAM.project.dataset}
\alias{OPAM.project.dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.project.dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.project.dataset(input.ds, output.ds, gene.set.databases, gene.set.selection = "ALL", sample.norm.type = "rank", output.score.type = "ES")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{gene.set.databases}{
%%     ~~Describe \code{gene.set.databases} here~~
}
  \item{gene.set.selection}{
%%     ~~Describe \code{gene.set.selection} here~~
}
  \item{sample.norm.type}{
%%     ~~Describe \code{sample.norm.type} here~~
}
  \item{output.score.type}{
%%     ~~Describe \code{output.score.type} here~~
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
function (input.ds, output.ds, gene.set.databases, gene.set.selection = "ALL", 
    sample.norm.type = "rank", output.score.type = "ES") 
{
    library(gtools)
    library(verification)
    library(RColorBrewer)
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gene.names <- dataset$row.names
    gene.descs <- dataset$descs
    sample.names <- dataset$names
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    temp <- strsplit(input.ds, split = "/")
    s <- length(temp[[1]])
    input.file.name <- temp[[1]][s]
    temp <- strsplit(input.file.name, split = ".gct")
    input.file.prefix <- temp[[1]][1]
    if (sample.norm.type == "rank") {
        for (j in 1:Ns) {
            m[, j] <- rank(m[, j], ties.method = "average")
        }
        m <- 10000 * m/Ng
    }
    else if (sample.norm.type == "log.rank") {
        for (j in 1:Ns) {
            m[, j] <- rank(m[, j], ties.method = "average")
        }
        m <- log(10000 * m/Ng + exp(1))
    }
    else if (sample.norm.type == "log") {
        m[m < 1] <- 1
        m <- log(m + exp(1))
    }
    max.G <- 0
    max.N <- 0
    for (gsdb in gene.set.databases) {
        GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, 
            gene.names = NULL)
        max.G <- max(max.G, max(GSDB$size.G))
        max.N <- max.N + GSDB$N.gs
    }
    N.gs <- 0
    gs <- matrix("null", nrow = max.N, ncol = max.G)
    gs.names <- vector(length = max.N, mode = "character")
    gs.desc <- vector(length = max.N, mode = "character")
    size.G <- vector(length = max.N, mode = "numeric")
    start <- 1
    for (gsdb in gene.set.databases) {
        GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, 
            gene.names = NULL)
        N.gs <- GSDB$N.gs
        gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
        gs.desc[start:(start + N.gs - 1)] <- GSDB$gs.desc
        size.G[start:(start + N.gs - 1)] <- GSDB$size.G
        gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 
            1:max(GSDB$size.G)]
        start <- start + N.gs
    }
    N.gs <- max.N
    if (gene.set.selection[1] != "ALL") {
        locs <- match(gene.set.selection, gs.names)
        N.gs <- sum(!is.na(locs))
        gs <- gs[locs, ]
        gs.names <- gs.names[locs]
        gs.desc <- gs.desc[locs]
        size.G <- size.G[locs]
    }
    score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
    for (gs.i in 1:N.gs) {
        gene.set <- gs[gs.i, 1:size.G[gs.i]]
        gene.overlap <- intersect(gene.set, gene.names)
        print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", 
            length(gene.overlap)))
        if (length(gene.overlap) == 0) {
            score.matrix[gs.i, ] <- runif(Ns, min = 1e-06, max = 1.1e-06)
            next
        }
        else {
            gene.set.locs <- match(gene.overlap, gene.set)
            gene.names.locs <- match(gene.overlap, gene.names)
            msig <- m[gene.names.locs, ]
            msig.names <- gene.names[gene.names.locs]
            if (output.score.type == "ES") {
                OPAM <- OPAM.Projection(data.array = m, gene.names = gene.names, 
                  n.cols = Ns, n.rows = Ng, weight = 0.25, statistic = "area.under.RES", 
                  gene.set = gene.overlap, nperm = 1)
                score.matrix[gs.i, ] <- OPAM$ES.vector
            }
            else if (output.score.type == "NES") {
                OPAM <- OPAM.Projection(data.array = m, gene.names = gene.names, 
                  n.cols = Ns, n.rows = Ng, weight = 0.25, statistic = "area.under.RES", 
                  gene.set = gene.overlap, nperm = 100)
                score.matrix[gs.i, ] <- OPAM$NES.vector
            }
        }
    }
    V.GCT <- data.frame(score.matrix)
    names(V.GCT) <- sample.names
    row.names(V.GCT) <- gs.names
    write.gct(gct.data.frame = V.GCT, descs = gs.desc, filename = output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
