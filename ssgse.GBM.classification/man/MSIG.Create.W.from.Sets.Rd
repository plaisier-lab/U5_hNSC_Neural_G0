\name{MSIG.Create.W.from.Sets}
\alias{MSIG.Create.W.from.Sets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Create.W.from.Sets
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Create.W.from.Sets(gs.db, output.file, non.interactive.run = F, gs.size.threshold.min = 5, gs.size.threshold.max = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gs.db}{
%%     ~~Describe \code{gs.db} here~~
}
  \item{output.file}{
%%     ~~Describe \code{output.file} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{gs.size.threshold.min}{
%%     ~~Describe \code{gs.size.threshold.min} here~~
}
  \item{gs.size.threshold.max}{
%%     ~~Describe \code{gs.size.threshold.max} here~~
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
function (gs.db, output.file, non.interactive.run = F, gs.size.threshold.min = 5, 
    gs.size.threshold.max = 1e+05) 
{
    print(" *** Running MSIG.Create.W.from.Sets...")
    if (.Platform$OS.type == "windows") {
        memory.limit(6e+09)
        memory.limit()
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
        set.size <- length(gene.set.tags)
        if ((set.size < gs.size.threshold.min) || (set.size > 
            gs.size.threshold.max)) 
            next
        temp.size.G[gs.count] <- set.size
        gs[gs.count, ] <- c(gene.set.tags, rep("null", max.size.G - 
            temp.size.G[gs.count]))
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
    print(c("Number of Gene Sets:", Ng))
    print(c("Original number of Gene Sets:", max.Ng))
    print(c("Maximum gene set size:", max.size.G))
    all.genes <- NULL
    for (i in 1:Ng) {
        gene.set <- gs[i, gs[i, ] != "null"]
        all.genes <- union(all.genes, gene.set)
    }
    W <- matrix(0, nrow = length(all.genes), ncol = Ng)
    for (i in 1:Ng) {
        print(paste("Computing W column for gene set:", i, sep = " "))
        gene.set <- gs[i, gs[i, ] != "null"]
        W[, i] <- as.real(sign(match(all.genes, gene.set, nomatch = 0)))
    }
    V <- data.frame(W)
    names(V) <- gs.names
    row.names(V) <- all.genes
    write.gct(gct.data.frame = V, filename = output.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
