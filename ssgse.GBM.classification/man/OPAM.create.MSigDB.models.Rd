\name{OPAM.create.MSigDB.models}
\alias{OPAM.create.MSigDB.models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.create.MSigDB.models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.create.MSigDB.models(gs.db, models.dir, gs.size.threshold.min = 10, gs.size.threshold.max = 500, source = "MSigDB")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gs.db}{
%%     ~~Describe \code{gs.db} here~~
}
  \item{models.dir}{
%%     ~~Describe \code{models.dir} here~~
}
  \item{gs.size.threshold.min}{
%%     ~~Describe \code{gs.size.threshold.min} here~~
}
  \item{gs.size.threshold.max}{
%%     ~~Describe \code{gs.size.threshold.max} here~~
}
  \item{source}{
%%     ~~Describe \code{source} here~~
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
function (gs.db, models.dir, gs.size.threshold.min = 10, gs.size.threshold.max = 500, 
    source = "MSigDB") 
{
    temp <- readLines(gs.db)
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
        gene.set.desc <- gs.line[2]
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
    for (i in 1:Ng) {
        print(paste("Creating model for gene set:", i, gs.names[i], 
            sep = " "))
        gene.set <- gs[i, gs[i, ] != "null"]
        gene.set.string <- paste("c('", paste(gene.set, collapse = "','"), 
            "')\n", sep = "")
        gene.set.name <- gs.names[i]
        gene.set.desc <- gs.desc[i]
        m.file <- paste(models.dir, gene.set.name, ".mod", sep = "")
        cat("model.creation.date", paste("'", date(), "'\n", 
            sep = ""), file = m.file, append = FALSE, sep = "\t")
        cat("model.name", paste("'", gene.set.name, "'\n", sep = ""), 
            file = m.file, append = TRUE, sep = "\t")
        cat("model.description", paste("'", source, "'\n", sep = ""), 
            file = m.file, append = TRUE, sep = "\t")
        cat("sample.norm.type", "'rank'\n", file = m.file, append = TRUE, 
            sep = "\t")
        cat("statistic", "'area.under.RES'\n", file = m.file, 
            append = TRUE, sep = "\t")
        cat("weight", paste(0.25, "\n", sep = ""), file = m.file, 
            append = TRUE, sep = "\t")
        cat("random.seed", paste(12345, "\n", sep = ""), file = m.file, 
            append = TRUE, sep = "\t")
        cat("nperm", paste(0, "\n", sep = ""), file = m.file, 
            append = TRUE, sep = "\t")
        cat("link.function", "'logit'\n", file = m.file, append = TRUE, 
            sep = "\t")
        cat("c1", "c('black','lightgrey')\n", file = m.file, 
            append = TRUE, sep = "\t")
        cat("msig.up.genes", gene.set.string, file = m.file, 
            append = TRUE, sep = "\t")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
