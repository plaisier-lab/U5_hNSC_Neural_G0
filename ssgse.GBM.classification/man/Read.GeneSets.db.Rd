\name{Read.GeneSets.db}
\alias{Read.GeneSets.db}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read.GeneSets.db
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Read.GeneSets.db(gs.db, thres.min = 2, thres.max = 2000, gene.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gs.db}{
%%     ~~Describe \code{gs.db} here~~
}
  \item{thres.min}{
%%     ~~Describe \code{thres.min} here~~
}
  \item{thres.max}{
%%     ~~Describe \code{thres.max} here~~
}
  \item{gene.names}{
%%     ~~Describe \code{gene.names} here~~
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
function (gs.db, thres.min = 2, thres.max = 2000, gene.names = NULL) 
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
        if (is.null(gene.names)) {
            existing.set <- rep(TRUE, length(gene.set.tags))
        }
        else {
            existing.set <- is.element(gene.set.tags, gene.names)
        }
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < thres.min) || (set.size > thres.max)) 
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
    return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, 
        size.G = size.G, max.N.gs = max.Ng))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
