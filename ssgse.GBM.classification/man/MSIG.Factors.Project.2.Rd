\name{MSIG.Factors.Project.2}
\alias{MSIG.Factors.Project.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Factors.Project.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Factors.Project.2(input.ds, factors.ds, inverse.type = "transpose", weighting.type = -1, postprojnorm = FALSE, output.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{factors.ds}{
%%     ~~Describe \code{factors.ds} here~~
}
  \item{inverse.type}{
%%     ~~Describe \code{inverse.type} here~~
}
  \item{weighting.type}{
%%     ~~Describe \code{weighting.type} here~~
}
  \item{postprojnorm}{
%%     ~~Describe \code{postprojnorm} here~~
}
  \item{output.file}{
%%     ~~Describe \code{output.file} here~~
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
function (input.ds, factors.ds, inverse.type = "transpose", weighting.type = -1, 
    postprojnorm = FALSE, output.file) 
{
    library(MASS)
    print("Running MSIG.Factors.Project..")
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dataset <- MSIG.Gct2Frame(filename = factors.ds)
    W <- data.matrix(dataset$ds)
    W.row.names <- dataset$row.names
    W.row.descs <- dataset$descs
    W.names <- dataset$names
    overlap <- intersect(gs.names, W.row.names)
    locations <- match(overlap, W.row.names, nomatch = 0)
    W2 <- W[locations, ]
    if (inverse.type == "transpose") {
        if (weighting.type == -1) {
            m2 <- m[locations, ]
        }
        else if (weighting.type == 0) {
            m.rank <- MSIG.NormalizeCols.Rank(m)
            m2 <- m.rank[locations, ]
        }
        else if (weighting.type == 0) {
            m.rank <- MSIG.NormalizeCols.Rank(m)
            maxm <- max(m)
            minm <- min(m)
            m.weight <- (m - min(m))/(maxm - minm)
            m2 <- m.rank[locations, ] * m.weight[locations, ]
        }
    }
    if (inverse.type == "Penrose-Moore") {
        H <- ginv(W2) \%*\% m2
    }
    else if (inverse.type == "transpose") {
        H <- t(W2) \%*\% m2
    }
    if (postprojnorm == TRUE) {
        n.col <- length(H[1, ])
        for (i in 1:n.col) {
            S.2 <- sqrt(sum(H[, i] * H[, i]))
            norm <- 1/S.2
            H[, i] <- H[, i]/norm
        }
    }
    V <- data.frame(H)
    names(V) <- sample.names
    row.names(V) <- W.names
    write.gct(gct.data.frame = V, descs = gs.descs, filename = output.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
