\name{GSEA.VarFilter}
\alias{GSEA.VarFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.VarFilter
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.VarFilter(V, fold, delta, gene.names = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{fold}{
%%     ~~Describe \code{fold} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
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
function (V, fold, delta, gene.names = "") 
{
    cols <- length(V[1, ])
    rows <- length(V[, 1])
    row.max <- apply(V, MARGIN = 1, FUN = max)
    row.min <- apply(V, MARGIN = 1, FUN = min)
    flag <- array(dim = rows)
    flag <- (row.max/row.min >= fold) & (row.max - row.min >= 
        delta)
    size <- sum(flag)
    B <- matrix(0, nrow = size, ncol = cols)
    j <- 1
    if (length(gene.names) == 1) {
        for (i in 1:rows) {
            if (flag[i]) {
                B[j, ] <- V[i, ]
                j <- j + 1
            }
        }
        return(B)
    }
    else {
        new.list <- vector(mode = "character", length = size)
        for (i in 1:rows) {
            if (flag[i]) {
                B[j, ] <- V[i, ]
                new.list[j] <- gene.names[i]
                j <- j + 1
            }
        }
        return(list(V = B, new.list = new.list))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
