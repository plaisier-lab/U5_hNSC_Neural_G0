\name{MSIG.Factors.Project.old}
\alias{MSIG.Factors.Project.old}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Factors.Project.old
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Factors.Project.old(input.ds, factors.ds, postprojnorm = TRUE, output.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{factors.ds}{
%%     ~~Describe \code{factors.ds} here~~
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
function (input.ds, factors.ds, postprojnorm = TRUE, output.file) 
{
    library(MASS)
    print(c("Running MSIG.Factors.Project... on: ", input.ds))
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
    print(c("Size of Input dataset=", length(gs.names), " genes"))
    print(c("Size of W matrix (rows)=", length(W.row.names), 
        " genes"))
    print(c("Size of overlap=", length(overlap), " genes"))
    locations.m <- match(overlap, gs.names, nomatch = 0)
    m2 <- m[locations.m, ]
    locations.W <- match(overlap, W.row.names, nomatch = 0)
    W2 <- W[locations.W, ]
    H <- ginv(W2) \%*\% m2
    max.H <- max(H)
    min.H <- min(H)
    H <- (H - min.H)/(max.H - min.H)
    if (postprojnorm == TRUE) {
        n.col <- length(H[1, ])
        for (i in 1:n.col) {
            S.2 <- sqrt(sum(H[, i] * H[, i]))
            H[, i] <- H[, i]/S.2
        }
    }
    V <- data.frame(H)
    names(V) <- sample.names
    row.names(V) <- W.names
    write.gct(gct.data.frame = V, filename = output.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
