\name{MSIG.Extract.Features}
\alias{MSIG.Extract.Features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Extract.Features
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Extract.Features(input.ds, features, output.ds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{features}{
%%     ~~Describe \code{features} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
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
function (input.ds, features, output.ds) 
{
    print(c("Running MSIG.Extract.Features... on dataset:", input.ds))
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    A <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dim(A)
    cols <- length(A[1, ])
    rows <- length(A[, 1])
    temp <- readLines(features)
    gene.set.size <- length(unlist(strsplit(temp, "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp, "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[1]
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
        gene.set.tags[j] <- gs.line[j + 2]
    }
    locations <- match(gene.set.tags, gs.names, nomatch = 0)
    A2 <- A[locations, ]
    A2 <- data.frame(A2)
    names(A2) <- sample.names
    row.names(A2) <- gene.set.tags
    write.gct(A2, filename = output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
