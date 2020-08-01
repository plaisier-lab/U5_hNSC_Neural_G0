\name{MSIG.Reorder.Dataset}
\alias{MSIG.Reorder.Dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Reorder.Dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Reorder.Dataset(input.ds, input.cls, phen.reorder, output.ds, output.cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{phen.reorder}{
%%     ~~Describe \code{phen.reorder} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{output.cls}{
%%     ~~Describe \code{output.cls} here~~
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
function (input.ds, input.cls, phen.reorder, output.ds, output.cls) 
{
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    CLS <- ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    new.phen <- class.phen[phen.reorder]
    new.labels <- vector(length(class.labels), mode = "numeric")
    for (i in 1:length(class.labels)) {
        new.labels[i] <- phen.reorder[class.labels[i]]
    }
    col.index <- order(new.labels, decreasing = F)
    new.labels <- new.labels[col.index]
    sample.names <- sample.names[col.index]
    rows <- length(m[, 1])
    for (j in 1:rows) {
        m[j, ] <- m[j, col.index]
    }
    names(m) <- sample.names
    V <- data.frame(m)
    names(V) <- sample.names
    row.names(V) <- gs.names
    write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
    write.cls(class.v = new.labels, phen = new.phen, filename = output.cls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
