\name{MSIG.Resort.Rows}
\alias{MSIG.Resort.Rows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Resort.Rows
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Resort.Rows(input.ds, new.order, new.row.labs = NULL, new.row.descs = NULL, output.ds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{new.order}{
%%     ~~Describe \code{new.order} here~~
}
  \item{new.row.labs}{
%%     ~~Describe \code{new.row.labs} here~~
}
  \item{new.row.descs}{
%%     ~~Describe \code{new.row.descs} here~~
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
function (input.ds, new.order, new.row.labs = NULL, new.row.descs = NULL, 
    output.ds) 
{
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    m <- m[new.order, ]
    if (length(new.row.labs) == 1) {
        gs.names <- gs.names[new.order]
        gs.descs <- gs.descs[new.order]
    }
    else {
        gs.names <- new.row.labs
        gs.names <- new.row.descs
    }
    V <- data.frame(m)
    names(V) <- sample.names
    row.names(V) <- gs.names
    write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
