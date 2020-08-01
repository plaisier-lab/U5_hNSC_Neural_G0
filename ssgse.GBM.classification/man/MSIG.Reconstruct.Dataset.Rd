\name{MSIG.Reconstruct.Dataset}
\alias{MSIG.Reconstruct.Dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Reconstruct.Dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Reconstruct.Dataset(input.H.ds, input.W.ds, output.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.H.ds}{
%%     ~~Describe \code{input.H.ds} here~~
}
  \item{input.W.ds}{
%%     ~~Describe \code{input.W.ds} here~~
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
function (input.H.ds, input.W.ds, output.file) 
{
    library(MASS)
    print(c("Running MSIG.Reconstruct.Dataset... on: ", input.H.ds, 
        input.W.ds))
    dataset <- MSIG.Gct2Frame(filename = input.W.ds)
    W <- data.matrix(dataset$ds)
    W.row.names <- dataset$row.names
    W.row.descs <- dataset$descs
    W.names <- dataset$names
    dataset <- MSIG.Gct2Frame(filename = input.H.ds)
    H <- data.matrix(dataset$ds)
    H.row.names <- dataset$row.names
    H.row.descs <- dataset$descs
    H.names <- dataset$names
    A <- W \%*\% H
    V <- data.frame(A)
    names(V) <- H.names
    row.names(V) <- W.row.names
    write.gct(gct.data.frame = V, filename = output.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
