\name{MSIG.Replace.Missing.Values}
\alias{MSIG.Replace.Missing.Values}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Replace.Missing.Values
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Replace.Missing.Values(input.ds, output.ds, miss.val.replacement)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{miss.val.replacement}{
%%     ~~Describe \code{miss.val.replacement} here~~
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
function (input.ds, output.ds, miss.val.replacement) 
{
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    data <- dataset$ds
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    N <- length(data[, 1])
    M <- length(data[1, ])
    for (i in 1:N) {
        for (j in 1:M) {
            if (is.na(data[i, j])) {
                data[i, j] <- miss.val.replacement
            }
        }
    }
    m <- data.matrix(data)
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
