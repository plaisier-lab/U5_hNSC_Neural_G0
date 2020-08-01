\name{MSIG.replace.row.col}
\alias{MSIG.replace.row.col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.replace.row.col
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.replace.row.col(input.ds, output.ds, mode = "row", number, values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
}
  \item{number}{
%%     ~~Describe \code{number} here~~
}
  \item{values}{
%%     ~~Describe \code{values} here~~
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
function (input.ds, output.ds, mode = "row", number, values) 
{
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    if (mode == "row") {
        m[number, ] <- ifelse(length(values) == 1, rep(values, 
            length(m[number, ])), values)
    }
    else if (mode == "col") {
        m[, number] <- ifelse(length(values) == 1, rep(values, 
            length(m[, number])), values)
    }
    else {
        stop(c("unknown mode:", mode))
    }
    V <- data.frame(m)
    names(V) <- dataset$names
    row.names(V) <- dataset$row.names
    write.gct(gct.data.frame = V, descs = dataset$descs, filename = output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
