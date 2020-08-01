\name{MSIG.StJude2gct}
\alias{MSIG.StJude2gct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.StJude2gct
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.StJude2gct(input.ds, output.ds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
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
function (input.ds, output.ds) 
{
    header.cont <- readLines(input.ds, n = 1)
    temp <- unlist(strsplit(header.cont, "\t"))
    colst <- length(temp)
    header.labels <- temp[seq(3, colst, 3)]
    ds <- read.delim(input.ds, header = F, row.names = 1, sep = "\t", 
        skip = 2, blank.lines.skip = T, comment.char = "", as.is = T)
    colst <- length(ds[1, ])
    cols <- (colst - 1)/3
    rows <- length(ds[, 1])
    A <- matrix(nrow = rows - 1, ncol = cols)
    A <- ds[1:rows, seq(2, colst, 3)]
    descs <- ds[, 1]
    table1 <- data.frame(A)
    names(table1) <- header.labels
    write.gct(gct.data.frame = table1, descs = descs, output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
