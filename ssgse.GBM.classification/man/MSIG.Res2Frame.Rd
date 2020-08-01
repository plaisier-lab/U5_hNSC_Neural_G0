\name{MSIG.Res2Frame}
\alias{MSIG.Res2Frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Res2Frame
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Res2Frame(filename = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
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
function (filename = "NULL") 
{
    header.cont <- readLines(filename, n = 1)
    temp <- unlist(strsplit(header.cont, "\t"))
    colst <- length(temp)
    header.labels <- temp[seq(3, colst, 2)]
    ds <- read.delim(filename, header = F, row.names = 2, sep = "\t", 
        skip = 3, blank.lines.skip = T, comment.char = "", as.is = T)
    colst <- length(ds[1, ])
    cols <- (colst - 1)/2
    rows <- length(ds[, 1])
    A <- matrix(nrow = rows - 1, ncol = cols)
    A <- ds[1:rows, seq(2, colst, 2)]
    table1 <- data.frame(A)
    names(table1) <- header.labels
    names <- header.labels
    row.names <- row.names(A)
    descs <- ds[, 1]
    return(list(ds = table1, row.names = row.names, descs = descs, 
        names = names))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
