\name{MSIG.ReadPhenFile}
\alias{MSIG.ReadPhenFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.ReadPhenFile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.ReadPhenFile(file = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
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
function (file = "NULL") 
{
    cls.cont <- readLines(file)
    num.lines <- length(cls.cont)
    temp <- unlist(strsplit(cls.cont[[1]], " "))
    if (length(temp) == 3) {
        phen.names <- NULL
        col.phen <- NULL
    }
    else {
        l.phen.names <- match("phen.names:", temp)
        l.col.phen <- match("col.phen:", temp)
        phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
        col.phen <- temp[(l.col.phen + 1):length(temp)]
    }
    temp <- unlist(strsplit(cls.cont[[2]], " "))
    phen.list <- temp[2:length(temp)]
    for (k in 1:(num.lines - 2)) {
        temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
        if (k == 1) {
            len <- length(temp)
            class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
            class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
            phen <- NULL
        }
        class.list[k, ] <- temp
        classes <- unique(temp)
        class.v[k, ] <- match(temp, classes)
        phen[[k]] <- classes
    }
    if (num.lines == 3) {
        class.list <- as.vector(class.list)
        class.v <- as.vector(class.v)
        phen <- unlist(phen)
    }
    return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, 
        col.phen = col.phen, class.v = class.v, class.list = class.list))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
