\name{MSIG.ReadClsFile}
\alias{MSIG.ReadClsFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.ReadClsFile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.ReadClsFile(file = "NULL")
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
    class.list <- unlist(strsplit(cls.cont[[3]], " "))
    s <- length(class.list)
    t <- table(class.list)
    l <- length(t)
    phen <- vector(length = l, mode = "character")
    class.v <- vector(length = s, mode = "numeric")
    current.label <- class.list[1]
    current.number <- 1
    class.v[1] <- current.number
    phen[1] <- current.label
    phen.count <- 1
    if (length(class.list) > 1) {
        for (i in 2:s) {
            if (class.list[i] == current.label) {
                class.v[i] <- current.number
            }
            else {
                phen.count <- phen.count + 1
                current.number <- current.number + 1
                current.label <- class.list[i]
                phen[phen.count] <- current.label
                class.v[i] <- current.number
            }
        }
    }
    return(list(phen = phen, class.v = class.v, class.list = class.list))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
