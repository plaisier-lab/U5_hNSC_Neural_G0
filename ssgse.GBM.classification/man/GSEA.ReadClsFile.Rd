\name{GSEA.ReadClsFile}
\alias{GSEA.ReadClsFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.ReadClsFile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.ReadClsFile(file = "NULL")
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
    phen.label <- vector(length = l, mode = "numeric")
    class.v <- vector(length = s, mode = "numeric")
    for (i in 1:l) {
        phen[i] <- noquote(names(t)[i])
        phen.label[i] <- i - 1
    }
    for (i in 1:s) {
        for (j in 1:l) {
            if (class.list[i] == phen[j]) {
                class.v[i] <- phen.label[j]
            }
        }
    }
    return(list(phen = phen, class.v = class.v))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
