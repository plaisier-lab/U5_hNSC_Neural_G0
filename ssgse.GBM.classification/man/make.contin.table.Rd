\name{make.contin.table}
\alias{make.contin.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make.contin.table
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.contin.table(pred, actual)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{
%%     ~~Describe \code{pred} here~~
}
  \item{actual}{
%%     ~~Describe \code{actual} here~~
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
function (pred, actual) 
{
    tab <- table(pred, actual)
    pred.classes <- unique(pred)
    actual.classes <- unique(actual)
    missing.classes <- setdiff(actual.classes, pred.classes)
    tab.rows <- row.names(tab)
    tab <- rbind(tab, matrix(0, nrow = length(missing.classes), 
        ncol = length(tab[1, ])))
    row.names(tab) <- c(tab.rows, missing.classes)
    row.index <- order(row.names(tab), decreasing = F)
    col.index <- order(colnames(tab), decreasing = F)
    tab2 <- tab[row.index, col.index]
    return(tab2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
