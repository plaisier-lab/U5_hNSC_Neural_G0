\name{fold.changes}
\alias{fold.changes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fold.changes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fold.changes(dataset, class.labels, method = "MEAN.DIFF", thres = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{class.labels}{
%%     ~~Describe \code{class.labels} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
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
function (dataset, class.labels, method = "MEAN.DIFF", thres = 1) 
{
    rows <- length(dataset[, 1])
    cols <- length(dataset[1, ])
    fold.changes <- vector(length = rows, mode = "numeric")
    for (i in 1:rows) {
        y <- dataset[i, ]
        y[y < thres] <- thres
        x <- split(y, class.labels)
        if (method == "MEAN.DIFF") {
            m1 <- mean(x[[1]])
            m2 <- mean(x[[2]])
        }
        else if (method == "MEDIAN.DIFF") {
            m1 <- median(x[[1]])
            m2 <- median(x[[2]])
        }
        else {
            stop("ERROR in fold.changes(): unknown method")
        }
        if (abs(m2) > 1e-06) {
            fold.changes[i] <- m1/m2
        }
        else {
            fold.changes[i] <- 1
        }
    }
    return(fold.changes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
