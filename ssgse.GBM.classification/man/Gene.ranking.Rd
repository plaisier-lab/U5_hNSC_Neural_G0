\name{Gene.ranking}
\alias{Gene.ranking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene.ranking
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Gene.ranking(dataset, class.labels, method = "S2N")
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
function (dataset, class.labels, method = "S2N") 
{
    rows <- length(dataset[, 1])
    cols <- length(dataset[1, ])
    values.vector <- vector(length = rows, mode = "numeric")
    if (method == "S2N") {
        for (i in 1:rows) {
            values.vector[i] <- S2N(dataset[i, ], class.labels)
        }
    }
    else if (method == "RS2N") {
        for (i in 1:rows) {
            values.vector[i] <- RS2N(dataset[i, ], class.labels)
        }
    }
    else if (method == "MEAN.DIFF") {
        for (i in 1:rows) {
            values.vector[i] <- MEAN.DIFF(dataset[i, ], class.labels)
        }
    }
    else if (method == "MEDIAN.DIFF") {
        for (i in 1:rows) {
            values.vector[i] <- MEDIAN.DIFF(dataset[i, ], class.labels)
        }
    }
    else if (method == "ROC") {
        status <- ifelse(class.labels == class.labels[1], 1, 
            0)
        for (i in 1:rows) {
            m.score <- as.numeric(dataset[i, ])
            m.score.norm <- (m.score - min(m.score))/(max(m.score) - 
                min(m.score))
            values.vector[i] <- signif(roc.area(status, m.score.norm)$A, 
                digits = 3)
        }
    }
    else {
        stop(c("Gene.ranking -- unknown method:", method))
    }
    return(values.vector)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
