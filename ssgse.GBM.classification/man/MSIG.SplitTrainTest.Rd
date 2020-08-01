\name{MSIG.SplitTrainTest}
\alias{MSIG.SplitTrainTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.SplitTrainTest
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.SplitTrainTest(V, class.v = 0, fraction = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{class.v}{
%%     ~~Describe \code{class.v} here~~
}
  \item{fraction}{
%%     ~~Describe \code{fraction} here~~
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
function (V, class.v = 0, fraction = 0.5) 
{
    cols <- length(V[1, ])
    rows <- length(V[, 1])
    train.size <- ceiling(cols * fraction)
    test.size <- cols - train.size
    index <- sample(1:cols)
    train <- matrix(0, nrow = rows, ncol = train.size)
    test <- matrix(0, nrow = rows, ncol = test.size)
    for (i in 1:train.size) {
        train[, i] <- V[, index[i]]
    }
    for (i in 1:test.size) {
        test[, i] <- V[, index[train.size + i]]
    }
    if (length(class.v) > 0) {
        cls.train <- vector(mode = typeof(class.v), length = train.size)
        cls.test <- vector(mode = typeof(class.v), length = test.size)
        for (i in 1:train.size) {
            cls.train[i] <- class.v[index[i]]
        }
        for (i in 1:test.size) {
            cls.test[i] <- class.v[index[train.size + i]]
        }
        return(list(train = train, test = test, cls.train = cls.train, 
            cls.test = cls.test))
    }
    else {
        return(list(train = train, test = test))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
