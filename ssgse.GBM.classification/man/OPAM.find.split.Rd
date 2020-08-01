\name{OPAM.find.split}
\alias{OPAM.find.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.find.split
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.find.split(msig, status)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msig}{
%%     ~~Describe \code{msig} here~~
}
  \item{status}{
%%     ~~Describe \code{status} here~~
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
function (msig, status) 
{
    index <- order(msig, decreasing = F)
    msig.s <- msig[index]
    status.s <- status[index]
    error.1 <- cumsum(status.s == 1)
    error.0 <- rev(cumsum(rev(status.s) == 0))
    min.i <- which.min(abs(error.1 - error.0))
    thres <- msig.s[min.i]
    return(thres)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
