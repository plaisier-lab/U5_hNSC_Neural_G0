\name{S2N_GC}
\alias{S2N_GC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S2N_GC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
S2N_GC(A, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{C}{
%%     ~~Describe \code{C} here~~
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
function (A, C) 
{
    x <- split(A, C)
    mean.val <- sapply(x, mean)
    std.val <- sapply(x, sd)
    s2n <- (mean.val[1] - mean.val[2])/(fix.sd(std.val[1], mean.val[1]) + 
        fix.sd(std.val[2], mean.val[2]))
    return(s2n)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
