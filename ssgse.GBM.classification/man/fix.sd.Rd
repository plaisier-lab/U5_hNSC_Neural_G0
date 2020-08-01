\name{fix.sd}
\alias{fix.sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fix.sd
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fix.sd(s, m, s.percent = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{s.percent}{
%%     ~~Describe \code{s.percent} here~~
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
function (s, m, s.percent = 0.2) 
{
    min.s <- s.percent * abs(m)
    if (min.s < s) {
        min.s <- s
    }
    if (min.s == 0) {
        min.s <- 0.2
    }
    return(min.s)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
