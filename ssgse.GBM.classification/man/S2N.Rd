\name{S2N}
\alias{S2N}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S2N
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
S2N(A, C)
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
    m1 <- mean(x[[1]])
    m2 <- mean(x[[2]])
    s1 <- ifelse(length(x[[1]]) > 1, sd(x[[1]]), 0)
    s2 <- ifelse(length(x[[2]]) > 1, sd(x[[2]]), 0)
    s1 <- ifelse(s1 < 0.1 * abs(m1), 0.1 * abs(m1), s1)
    s2 <- ifelse(s2 < 0.1 * abs(m2), 0.1 * abs(m2), s2)
    s2n <- (m1 - m2)/(s1 + s2 + 0.1)
    return(s2n)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
