\name{metagene.plot}
\alias{metagene.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
metagene.plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metagene.plot(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{
%%     ~~Describe \code{H} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
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
function (H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") 
{
    k <- length(H[, 1])
    S <- length(H[1, ])
    index <- 1:S
    maxval <- max(H)
    minval <- min(H)
    plot(index, H[1, ], xlim = c(1, S), ylim = c(minval, maxval), 
        main = main, sub = sub, ylab = ylab, xlab = xlab, type = "n")
    for (i in 1:k) {
        lines(index, H[i, ], type = "l", col = i, lwd = 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
