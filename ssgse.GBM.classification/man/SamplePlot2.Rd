\name{SamplePlot2}
\alias{SamplePlot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SamplePlot2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SamplePlot2(S1 = "NULL", S2 = "NULL", title = " ", phen = "NULL", class.v = "NULL", second.phen = "NULL", second.class.v = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S1}{
%%     ~~Describe \code{S1} here~~
}
  \item{S2}{
%%     ~~Describe \code{S2} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{phen}{
%%     ~~Describe \code{phen} here~~
}
  \item{class.v}{
%%     ~~Describe \code{class.v} here~~
}
  \item{second.phen}{
%%     ~~Describe \code{second.phen} here~~
}
  \item{second.class.v}{
%%     ~~Describe \code{second.class.v} here~~
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
function (S1 = "NULL", S2 = "NULL", title = " ", phen = "NULL", 
    class.v = "NULL", second.phen = "NULL", second.class.v = "NULL") 
{
    num.samples <- length(S1)
    min.S1 <- min(S1)
    max.S1 <- max(S1)
    min.S2 <- min(S2)
    max.S2 <- max(S2)
    limit.S1 <- max(S1) + 0.35 * (max(S1) - min(S1))
    limit.S2 <- max(S2) + 0.35 * (max(S2) - min(S2))
    mid.S1 <- min.S1 + 0.5 * (max.S1 - min.S1)
    plot(S1, S2, ylim = c(min.S2, limit.S2), type = "n", main = title)
    offset <- 20
    for (j in 1:num.samples) {
        color.code <- class.v[j] + 2
        pch.code <- second.class.v[j] + offset
        points(S1[j], S2[j], pch = pch.code, type = "p", cex = 2, 
            bg = color.code, col = 1)
    }
    leg.txt <- phen
    n.phen <- length(phen)
    p.vec <- rep(22, n.phen)
    c.vec <- 2:(n.phen + 1)
    legend(x = min.S1, y = limit.S2, legend = leg.txt, bty = "n", 
        pch = p.vec, bg = "white", pt.bg = c.vec, col = 1, cex = 2)
    leg2.txt <- second.phen
    n.phen2 <- length(second.phen)
    p.vec2 <- 1:n.phen2
    p.vec2 <- p.vec2 + offset
    c.vec2 <- rep(1, n.phen2)
    legend(x = mid.S1, y = limit.S2, legend = leg2.txt, bty = "n", 
        pch = p.vec2, bg = "white", pt.bg = c.vec2, col = c.vec2, 
        cex = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
