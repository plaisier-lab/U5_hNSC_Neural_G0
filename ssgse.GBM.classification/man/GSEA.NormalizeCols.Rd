\name{GSEA.NormalizeCols}
\alias{GSEA.NormalizeCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.NormalizeCols
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.NormalizeCols(V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
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
function (V) 
{
    col.mean <- apply(V, MARGIN = 2, FUN = mean)
    col.sd <- apply(V, MARGIN = 2, FUN = sd)
    col.n <- length(V[1, ])
    for (i in 1:col.n) {
        if (col.sd[i] == 0) {
            V[i, ] <- 0
        }
        else {
            V[, i] <- (V[, i] - col.mean[i])/col.sd[i]
        }
    }
    return(V)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
