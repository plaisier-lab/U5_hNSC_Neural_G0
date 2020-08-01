\name{GSEA.EnrichmentScore3}
\alias{GSEA.EnrichmentScore3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.EnrichmentScore3
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.EnrichmentScore3(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.list}{
%%     ~~Describe \code{gene.list} here~~
}
  \item{gene.set}{
%%     ~~Describe \code{gene.set} here~~
}
  \item{weighted.score.type}{
%%     ~~Describe \code{weighted.score.type} here~~
}
  \item{correl.vector}{
%%     ~~Describe \code{correl.vector} here~~
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
function (gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) 
{
    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    if (weighted.score.type == 0) {
        correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector^alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
    norm.tag <- 1/sum.correl.tag
    norm.no.tag <- 1/Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag - 
        no.tag.indicator * norm.no.tag)
    ES <- sum(RES)
    return(list(ES = ES))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
