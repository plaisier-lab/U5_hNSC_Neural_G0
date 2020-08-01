\name{GSEA.EnrichmentScore2}
\alias{GSEA.EnrichmentScore2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.EnrichmentScore2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.EnrichmentScore2(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL)
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
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    loc.vector <- vector(length = N, mode = "numeric")
    peak.res.vector <- vector(length = Nh, mode = "numeric")
    valley.res.vector <- vector(length = Nh, mode = "numeric")
    tag.correl.vector <- vector(length = Nh, mode = "numeric")
    tag.diff.vector <- vector(length = Nh, mode = "numeric")
    tag.loc.vector <- vector(length = Nh, mode = "numeric")
    loc.vector[gene.list] <- seq(1, N)
    tag.loc.vector <- loc.vector[gene.set]
    tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
    if (weighted.score.type == 0) {
        tag.correl.vector <- rep(1, Nh)
    }
    else if (weighted.score.type == 1) {
        tag.correl.vector <- correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    }
    else if (weighted.score.type == 2) {
        tag.correl.vector <- correl.vector[tag.loc.vector] * 
            correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    }
    else {
        tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
        tag.correl.vector <- abs(tag.correl.vector)
    }
    norm.tag <- 1/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 
        1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- ifelse(max.ES > -min.ES, max.ES, min.ES)
    return(list(ES = ES))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
