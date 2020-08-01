\name{Wilcox.Score}
\alias{Wilcox.Score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wilcox.Score
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Wilcox.Score(gene.list, gene.set)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.list}{
%%     ~~Describe \code{gene.list} here~~
}
  \item{gene.set}{
%%     ~~Describe \code{gene.set} here~~
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
function (gene.list, gene.set) 
{
    library(exactRankTests)
    N <- length(gene.list)
    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
    seq.index <- seq(1, N)
    gene.set.ranks <- seq.index[tag.indicator == 1]
    gene.set.comp.ranks <- seq.index[tag.indicator == 0]
    W <- wilcox.exact(x = gene.set.ranks, y = gene.set.comp.ranks, 
        alternative = "two.sided", mu = 0, paired = FALSE, exact = F, 
        conf.int = T, conf.level = 0.95)
    ES <- log(1/W$p.value)
    return(list(ES = ES))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
