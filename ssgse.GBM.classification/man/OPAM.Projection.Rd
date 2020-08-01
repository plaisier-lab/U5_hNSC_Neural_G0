\name{OPAM.Projection}
\alias{OPAM.Projection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.Projection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.Projection(data.array, gene.names, n.cols, n.rows, weight = 0, statistic = "Kolmogorov-Smirnov", gene.set, nperm = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.array}{
%%     ~~Describe \code{data.array} here~~
}
  \item{gene.names}{
%%     ~~Describe \code{gene.names} here~~
}
  \item{n.cols}{
%%     ~~Describe \code{n.cols} here~~
}
  \item{n.rows}{
%%     ~~Describe \code{n.rows} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{statistic}{
%%     ~~Describe \code{statistic} here~~
}
  \item{gene.set}{
%%     ~~Describe \code{gene.set} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
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
function (data.array, gene.names, n.cols, n.rows, weight = 0, 
    statistic = "Kolmogorov-Smirnov", gene.set, nperm = 200) 
{
    ES.vector <- vector(length = n.cols)
    NES.vector <- vector(length = n.cols)
    p.val.vector <- vector(length = n.cols)
    correl.vector <- vector(length = n.rows, mode = "numeric")
    phi <- array(0, c(n.cols, nperm))
    for (sample.index in 1:n.cols) {
        gene.list <- order(data.array[, sample.index], decreasing = T)
        gene.set2 <- match(gene.set, gene.names)
        if (weight == 0) {
            correl.vector <- rep(1, n.rows)
        }
        else if (weight > 0) {
            correl.vector <- data.array[gene.list, sample.index]
        }
        GSEA.results <- GSEA.EnrichmentScore5(gene.list = gene.list, 
            gene.set = gene.set2, statistic = statistic, alpha = weight, 
            correl.vector = correl.vector)
        ES.vector[sample.index] <- GSEA.results$ES
        if (nperm == 0) {
            NES.vector[sample.index] <- ES.vector[sample.index]
            p.val.vector[sample.index] <- 1
        }
        else {
            for (r in 1:nperm) {
                reshuffled.gene.labels <- sample(1:n.rows)
                if (weight == 0) {
                  correl.vector <- rep(1, n.rows)
                }
                else if (weight > 0) {
                  correl.vector <- data.array[reshuffled.gene.labels, 
                    sample.index]
                }
                GSEA.results <- GSEA.EnrichmentScore5(gene.list = reshuffled.gene.labels, 
                  gene.set = gene.set2, statistic = statistic, 
                  alpha = weight, correl.vector = correl.vector)
                phi[sample.index, r] <- GSEA.results$ES
            }
            if (ES.vector[sample.index] >= 0) {
                pos.phi <- phi[sample.index, phi[sample.index, 
                  ] >= 0]
                if (length(pos.phi) == 0) 
                  pos.phi <- 0.5
                pos.m <- mean(pos.phi)
                NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
                s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
                p.val.vector[sample.index] <- ifelse(s == 0, 
                  1/nperm, s)
            }
            else {
                neg.phi <- phi[sample.index, phi[sample.index, 
                  ] < 0]
                if (length(neg.phi) == 0) 
                  neg.phi <- 0.5
                neg.m <- mean(neg.phi)
                NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
                s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
                p.val.vector[sample.index] <- ifelse(s == 0, 
                  1/nperm, s)
            }
        }
    }
    return(list(ES.vector = ES.vector, NES.vector = NES.vector, 
        p.val.vector = p.val.vector))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
