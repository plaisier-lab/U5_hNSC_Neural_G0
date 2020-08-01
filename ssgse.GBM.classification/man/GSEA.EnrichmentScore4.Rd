\name{GSEA.EnrichmentScore4}
\alias{GSEA.EnrichmentScore4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.EnrichmentScore4
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.EnrichmentScore4(gene.list, gene.set, statistic = "Kolmogorov-Smirnov", alpha = 1, correl.vector = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.list}{
%%     ~~Describe \code{gene.list} here~~
}
  \item{gene.set}{
%%     ~~Describe \code{gene.set} here~~
}
  \item{statistic}{
%%     ~~Describe \code{statistic} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (gene.list, gene.set, statistic = "Kolmogorov-Smirnov", 
    alpha = 1, correl.vector = NULL) 
{
    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    orig.correl.vector <- correl.vector
    if (alpha == 0) 
        correl.vector <- rep(1, N)
    correl.vector <- abs(correl.vector)^alpha
    sum.correl <- sum(correl.vector[tag.indicator == 1])
    P0 <- no.tag.indicator/Nm
    F0 <- cumsum(P0)
    Pn <- tag.indicator * correl.vector/sum.correl
    Fn <- cumsum(Pn)
    if (statistic == "Kolmogorov-Smirnov") {
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
            ES <- signif(max.ES, digits = 5)
            arg.ES <- which.max(RES)
        }
        else {
            ES <- signif(min.ES, digits = 5)
            arg.ES <- which.min(RES)
        }
        return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
    }
    else if (statistic == "Cramer-von-Mises") {
        RES <- Fn - F0
        X <- RES^2 * P0
        X_p <- X[RES >= 0]
        X_n <- X[RES < 0]
        ES_p <- sqrt(sum(X_p)/N)
        ES_n <- sqrt(sum(X_n)/N)
        if (ES_p > ES_n) {
            ES <- signif(ES_p, digits = 5)
            arg.ES <- which.min(abs(X - max(X_p)))
        }
        else {
            ES <- -signif(ES_n, digits = 5)
            arg.ES <- which.min(abs(X - max(X_n)))
        }
        return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
    }
    else if (statistic == "Anderson-Darling") {
        RES <- Fn - F0
        F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1, 
            N), F0 * (1 - F0))
        X <- RES^2 * P0/F0_factor
        X_p <- X[RES >= 0]
        X_n <- X[RES < 0]
        ES_p <- sqrt(sum(X_p)/N)
        ES_n <- sqrt(sum(X_n)/N)
        if (ES_p > ES_n) {
            ES <- signif(ES_p, digits = 5)
            arg.ES <- which.min(abs(X - max(X_p)))
        }
        else {
            ES <- -signif(ES_n, digits = 5)
            arg.ES <- which.min(abs(X - max(X_n)))
        }
        return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
    }
    else if (statistic == "Zhang_A") {
        RES <- Fn - F0
        Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * 
            log(Fn/F0))
        Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 
            1)/sum.correl, 0, (1 - Fn) * log((1 - Fn)/(1 - F0)))
        Fn_factor <- ifelse(Fn < 1/sum.correl | Fn > (sum.correl - 
            1)/sum.correl, rep(1, N), Fn * (1 - Fn))
        G <- (Fact1 + Fact2) * Pn/Fn_factor
        G_p <- G[RES >= 0]
        G_n <- G[RES < 0]
        ES_p <- sum(G_p)/N
        ES_n <- sum(G_n)/N
        if (ES_p > ES_n) {
            ES <- signif(ES_p, digits = 5)
            arg.ES <- which.min(abs(G - max(G_p)))
        }
        else {
            ES <- -signif(ES_n, digits = 5)
            arg.ES <- which.min(abs(G - max(G_n)))
        }
        return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
    }
    else if (statistic == "Zhang_C") {
        RES <- Fn - F0
        Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * 
            log(Fn/F0))
        Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 
            1)/sum.correl, 0, (1 - Fn) * log((1 - Fn)/(1 - F0)))
        F0_factor <- ifelse(F0 < 1/Nm | F0 > (Nm - 1)/Nm, rep(1, 
            N), F0 * (1 - F0))
        G <- (Fact1 + Fact2) * P0/F0_factor
        G_p <- G[RES >= 0]
        G_n <- G[RES < 0]
        ES_p <- sum(G_p)/N
        ES_n <- sum(G_n)/N
        if (ES_p > ES_n) {
            ES <- signif(ES_p, digits = 5)
            arg.ES <- which.min(abs(G - max(G_p)))
        }
        else {
            ES <- -signif(ES_n, digits = 5)
            arg.ES <- which.min(abs(G - max(G_n)))
        }
        return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
    }
    else if (statistic == "Zhang_K") {
        RES <- Fn - F0
        Fact1 <- ifelse(F0 < 1/Nm | Fn < 1/sum.correl, 0, Fn * 
            log(Fn/F0))
        Fact2 <- ifelse(F0 > (Nm - 1)/Nm | Fn > (sum.correl - 
            1)/sum.correl, 0, (1 - Fn) * log((1 - Fn)/(1 - F0)))
        G <- Fact1 + Fact2
        G_p <- G[RES >= 0]
        G_n <- G[RES < 0]
        ES_p <- max(G_p)
        ES_n <- max(G_n)
        if (ES_p > ES_n) {
            ES <- signif(ES_p, digits = 5)
            arg.ES <- which.min(abs(G - ES_p))
        }
        else {
            ES <- -signif(ES_n, digits = 5)
            arg.ES <- which.min(abs(G - ES_n))
        }
        return(list(ES = ES, RES = RES, arg.ES = arg.ES, indicator = tag.indicator))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
