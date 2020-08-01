\name{GSEA.GeneRanking}
\alias{GSEA.GeneRanking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.GeneRanking
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.GeneRanking(A, class.labels, gene.labels, nperm, permutation.type = 0, sigma.correction = "GeneCluster", fraction = 1, replace = F, reverse.sign = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{class.labels}{
%%     ~~Describe \code{class.labels} here~~
}
  \item{gene.labels}{
%%     ~~Describe \code{gene.labels} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{permutation.type}{
%%     ~~Describe \code{permutation.type} here~~
}
  \item{sigma.correction}{
%%     ~~Describe \code{sigma.correction} here~~
}
  \item{fraction}{
%%     ~~Describe \code{fraction} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
}
  \item{reverse.sign}{
%%     ~~Describe \code{reverse.sign} here~~
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
function (A, class.labels, gene.labels, nperm, permutation.type = 0, 
    sigma.correction = "GeneCluster", fraction = 1, replace = F, 
    reverse.sign = F) 
{
    A <- A + 1e-08
    N <- length(A[, 1])
    Ns <- length(A[1, ])
    subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
    reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
    reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
    class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
    class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
    order.matrix <- matrix(0, nrow = N, ncol = nperm)
    obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
    s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
    obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
    obs.gene.labels <- vector(length = N, mode = "character")
    obs.gene.descs <- vector(length = N, mode = "character")
    obs.gene.symbols <- vector(length = N, mode = "character")
    M1 <- matrix(0, nrow = N, ncol = nperm)
    M2 <- matrix(0, nrow = N, ncol = nperm)
    S1 <- matrix(0, nrow = N, ncol = nperm)
    S2 <- matrix(0, nrow = N, ncol = nperm)
    gc()
    C <- split(class.labels, class.labels)
    class1.size <- length(C[[1]])
    class2.size <- length(C[[2]])
    class1.index <- seq(1, class1.size, 1)
    class2.index <- seq(class1.size + 1, class1.size + class2.size, 
        1)
    for (r in 1:nperm) {
        class1.subset <- sample(class1.index, size = ceiling(class1.size * 
            fraction), replace = replace)
        class2.subset <- sample(class2.index, size = ceiling(class2.size * 
            fraction), replace = replace)
        class1.subset.size <- length(class1.subset)
        class2.subset.size <- length(class2.subset)
        subset.class1 <- rep(0, class1.size)
        for (i in 1:class1.size) {
            if (is.element(class1.index[i], class1.subset)) {
                subset.class1[i] <- 1
            }
        }
        subset.class2 <- rep(0, class2.size)
        for (i in 1:class2.size) {
            if (is.element(class2.index[i], class2.subset)) {
                subset.class2[i] <- 1
            }
        }
        subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
        fraction.class1 <- class1.size/Ns
        fraction.class2 <- class2.size/Ns
        if (permutation.type == 0) {
            full.subset <- c(class1.subset, class2.subset)
            label1.subset <- sample(full.subset, size = Ns * 
                fraction.class1)
            reshuffled.class.labels1[, r] <- rep(0, Ns)
            reshuffled.class.labels2[, r] <- rep(0, Ns)
            class.labels1[, r] <- rep(0, Ns)
            class.labels2[, r] <- rep(0, Ns)
            for (i in 1:Ns) {
                m1 <- sum(!is.na(match(label1.subset, i)))
                m2 <- sum(!is.na(match(full.subset, i)))
                reshuffled.class.labels1[i, r] <- m1
                reshuffled.class.labels2[i, r] <- m2 - m1
                if (i <= class1.size) {
                  class.labels1[i, r] <- m2
                  class.labels2[i, r] <- 0
                }
                else {
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
                }
            }
        }
        else if (permutation.type == 1) {
            class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * 
                fraction.class1))
            class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * 
                fraction.class1))
            reshuffled.class.labels1[, r] <- rep(0, Ns)
            reshuffled.class.labels2[, r] <- rep(0, Ns)
            class.labels1[, r] <- rep(0, Ns)
            class.labels2[, r] <- rep(0, Ns)
            for (i in 1:Ns) {
                if (i <= class1.size) {
                  m1 <- sum(!is.na(match(class1.label1.subset, 
                    i)))
                  m2 <- sum(!is.na(match(class1.subset, i)))
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- m2
                  class.labels2[i, r] <- 0
                }
                else {
                  m1 <- sum(!is.na(match(class2.label1.subset, 
                    i)))
                  m2 <- sum(!is.na(match(class2.subset, i)))
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
                }
            }
        }
    }
    P <- reshuffled.class.labels1 * subset.mask
    n1 <- sum(P[, 1])
    M1 <- A \%*\% P
    M1 <- M1/n1
    gc()
    A2 <- A * A
    S1 <- A2 \%*\% P
    S1 <- S1/n1 - M1 * M1
    S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
    gc()
    P <- reshuffled.class.labels2 * subset.mask
    n2 <- sum(P[, 1])
    M2 <- A \%*\% P
    M2 <- M2/n2
    gc()
    A2 <- A * A
    S2 <- A2 \%*\% P
    S2 <- S2/n2 - M2 * M2
    S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
    rm(P)
    rm(A2)
    gc()
    if (sigma.correction == "GeneCluster") {
        S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        gc()
    }
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- S1 + S2
    rm(S2)
    gc()
    s2n.matrix <- M1/S1
    if (reverse.sign == T) {
        s2n.matrix <- -s2n.matrix
    }
    gc()
    P <- class.labels1 * subset.mask
    n1 <- sum(P[, 1])
    M1 <- A \%*\% P
    M1 <- M1/n1
    gc()
    A2 <- A * A
    S1 <- A2 \%*\% P
    S1 <- S1/n1 - M1 * M1
    S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
    gc()
    P <- class.labels2 * subset.mask
    n2 <- sum(P[, 1])
    M2 <- A \%*\% P
    M2 <- M2/n2
    gc()
    A2 <- A * A
    S2 <- A2 \%*\% P
    S2 <- S2/n2 - M2 * M2
    S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
    rm(P)
    rm(A2)
    gc()
    if (sigma.correction == "GeneCluster") {
        S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        gc()
    }
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- S1 + S2
    rm(S2)
    gc()
    obs.s2n.matrix <- M1/S1
    gc()
    if (reverse.sign == T) {
        obs.s2n.matrix <- -obs.s2n.matrix
    }
    for (r in 1:nperm) {
        obs.order.matrix[, r] <- order(obs.s2n.matrix[, r], decreasing = T)
    }
    return(list(s2n.matrix = s2n.matrix, obs.s2n.matrix = obs.s2n.matrix, 
        order.matrix = order.matrix, obs.order.matrix = obs.order.matrix))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
