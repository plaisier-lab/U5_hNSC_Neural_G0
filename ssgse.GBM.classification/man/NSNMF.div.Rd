\name{NSNMF.div}
\alias{NSNMF.div}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NSNMF.div
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NSNMF.div(V, k, theta, maxniter = 2000, seed = 123456, stopconv = 40, stopfreq = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{maxniter}{
%%     ~~Describe \code{maxniter} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{stopconv}{
%%     ~~Describe \code{stopconv} here~~
}
  \item{stopfreq}{
%%     ~~Describe \code{stopfreq} here~~
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
function (V, k, theta, maxniter = 2000, seed = 123456, stopconv = 40, 
    stopfreq = 10) 
{
    N <- length(V[, 1])
    M <- length(V[1, ])
    set.seed(seed)
    W <- matrix(runif(N * k), nrow = N, ncol = k)
    H <- matrix(runif(k * M), nrow = k, ncol = M)
    VP <- matrix(nrow = N, ncol = M)
    error.v <- vector(mode = "numeric", length = maxniter)
    new.membership <- vector(mode = "numeric", length = M)
    old.membership <- vector(mode = "numeric", length = M)
    no.change.count <- 0
    eps <- .Machine$double.eps
    SS <- theta/k * (rep(1, k) \%*\% t(rep(1, k)))
    for (i in 1:k) {
        SS[i, i] <- SS[i, i] + (1 - theta)
    }
    for (t in 1:maxniter) {
        W <- W \%*\% SS
        VP = W \%*\% H
        W.t <- t(W)
        H <- H * (W.t \%*\% (V/VP)) + eps
        norm <- apply(W, MARGIN = 2, FUN = sum)
        for (i in 1:k) {
            H[i, ] <- H[i, ]/norm[i]
        }
        H <- SS \%*\% H
        VP = W \%*\% H
        H.t <- t(H)
        W <- W * ((V/VP) \%*\% H.t) + eps
        norm <- apply(H, MARGIN = 1, FUN = sum)
        for (i in 1:k) {
            W[, i] <- W[, i]/norm[i]
        }
        error.v[t] <- sum(V * log((V + eps)/(VP + eps)) - V + 
            VP)/(M * N)
        if (t\%\%stopfreq == 0) {
            for (j in 1:M) {
                class <- order(H[, j], decreasing = T)
                new.membership[j] <- class[1]
            }
            if (sum(new.membership == old.membership) == M) {
                no.change.count <- no.change.count + 1
            }
            else {
                no.change.count <- 0
            }
            if (no.change.count == stopconv) 
                break
            old.membership <- new.membership
        }
    }
    return(list(W = W, H = H, t = t, error.v = error.v))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
