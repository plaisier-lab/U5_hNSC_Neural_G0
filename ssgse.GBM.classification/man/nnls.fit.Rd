\name{nnls.fit}
\alias{nnls.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
nnls.fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nnls.fit(x, y, wsqrt = 1, eps = 0, rank.tol = 1e-07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{wsqrt}{
%%     ~~Describe \code{wsqrt} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{rank.tol}{
%%     ~~Describe \code{rank.tol} here~~
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
function (x, y, wsqrt = 1, eps = 0, rank.tol = 1e-07) 
{
    require("quadprog")
    m <- NCOL(x)
    if (length(eps) == 1) 
        eps <- rep(eps, m)
    x <- x * wsqrt
    y <- y * wsqrt
    xscale <- apply(abs(x), 2, mean)
    yscale <- mean(abs(y))
    x <- t(t(x)/xscale)
    y <- y/yscale
    Rinv <- backsolve(qr.R(qr(x)), diag(m))
    cf <- solve.QP.compact(Dmat = Rinv, dvec = t(x) \%*\% y, Amat = rbind(rep(1, 
        m)), Aind = rbind(rep(1, m), 1:m), bvec = eps * xscale/yscale, 
        factorized = TRUE)$sol
    cf <- cf * yscale/xscale
    cf
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
