\name{CNMF.matrix.abs.plot}
\alias{CNMF.matrix.abs.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CNMF.matrix.abs.plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNMF.matrix.abs.plot(V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{axes}{
%%     ~~Describe \code{axes} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
}
  \item{transpose}{
%%     ~~Describe \code{transpose} here~~
}
  \item{matrix.order}{
%%     ~~Describe \code{matrix.order} here~~
}
  \item{max.v}{
%%     ~~Describe \code{max.v} here~~
}
  \item{min.v}{
%%     ~~Describe \code{min.v} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{sub}{
%%     ~~Describe \code{sub} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
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
function (V, axes = F, log = F, norm = T, transpose = T, matrix.order = T, 
    max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", 
    ylab = "  ") 
{
    rows <- length(V[, 1])
    cols <- length(V[1, ])
    if (log == T) {
        V <- log(V)
    }
    B <- matrix(0, nrow = rows, ncol = cols)
    for (i in 1:rows) {
        for (j in 1:cols) {
            if (matrix.order == T) {
                k <- rows - i + 1
            }
            else {
                k <- i
            }
            if (norm == T) {
                if ((max.v == 1) && (min.v == 0)) {
                  max.val <- max(V)
                  min.val <- min(V)
                }
                else {
                  max.val = max.v
                  min.val = min.v
                }
            }
            B[k, j] <- max.val - V[i, j] + min.val
        }
    }
    if (transpose == T) {
        B <- t(B)
    }
    if (norm == T) {
        image(z = B, zlim = c(min.val, max.val), axes = axes, 
            col = rainbow(100, s = 1, v = 0.75, start = 0, end = 0.75, 
                gamma = 1.5), main = main, sub = sub, xlab = xlab, 
            ylab = ylab)
    }
    else {
        image(z = B, axes = axes, col = rainbow(100, s = 1, v = 0.6, 
            start = 0.1, end = 0.9, gamma = 1), main = main, 
            sub = sub, xlab = xlab, ylab = ylab)
    }
    return(list(B, max.val, min.val))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
