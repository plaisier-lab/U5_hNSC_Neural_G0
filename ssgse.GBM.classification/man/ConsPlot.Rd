\name{ConsPlot}
\alias{ConsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ConsPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ConsPlot(V, col.names, main = " ", sub = " ", xlab = " ", ylab = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{col.names}{
%%     ~~Describe \code{col.names} here~~
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
function (V, col.names, main = " ", sub = " ", xlab = " ", ylab = " ") 
{
    cols <- length(V[1, ])
    B <- matrix(0, nrow = cols, ncol = cols)
    max.val <- max(V)
    min.val <- min(V)
    for (i in 1:cols) {
        for (j in 1:cols) {
            k <- cols - i + 1
            B[k, j] <- max.val - V[i, j] + min.val
        }
    }
    col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", 
        "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", 
        "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))
    par(mar = c(5, 20, 20, 1))
    image(1:cols, 1:cols, t(B), col = col.map, axes = FALSE, 
        main = main, sub = sub, xlab = xlab, ylab = ylab)
    for (i in 1:cols) {
        col.names[i] <- substr(col.names[i], 1, 25)
    }
    col.names2 <- rev(col.names)
    size.col.char <- ifelse(cols < 20, 1.25, 1.25 * sqrt(20/cols))
    axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, 
        las = 1, cex.axis = size.col.char, font.axis = 1, line = -1)
    axis(3, at = 1:cols, labels = col.names, adj = 1, tick = FALSE, 
        las = 3, cex.axis = size.col.char, font.axis = 1, line = -1)
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
