\name{CNMF.ConsPlot}
\alias{CNMF.ConsPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CNMF.ConsPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNMF.ConsPlot(V, col.labels, col.names, main = " ", sub = " ", xlab = " ", ylab = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{col.labels}{
%%     ~~Describe \code{col.labels} here~~
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
function (V, col.labels, col.names, main = " ", sub = " ", xlab = " ", 
    ylab = " ") 
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
    col.names2 <- rev(col.names)
    col.labels2 <- rev(col.labels)
    D <- matrix(0, nrow = (cols + 1), ncol = (cols + 1))
    col.tag <- vector(length = cols, mode = "numeric")
    current.tag <- 0
    col.tag[1] <- current.tag
    for (i in 2:cols) {
        if (col.labels[i] != col.labels[i - 1]) {
            current.tag <- 1 - current.tag
        }
        col.tag[i] <- current.tag
    }
    col.tag2 <- rev(col.tag)
    D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag\%\%2 == 0, 1.02, 
        1.01)
    D[1:cols, 1] <- ifelse(col.tag2\%\%2 == 0, 1.02, 1.01)
    D[(cols + 1), 1] <- 1.03
    D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]
    col.map <- c(rainbow(100, s = 1, v = 0.75, start = 0, end = 0.75, 
        gamma = 1.5), "#BBBBBB", "#333333", "#FFFFFF")
    image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes = FALSE, 
        main = main, sub = sub, xlab = xlab, ylab = ylab)
    for (i in 1:cols) {
        col.names[i] <- paste("      ", substr(col.names[i], 
            1, 12), sep = "")
        col.names2[i] <- paste(substr(col.names2[i], 1, 12), 
            "     ", sep = "")
    }
    axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, 
        las = 1, cex.axis = 0.5, font.axis = 1, line = -1)
    axis(2, at = 1:cols, labels = col.labels2, adj = 0.5, tick = FALSE, 
        las = 1, cex.axis = 0.65, font.axis = 1, line = -1)
    axis(3, at = 2:(cols + 1), labels = col.names, adj = 1, tick = FALSE, 
        las = 3, cex.axis = 0.5, font.axis = 1, line = -1)
    axis(3, at = 2:(cols + 1), labels = as.character(col.labels), 
        adj = 1, tick = FALSE, las = 1, cex.axis = 0.65, font.axis = 1, 
        line = -1)
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
