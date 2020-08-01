\name{GSEA.HeatMapPlot2}
\alias{GSEA.HeatMapPlot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.HeatMapPlot2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.HeatMapPlot2(V, row.names = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{row.names}{
%%     ~~Describe \code{row.names} here~~
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
function (V, row.names = "NA", col.names = "NA", main = " ", 
    sub = " ", xlab = " ", ylab = " ") 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    mycol <- rev(rainbow(100, s = 1, v = 0.75, start = 0, end = 0.75, 
        gamma = 1.5))
    heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
    heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
    par(mar = c(5, 20, 5, 2))
    image(1:n.cols, 1:n.rows, t(heatm), col = mycol, axes = FALSE, 
        main = main, sub = sub, xlab = xlab, ylab = ylab)
    if (length(row.names) > 1) {
        size.row.char <- ifelse(n.rows < 15, 1, sqrt(15/n.rows))
        size.col.char <- ifelse(n.cols < 25, 1, sqrt(25/n.cols))
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 40)
        }
        row.names <- row.names[seq(n.rows, 1, -1)]
        axis(2, at = 1:n.rows, labels = row.names, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 1, line = -1)
    }
    if (length(col.names) > 1) {
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, 
            las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
