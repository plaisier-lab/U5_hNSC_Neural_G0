\name{GSEA.HeatMapPlot}
\alias{GSEA.HeatMapPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA.HeatMapPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA.HeatMapPlot(V, row.names = F, col.labels, col.classes, col.names = F, main = " ", xlab = " ", ylab = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{row.names}{
%%     ~~Describe \code{row.names} here~~
}
  \item{col.labels}{
%%     ~~Describe \code{col.labels} here~~
}
  \item{col.classes}{
%%     ~~Describe \code{col.classes} here~~
}
  \item{col.names}{
%%     ~~Describe \code{col.names} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
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
function (V, row.names = F, col.labels, col.classes, col.names = F, 
    main = " ", xlab = " ", ylab = " ") 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
        if (row.sd[i] == 0) {
            V[i, ] <- 0
        }
        else {
            V[i, ] <- (V[i, ] - row.mean[i])/(0.5 * row.sd[i])
        }
        V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
        V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
    }
    mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", 
        "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", 
        "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")
    mid.range.V <- mean(range(V)) - 0.1
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
    heatm[n.rows + 1, ] <- ifelse(col.labels == 0, 7, -7)
    image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, 
        main = main, xlab = xlab, ylab = ylab)
    if (length(row.names) > 1) {
        numC <- nchar(row.names)
        size.row.char <- 35/(n.rows + 5)
        size.col.char <- 25/(n.cols + 5)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 10)
        }
        row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
        axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (length(col.names) > 1) {
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, 
            las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    C <- split(col.labels, col.labels)
    class1.size <- length(C[[1]])
    class2.size <- length(C[[2]])
    axis(3, at = c(floor(class1.size/2), class1.size + floor(class2.size/2)), 
        labels = col.classes, tick = FALSE, las = 1, cex.axis = 1.25, 
        font.axis = 2, line = -1)
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
