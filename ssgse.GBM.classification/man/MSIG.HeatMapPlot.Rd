\name{MSIG.HeatMapPlot}
\alias{MSIG.HeatMapPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.HeatMapPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.HeatMapPlot(V, row.names = "NA", col.labels = "NA", col.classes = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", row.norm = TRUE, cmap.type = 1)
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
  \item{sub}{
%%     ~~Describe \code{sub} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{row.norm}{
%%     ~~Describe \code{row.norm} here~~
}
  \item{cmap.type}{
%%     ~~Describe \code{cmap.type} here~~
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
function (V, row.names = "NA", col.labels = "NA", col.classes = "NA", 
    col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", 
    row.norm = TRUE, cmap.type = 1) 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    if (row.norm == TRUE) {
        row.mean <- apply(V, MARGIN = 1, FUN = mean)
        row.sd <- apply(V, MARGIN = 1, FUN = sd)
        row.n <- length(V[, 1])
        for (i in 1:n.rows) {
            if (row.sd[i] == 0) {
                V[i, ] <- 0
            }
            else {
                V[i, ] <- (V[i, ] - row.mean[i])/(0.333 * row.sd[i])
            }
            V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
            V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
        }
    }
    if (cmap.type == 1) {
        mycol <- c("#BBBBBB", "#0000FF", "#4040FF", "#7070FF", 
            "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", 
            "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", 
            "#333333")
        ncolors <- length(mycol) - 2
    }
    else if (cmap.type >= 2) {
        library("RColorBrewer")
        mycol <- c("#FFA500", "#FCFBFD", "#F4F2F8", "#F8F7FB", 
            "#EFEDF5", "#E1E1EF", "#E8E7F2", "#DADAEB", "#C6C7E1", 
            "#D0D1E6", "#BCBDDC", "#A8A6CF", "#B2B2D6", "#9E9AC8", 
            "#8A87BF", "#9491C4", "#807DBA", "#7260AB", "#796FB3", 
            "#6A51A3", "#5C3596", "#63439D", "#54278F", "#460D83", 
            "#4D1A89", "#3F007D", "#2E8657")
        ncolors <- length(mycol) - 2
    }
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
    maxv <- max(V)
    minv <- min(V)
    rangev <- maxv - minv
    heatm[n.rows + 1, ] <- ifelse(col.labels\%\%2 == 0, maxv + 
        (rangev/(ncolors - 1)), minv - (rangev/(ncolors - 1)))
    par(mar = c(10, 20, 10, 4))
    image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, 
        main = main, sub = sub, xlab = xlab, ylab = ylab)
    if (length(row.names) > 1) {
        numC <- nchar(row.names)
        size.row.char <- 20/(n.rows + 15)
        size.col.char <- 35/(n.cols + 15)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 35)
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
    if (length(col.labels) > 1) {
        tab <- table(col.labels)
        zinit <- cumsum(tab) - tab + 1
        zfinal <- cumsum(tab)
        locations <- ceiling((zinit + zfinal)/2)
        size.class.char <- 13/(length(tab) + 5)
        axis(3, at = locations, labels = col.classes, tick = FALSE, 
            las = 1, cex.axis = size.class.char, font.axis = 2, 
            line = -1)
    }
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
