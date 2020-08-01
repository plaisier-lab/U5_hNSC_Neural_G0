\name{MSIG.HeatMapPlot.3}
\alias{MSIG.HeatMapPlot.3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.HeatMapPlot.3
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.HeatMapPlot.3(V, row.names = "NA", col.labels = "NA", col.classes = "NA", phen.cmap = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 1, cmap.type = 1, max.v = "NA")
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
  \item{phen.cmap}{
%%     ~~Describe \code{phen.cmap} here~~
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
  \item{char.rescale}{
%%     ~~Describe \code{char.rescale} here~~
}
  \item{cmap.type}{
%%     ~~Describe \code{cmap.type} here~~
}
  \item{max.v}{
%%     ~~Describe \code{max.v} here~~
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
    phen.cmap = "NA", col.names = "NA", main = " ", sub = " ", 
    xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 1, 
    cmap.type = 1, max.v = "NA") 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    if ((cmap.type == 3) | (cmap.type == 5)) {
        row.norm <- F
    }
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
        mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", 
            "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", 
            "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D")
    }
    else if (cmap.type == 2) {
        mycol <- c("#FCFBFD", "#F4F2F8", "#F8F7FB", "#EFEDF5", 
            "#E1E1EF", "#E8E7F2", "#DADAEB", "#C6C7E1", "#D0D1E6", 
            "#BCBDDC", "#A8A6CF", "#B2B2D6", "#9E9AC8", "#8A87BF", 
            "#9491C4", "#807DBA", "#7260AB", "#796FB3", "#6A51A3", 
            "#5C3596", "#63439D", "#54278F", "#460D83", "#4D1A89", 
            "#3F007D")
    }
    else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 
        5)) {
        mycol <- vector(length = 512, mode = "numeric")
        for (k in 1:256) {
            mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
        }
        for (k in 257:512) {
            mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, 
                maxColorValue = 255)
        }
        mycol <- rev(mycol)
    }
    else if (cmap.type == 6) {
        mycol <- c("#4500AD", "#2700D1", "#6B58EF", "#8888FF", 
            "#C7C1FF", "#D5D5FF", "#FFC0E5", "#FF8989", "#FF7080", 
            "#FF5A5A", "#EF4040", "#D60C00")
    }
    ncolors <- length(mycol)
    mycol <- c(mycol, phen.cmap[1:length(col.classes)])
    if (cmap.type == 5) {
        if (max.v == "NA") {
            max.v <- max(max(V), -min(V))
        }
        V <- ceiling(ncolors * (V - (-max.v))/(1.001 * (max.v - 
            (-max.v))))
    }
    else {
        V <- ceiling(ncolors * (V - min(V))/(1.001 * (max(V) - 
            min(V))))
    }
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
    heatm[n.rows + 1, ] <- ncolors + col.labels
    if (cmap.type == 2) {
        par(mar = c(3, 7, 3, 1))
    }
    else {
        par(mar = c(4, 15, 4, 1))
    }
    print(c("range=", range(V)))
    if (cmap.type == 5) {
        image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, 
            ncolors + max(col.labels)), col = mycol, axes = FALSE, 
            main = main, sub = sub, xlab = xlab, ylab = ylab)
    }
    else {
        image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, 
            axes = FALSE, main = main, sub = sub, xlab = xlab, 
            ylab = ylab)
    }
    if (row.names[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 30/(n.rows + 15)
        size.col.char <- char.rescale * 30/(n.cols + 15)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 30)
        }
        row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
        axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (col.names[1] != "NA") {
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
