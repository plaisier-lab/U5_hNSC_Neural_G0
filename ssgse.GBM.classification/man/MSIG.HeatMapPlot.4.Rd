\name{MSIG.HeatMapPlot.4}
\alias{MSIG.HeatMapPlot.4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.HeatMapPlot.4
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.HeatMapPlot.4(V, row.names = "NA", row.names2 = "NA", col.labels = "NA", col.labels2 = "NA", col.classes = "NA", phen.cmap = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 1, cmap.type = 1, max.v = "NA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
%%     ~~Describe \code{V} here~~
}
  \item{row.names}{
%%     ~~Describe \code{row.names} here~~
}
  \item{row.names2}{
%%     ~~Describe \code{row.names2} here~~
}
  \item{col.labels}{
%%     ~~Describe \code{col.labels} here~~
}
  \item{col.labels2}{
%%     ~~Describe \code{col.labels2} here~~
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
function (V, row.names = "NA", row.names2 = "NA", col.labels = "NA", 
    col.labels2 = "NA", col.classes = "NA", phen.cmap = "NA", 
    col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", 
    row.norm = TRUE, char.rescale = 1, cmap.type = 1, max.v = "NA") 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    V1 <- matrix(0, nrow = n.rows, ncol = n.cols)
    if ((cmap.type == 3) | (cmap.type == 5)) {
        row.norm <- F
    }
    if (row.norm == TRUE) {
        row.mean <- apply(V, MARGIN = 1, FUN = mean)
        row.sd <- apply(V, MARGIN = 1, FUN = sd)
        row.n <- length(V[, 1])
        for (i in 1:n.rows) {
            if (row.sd[i] == 0) {
                V1[i, ] <- 0
            }
            else {
                V1[i, ] <- (V[i, ] - row.mean[i])/(0.333 * row.sd[i])
            }
            V1[i, ] <- ifelse(V1[i, ] < -6, -6, V1[i, ])
            V1[i, ] <- ifelse(V1[i, ] > 6, 6, V1[i, ])
        }
    }
    else {
        V1 <- V
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
    else if (cmap.type == 6) {
        mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", 
            "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A")
    }
    else if (cmap.type == 7) {
        mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", 
            "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704")
    }
    else if (cmap.type == 8) {
        mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", 
            "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
    }
    else if (cmap.type == 9) {
        mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", 
            "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
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
    ncolors <- length(mycol)
    mycol <- c(mycol, phen.cmap[1:length(col.classes)])
    if (cmap.type == 5) {
        if (max.v == "NA") {
            max.v <- max(max(V1), -min(V1))
        }
        V2 <- ceiling(ncolors * (V1 - (-max.v))/(1.001 * (max.v - 
            (-max.v))))
    }
    else {
        V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001 * (max(V1) - 
            min(V1))))
    }
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
    heatm[n.rows + 1, ] <- ncolors + col.labels
    height <- ifelse(n.rows >= 9, 9, n.rows * 0.44 + 5)
    x11(width = 14, height = height)
    nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow = T), widths = c(10, 
        2), heights = c(5, 1), respect = FALSE)
    par(mar = c(3, 8, 3, 13))
    image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, ncolors + 
        max(col.labels)), col = mycol, axes = FALSE, main = main, 
        sub = sub, xlab = xlab, ylab = ylab)
    if (col.labels2[1] != "NA") {
        groups <- split(col.labels2, col.labels2)
        len.vec <- lapply(groups, length)
        plot.div <- c(0.51, cumsum(len.vec) + 0.5)
        for (i in plot.div) {
            lines(c(i, i), c(0, n.cols), lwd = 2, cex = 0.9, 
                col = "black")
        }
        lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, 
            cex = 0.9, col = "black")
        lines(c(0.51, n.cols + 0.49), c(n.rows + 1.48, n.rows + 
            1.48), lwd = 2, cex = 0.9, col = "black")
    }
    if (row.names[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 25/(n.rows + 15)
        size.col.char <- char.rescale * 35/(n.cols + 15)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 30)
            row.names[i] <- paste(row.names[i], " ", sep = "")
        }
        row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
        axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (row.names2[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 18/(n.rows + 15)
        for (i in 1:n.rows) {
            row.names2[i] <- substr(row.names2[i], 1, 50)
            row.names2[i] <- paste(" ", row.names2[i], sep = "")
        }
        row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "     ")
        axis(4, at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (col.names[1] != "NA") {
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, 
            las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    leg.txt <- col.classes
    p.vec <- rep(22, 22)
    c.vec <- phen.cmap
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = char.rescale * 1.2, pt.cex = char.rescale * 2)
    print(c("range V=", range(V)))
    print(c("range V1=", range(V1)))
    print(c("range V2=", range(V2)))
    par(mar = c(2, 8, 2, 8))
    num.v <- 20
    range.v <- range(V2)
    incr <- (range.v[1] - range.v[2])/(num.v - 1)
    heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), 
        nrow = num.v, ncol = 1)
    image(1:num.v, 1:1, heatm.v, zlim = c(0, ncolors + max(col.labels)), 
        col = mycol, axes = FALSE, main = " ", sub = " ", xlab = ylab, 
        ylab = xlab)
    range.v <- range(V1)
    incr <- (range.v[1] - range.v[2])/(num.v - 1)
    heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], 
        incr)), digits = 2), nrow = num.v, ncol = 1)
    print(c("heatm.v2=", heatm.v2))
    axis(3, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, 
        las = 1, cex.axis = 0.5 * char.rescale, font.axis = 1)
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
