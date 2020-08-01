\name{MSIG.HeatMapPlot.5}
\alias{MSIG.HeatMapPlot.5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.HeatMapPlot.5
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.HeatMapPlot.5(V, row.names = "NA", col.labels = "NA", col.classes = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 1, cmap.type = 1, max.v = "NA", rotated.col.labels = F, create.legend = T, create.window = T)
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
  \item{char.rescale}{
%%     ~~Describe \code{char.rescale} here~~
}
  \item{cmap.type}{
%%     ~~Describe \code{cmap.type} here~~
}
  \item{max.v}{
%%     ~~Describe \code{max.v} here~~
}
  \item{rotated.col.labels}{
%%     ~~Describe \code{rotated.col.labels} here~~
}
  \item{create.legend}{
%%     ~~Describe \code{create.legend} here~~
}
  \item{create.window}{
%%     ~~Describe \code{create.window} here~~
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
    row.norm = TRUE, char.rescale = 1, cmap.type = 1, max.v = "NA", 
    rotated.col.labels = F, create.legend = T, create.window = T) 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    n.phen <- length(col.classes)
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
        mycol <- vector(length = 256, mode = "numeric")
        for (k in 1:256) {
            red <- (k - 1) * 0.8 + 50
            green <- (k - 1) * 0.68 + 80
            blue <- (k - 1) * 0.6 + 100
            mycol[k] <- rgb(red, green, blue, maxColorValue = 255)
        }
        mycol <- rev(mycol)
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
        mycol <- vector(length = 256, mode = "numeric")
        max.pos.V1 <- max(V1)
        if (min(V1) < 0) {
            min.neg.V1 <- min(V1)
        }
        else {
            min.neg.V1 <- 0
        }
        neg.k <- ceiling(256 * (abs(min.neg.V1)/(max.pos.V1 - 
            min.neg.V1)))
        for (k in 1:neg.k) {
            max.red <- 255 - (255 - 50) * abs(min.neg.V1)
            min.red <- 255
            red <- max.red + (min.red - max.red) * (k - 1)/(neg.k - 
                1)
            max.green <- (255 - (255 - 50) * abs(min.neg.V1))
            min.green <- 255
            green <- max.green + (min.green - max.green) * (k - 
                1)/(neg.k - 1)
            max.blue <- (255 - (255 - 200) * abs(min.neg.V1))
            min.blue <- 255
            blue <- max.blue + (min.blue - max.blue) * (k - 1)/(neg.k - 
                1)
            mycol[k] <- rgb(red, green, blue, maxColorValue = 255)
        }
        for (k in (neg.k + 1):256) {
            max.red <- 205
            min.red <- 255
            red <- min.red - (min.red - max.red) * (k - (neg.k + 
                1))/(256 - (neg.k + 1))
            max.green <- 50
            min.green <- 255
            green <- min.green - (min.green - max.green) * (k - 
                (neg.k + 1))/(256 - (neg.k + 1))
            max.blue <- 50
            min.blue <- 255
            blue <- min.blue - (min.blue - max.blue) * (k - (neg.k + 
                1))/(256 - (neg.k + 1))
            mycol[k] <- rgb(red, green, blue, maxColorValue = 255)
        }
    }
    ncolors <- length(mycol)
    if (cmap.type == 5) {
        if (max.v == "NA") {
            max.v <- max(max(V1), -min(V1))
        }
        V2 <- ceiling(ncolors * (V1 - (-max.v))/(1.001 * (max.v - 
            (-max.v))))
    }
    else if (cmap.type == 3) {
        V2 <- ceiling(ncolors * (V1/1.001))
    }
    else {
        V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001 * (max(V1) - 
            min(V1))))
    }
    heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
    heatm[1:n.rows, ] <- V2[, seq(n.cols, 1, -1)]
    col.labels <- col.labels[seq(n.cols, 1, -1)]
    if (col.names[1] != "NA") {
        col.names <- col.names[seq(n.cols, 1, -1)]
    }
    height <- ifelse(n.rows >= 25, 25, n.rows * 0.8 + 2)
    if (create.window == T) {
        x11(width = 31, height = 19)
    }
    if ((create.window == T && create.legend == T)) {
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(7, 
            1), respect = FALSE)
    }
    n.rows2 <- ifelse(n.rows < 4, 4, n.rows)
    a <- -12/16
    b <- 27
    margin <- a * n.rows2 + b
    margin <- ifelse(margin < 2, 2, ifelse(margin > 24, 24, margin))
    if (rotated.col.labels == F) {
        par(mar = c(4, margin, 4, 10))
    }
    else {
        par(mar = c(4, margin, 10, 10))
    }
    image(1:n.rows, 1:n.cols, heatm, zlim = c(0, ncolors), col = mycol, 
        axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
    if (row.names[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 15/(n.rows + 15)
        size.col.char <- char.rescale * 30/(n.cols + 15)
        size.lab.char <- char.rescale * 30/(n.phen + 15)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 20)
        }
        if (rotated.col.labels == F) {
            axis(3, at = 1:n.rows, labels = row.names, adj = 1, 
                tick = FALSE, las = 1, cex.axis = size.row.char, 
                font.axis = 1, line = -0.75)
        }
        else {
            axis(3, at = 1:n.rows, labels = row.names, adj = 1, 
                tick = FALSE, las = 3, cex.axis = size.row.char, 
                font.axis = 1, line = -0.75)
        }
    }
    if (col.names[1] != "NA") {
        axis(2, at = 1:n.cols, labels = col.names, tick = FALSE, 
            las = 1, cex.axis = size.col.char, font.axis = 2.5, 
            line = -1)
    }
    for (i in 1:(n.rows)) {
        lines(x = c(i + 0.5, i + 0.5), y = c(0.5, n.cols + 0.5), 
            type = "l", lwd = 1, col = "black")
    }
    boundaries <- cumsum(sapply(split(rep(1, n.cols), col.labels), 
        sum))
    boundaries <- n.cols - boundaries
    loc.phen <- vector(length = n.phen, mode = "numeric")
    for (i in 1:(n.phen)) {
        lines(x = c(1 - 0.5, n.rows + 0.5), y = c(boundaries[i] + 
            0.5, boundaries[i] + 0.5), type = "l", lwd = 1, col = "black")
        if (i > 1) {
            loc.phen[i] <- mean(c(boundaries[i - 1], boundaries[i]))
        }
        else {
            loc.phen[i] <- mean(c(n.cols, boundaries[i]))
        }
    }
    axis(4, at = loc.phen, labels = col.classes, tick = FALSE, 
        las = 1, cex.axis = size.lab.char, font.axis = 1, line = -0.75)
    lines(x = c(0.5, n.rows + 0.5), y = c(0.5, 0.5), type = "l", 
        lwd = 1, col = "black")
    lines(x = c(0.5, n.rows + 0.5), y = c(n.cols + 0.5, n.cols + 
        0.5), type = "l", lwd = 1, col = "black")
    lines(x = c(0.5, 0.5), y = c(0.5, n.cols + 0.5), type = "l", 
        lwd = 1, col = "black")
    lines(x = c(n.rows + 0.5, n.rows + 0.5), y = c(0.5, n.cols + 
        0.5), type = "l", lwd = 1, col = "black")
    if (create.legend == T) {
        par(mar = c(10, 2, 10, 2))
        num.v <- 20
        if (cmap.type == 3) {
            range.v <- c(0, ncolors)
        }
        else {
            range.v <- range(V2)
        }
        incr <- (range.v[1] - range.v[2])/(num.v - 1)
        heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), 
            nrow = num.v, ncol = 1)
        image(1:1, 1:num.v, t(heatm.v), zlim = c(0, ncolors + 
            max(col.labels)), col = mycol, axes = FALSE, sub = "Color \n Legend ", 
            main = " ", xlab = xlab, ylab = ylab)
        if (cmap.type == 3) {
            range.v <- c(0, 1)
        }
        else {
            range.v <- range(V1)
        }
        incr <- (range.v[1] - range.v[2])/(num.v - 1)
        heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], 
            incr)), digits = 3), nrow = num.v, ncol = 1)
        axis(2, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, 
            las = 1, cex.axis = char.rescale * 0.6, font.axis = 1.25, 
            line = -0.8)
        lines(x = c(0.5, 1.5), y = c(0.5, 0.5), type = "l", lwd = 1, 
            col = "black")
        lines(x = c(0.5, 1.5), y = c(num.v + 0.5, num.v + 0.5), 
            type = "l", lwd = 1, col = "black")
        lines(x = c(0.6, 0.6), y = c(0.5, num.v + 0.5), type = "l", 
            lwd = 1, col = "black")
        lines(x = c(1.4, 1.4), y = c(0.5, num.v + 0.5), type = "l", 
            lwd = 1, col = "black")
    }
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
