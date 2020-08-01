\name{MSIG.HeatMapPlot.7}
\alias{MSIG.HeatMapPlot.7}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.HeatMapPlot.7
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.HeatMapPlot.7(V, row.names = "NA", row.names2 = "NA", col.labels = "NA", col.labels2 = "NA", col.classes = "NA", phen.cmap = "NA", col.names = "NA", phen.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 0.85, cmap.type = 1, max.v = "NA", legend = T)
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
  \item{phen.names}{
%%     ~~Describe \code{phen.names} here~~
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
  \item{legend}{
%%     ~~Describe \code{legend} here~~
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
    col.names = "NA", phen.names = "NA", main = " ", sub = " ", 
    xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 0.85, 
    cmap.type = 1, max.v = "NA", legend = T) 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    V1 <- matrix(0, nrow = n.rows, ncol = n.cols)
    if (cmap.type == 5) {
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
            V1[i, ] <- ifelse(V1[i, ] < -4, -4, V1[i, ])
            V1[i, ] <- ifelse(V1[i, ] > 4, 4, V1[i, ])
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
        violet.palette <- colorRampPalette(c("#400030", "white"), 
            space = "rgb")
        mycol <- rev(violet.palette(20))
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
    if (col.labels[1] == "NA") {
        heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
        heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
        tot.cols <- ncolors
        if (legend == T) {
            nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow = T), 
                widths = c(5, 1), heights = c(10, 1), respect = FALSE)
        }
        else {
            nf <- layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(8, 
                1), respect = FALSE)
        }
        par(mar = c(3, 16, 3, 16))
        mycol <- c(mycol, phen.cmap[1:length(col.classes)])
        image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), 
            col = mycol, axes = FALSE, main = main, sub = sub, 
            xlab = xlab, ylab = ylab)
        n.rows.phen <- 0
    }
    else {
        tot.cols <- ncolors
        if (is.vector(col.labels)) {
            heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
            heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
            n.rows.phen <- 1
            heatm[n.rows + 1, ] <- tot.cols + col.labels
            cols.row <- length(unique(col.labels))
            tot.cols <- tot.cols + cols.row
            phen.cmap <- phen.cmap[1:cols.row]
        }
        else {
            n.rows.phen <- length(col.labels[, 1])
            cols.row <- vector(length = n.rows.phen, mode = "numeric")
            heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
            heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
            for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
                heatm[k, ] <- tot.cols + col.labels[n.rows + 
                  n.rows.phen - k + 1, ]
                cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + 
                  n.rows.phen - k + 1, ]))
                tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - 
                  k + 1]
            }
            phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, 
                length)))]
        }
        if (legend == T) {
            nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow = T), 
                heights = c(8, 4, 1), respect = FALSE)
        }
        else {
            nf <- layout(matrix(c(1, 2), 2, 1, byrow = T), heights = c(5, 
                1), respect = FALSE)
        }
        par(mar = c(3, 16, 3, 16))
        mycol <- c(mycol, phen.cmap)
        image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, 
            tot.cols), col = mycol, axes = FALSE, main = main, 
            sub = sub, xlab = xlab, ylab = ylab)
    }
    if (col.labels2[1] != "NA") {
        groups <- split(col.labels2, col.labels2)
        len.vec <- lapply(groups, length)
        plot.div <- c(0.51, cumsum(len.vec) + 0.5)
        for (i in plot.div) {
            lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), 
                lwd = 2, cex = 0.9, col = "black")
        }
        lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, 
            cex = 0.9, col = "black")
        lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 
            0.48, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, 
            col = "black")
        lines(c(0.51, n.cols + 0.49), c(n.rows + 0.5, n.rows + 
            0.5), lwd = 2, cex = 0.9, col = "black")
    }
    if (row.names[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 25/(n.rows + 20)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 40)
            row.names[i] <- paste(row.names[i], " ", sep = "")
        }
        if (phen.names[1] == "NA") {
            head.names <- paste("Class", seq(n.rows.phen, 1, 
                -1))
        }
        else {
            head.names <- as.character(rev(phen.names))
        }
        row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
        axis(2, at = 1:(n.rows + n.rows.phen), labels = row.names, 
            adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (row.names2[1] != "NA") {
        numC <- nchar(row.names2)
        size.row.char <- char.rescale * 25/(n.rows + 20)
        for (i in 1:n.rows) {
            row.names2[i] <- substr(row.names2[i], 1, 40)
            row.names2[i] <- paste(" ", row.names2[i], sep = "")
        }
        row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "   ")
        axis(4, at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (col.names[1] != "NA") {
        size.col.char <- char.rescale * 20/(n.cols + 25)
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, 
            las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    if (legend == T) {
        leg.txt <- NULL
        p.vec <- NULL
        c.vec <- NULL
        c2.vec <- NULL
        ind <- 1
        par(mar = c(0, 0, 0, 0))
        plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(-0.05, 
            1.05), axes = F, type = "n", xlab = "", ylab = "")
        for (i in 1:n.rows.phen) {
            if (is.vector(col.labels)) {
                phen.v <- as.character(col.classes)
            }
            else {
                phen.v <- as.character(col.classes[[i]])
            }
            p.name <- paste(as.character(rev(head.names)[i]), 
                ":   ", sep = "")
            leg.txt <- c(p.name, phen.v)
            p.vec <- rep(22, cols.row[i] + 1)
            c.vec <- c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 
                1)])
            c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
            ind <- ind + cols.row[i]
            offset <- 0.07
            legend(x = 0, y = 1 - offset * i, horiz = T, x.intersp = 0.5, 
                legend = leg.txt, bty = "n", xjust = 0, yjust = 1, 
                pch = p.vec, pt.bg = c.vec, col = c2.vec, cex = 1.2, 
                pt.cex = 1.75)
        }
    }
    par(mar = c(2, 12, 2, 12))
    num.v <- 20
    range.v <- range(V2)
    incr <- (range.v[1] - range.v[2])/(num.v - 1)
    heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), 
        nrow = num.v, ncol = 1)
    image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col = mycol, 
        axes = FALSE, main = " ", sub = " ", xlab = ylab, ylab = xlab)
    range.v <- range(V1)
    incr <- (range.v[1] - range.v[2])/(num.v - 1)
    heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], 
        incr)), digits = 2), nrow = num.v, ncol = 1)
    axis(3, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, 
        las = 1, cex.axis = 0.5 * char.rescale, font.axis = 1)
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
