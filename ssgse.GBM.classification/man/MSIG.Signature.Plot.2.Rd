\name{MSIG.Signature.Plot.2}
\alias{MSIG.Signature.Plot.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Signature.Plot.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Signature.Plot.2(V, row.names = "NA", row.names2 = "NA", col.labels = "NA", col.classes = "NA", phen.cmap = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", char.rescale = 1, max.v = "NA", seed = 1729)
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
  \item{char.rescale}{
%%     ~~Describe \code{char.rescale} here~~
}
  \item{max.v}{
%%     ~~Describe \code{max.v} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
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
    col.classes = "NA", phen.cmap = "NA", col.names = "NA", main = " ", 
    sub = " ", xlab = " ", ylab = " ", char.rescale = 1, max.v = "NA", 
    seed = 1729) 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    V[V > 4] <- 4
    V[V < -4] <- -4
    mycol <- vector(length = 512, mode = "numeric")
    for (k in 1:256) {
        mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
        mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
    ncolors <- length(mycol)
    if (is.na(max(V))) {
        mycol <- c(mycol, phen.cmap[1:(length(col.classes) + 
            1)])
    }
    else {
        mycol <- c(mycol, phen.cmap[1:length(col.classes)])
    }
    V2 <- ceiling(ncolors * (V - -4)/(8 + 1e-04))
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
    for (i in 1:length(heatm[, 1])) {
        for (j in 1:length(heatm[1, ])) {
            heatm[i, j] <- ifelse(is.na(heatm[i, j]), ncolors + 
                max(col.labels) + 1, heatm[i, j])
        }
    }
    heatm[n.rows + 1, ] <- ncolors + col.labels
    height <- ifelse(n.rows >= 25, 25, n.rows * 0.8 + 5)
    x11(width = 31, height = 19)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(7, 
        1), respect = FALSE)
    par(mar = c(5, 7, 4, 7))
    if (is.na(max(V))) {
        image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, 
            ncolors + max(col.labels) + 1), col = mycol, axes = FALSE, 
            main = main, sub = sub, xlab = xlab, ylab = ylab)
    }
    else {
        image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, 
            ncolors + max(col.labels)), col = mycol, axes = FALSE, 
            main = main, sub = sub, xlab = xlab, ylab = ylab)
    }
    if (row.names[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 30/(n.rows + 15)
        size.col.char <- char.rescale * 20/(n.cols + 15)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 30)
        }
        row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
        axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (row.names2[1] != "NA") {
        numC <- nchar(row.names2)
        size.row.char <- char.rescale * 35/(n.rows + 15)
        size.col.char <- char.rescale * 20/(n.cols + 15)
        for (i in 1:n.rows) {
            row.names2[i] <- substr(row.names2[i], 1, 30)
        }
        row.names2 <- c(row.names2[seq(n.rows, 1, -1)], "Class")
        axis(4, at = 1:(n.rows + 1), labels = row.names2, adj = 0.5, 
            tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    if (col.names[1] != "NA") {
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, 
            las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    if (col.classes[1] != "NA") {
        C <- split(col.labels, col.labels)
        class1.size <- length(C[[1]])
        class2.size <- length(C[[2]])
        axis(3, at = c(class1.size/2, class1.size + class2.size/2), 
            labels = col.classes, tick = FALSE, las = 1, cex.axis = 1.25, 
            font.axis = 2, line = -1)
    }
    par(mar = c(10, 2, 10, 2))
    num.v <- 20
    range.v <- range(V2, na.rm = T)
    incr <- (range.v[1] - range.v[2])/(num.v - 1)
    heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), 
        nrow = num.v, ncol = 1)
    if (is.na(max(V))) {
        image(1:1, 1:num.v, t(heatm.v), zlim = c(0, ncolors + 
            max(col.labels) + 1), col = mycol, axes = FALSE, 
            sub = "Color \n Legend ", main = " ", xlab = xlab, 
            ylab = ylab)
    }
    else {
        image(1:1, 1:num.v, t(heatm.v), zlim = c(0, ncolors + 
            max(col.labels)), col = mycol, axes = FALSE, sub = "Color \n Legend ", 
            main = " ", xlab = xlab, ylab = ylab)
    }
    range.v <- range(V, na.rm = T)
    incr <- (range.v[1] - range.v[2])/(num.v - 1)
    heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], 
        incr)), digits = 3), nrow = num.v, ncol = 1)
    axis(2, at = 1:num.v, labels = heatm.v2, adj = 0.5, tick = FALSE, 
        las = 1, cex.axis = char.rescale * 0.6, font.axis = 1.25, 
        line = -0.8)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
