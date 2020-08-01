\name{MSIG.Signature.Plot}
\alias{MSIG.Signature.Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Signature.Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Signature.Plot(V, row.names = "NA", col.labels = "NA", col.classes = "NA", phen.cmap = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 1, max.v = "NA", seed = 1729)
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
function (V, row.names = "NA", col.labels = "NA", col.classes = "NA", 
    phen.cmap = "NA", col.names = "NA", main = " ", sub = " ", 
    xlab = " ", ylab = " ", row.norm = TRUE, char.rescale = 1, 
    max.v = "NA", seed = 1729) 
{
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    V1 <- matrix(0, nrow = n.rows, ncol = n.cols)
    if (row.norm == TRUE) {
        row.mean <- apply(V, MARGIN = 1, FUN = mean)
        row.sd <- apply(V, MARGIN = 1, FUN = sd)
        row.n <- length(V[, 1])
        for (i in 1:n.rows) {
            if (row.sd[i] == 0) {
                V1[i, ] <- 0
            }
            else {
                V1[i, ] <- (V[i, ] - row.mean[i])/row.sd[i]
            }
            V1[i, ] <- ifelse(V1[i, ] < -4, -4, V1[i, ])
            V1[i, ] <- ifelse(V1[i, ] > 4, 4, V1[i, ])
        }
    }
    else {
        V1 <- V
    }
    mycol <- vector(length = 512, mode = "numeric")
    for (k in 1:256) {
        mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue = 255)
    }
    for (k in 257:512) {
        mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue = 255)
    }
    mycol <- rev(mycol)
    ncolors <- length(mycol)
    mycol <- c(mycol, phen.cmap[1:length(col.classes)])
    V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001 * (max(V1) - 
        min(V1))))
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V2[seq(n.rows, 1, -1), ]
    heatm[n.rows + 1, ] <- ncolors + col.labels
    height <- ifelse(n.rows >= 25, 25, n.rows * 0.8 + 5)
    x11(width = 19, height = 11)
    nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), widths = c(4, 
        4, 1), heights = 1, respect = FALSE)
    par(mar = c(8, 9, 4, 4))
    image(1:n.cols, 1:(n.rows + 1), t(heatm), zlim = c(0, ncolors + 
        max(col.labels)), col = mycol, axes = FALSE, main = main, 
        sub = sub, xlab = xlab, ylab = ylab)
    if (row.names[1] != "NA") {
        numC <- nchar(row.names)
        size.row.char <- char.rescale * 35/(n.rows + 15)
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
    par(mar = c(8, 2, 4, 2))
    library(RColorBrewer)
    col.list <- c(brewer.pal(n = 9, name = "Oranges"), brewer.pal(n = 9, 
        name = "Reds"), brewer.pal(n = 9, name = "Blues"), brewer.pal(n = 9, 
        name = "Greens"), brewer.pal(n = 9, name = "Greys"), 
        brewer.pal(n = 9, name = "Purples"), brewer.pal(n = 9, 
            name = "Pastel1"), brewer.pal(n = 9, name = "Set1"), 
        brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 8, 
            name = "Accent"), brewer.pal(n = 10, name = "Spectral"), 
        brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"), 
        brewer.pal(n = 8, name = "Set3"), brewer.pal(n = 8, name = "BuGn"), 
        brewer.pal(n = 8, name = "GnBu"), brewer.pal(n = 8, name = "OrRd"), 
        brewer.pal(n = 8, name = "PuBu"), brewer.pal(n = 8, name = "PuBuGn"), 
        brewer.pal(n = 8, name = "RdPu"), brewer.pal(n = 8, name = "YlGn"), 
        brewer.pal(n = 8, name = "YlOrBr"), brewer.pal(n = 8, 
            name = "YlOrRd"))
    set.seed(seed)
    sig.col <- sample(col.list, size = 1000, replace = T)
    h.range <- range(apply(V2, MARGIN = 2, FUN = sum))
    h.size <- h.range[2]/30
    V3 <- rbind(rep(h.size, n.cols), apply(V2, MARGIN = 2, FUN = rev))
    barplot(V3, main = main, font.axis = 1.25, cex.lab = 1.5, 
        cex.axis = 1.25, cex.names = 1.25, width = 1, space = 0, 
        border = 1, col = sig.col)
    barplot(rep(h.size, n.cols), cex.lab = 1.5, cex.axis = 1.25, 
        cex.names = 1.25, width = 1, space = 0, border = 0, col = phen.cmap[col.labels], 
        add = T)
    leg.txt <- col.classes
    p.vec <- rep(21, 21)
    c.vec <- phen.cmap
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = char.rescale * 1.3, pt.cex = char.rescale * 2.5)
    return()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
