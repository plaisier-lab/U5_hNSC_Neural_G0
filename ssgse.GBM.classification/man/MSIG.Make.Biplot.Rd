\name{MSIG.Make.Biplot}
\alias{MSIG.Make.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Make.Biplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Make.Biplot(m.train, labels.train, phen.train, m.test = NULL, labels.test = NULL, phen.test = NULL, row.names = NULL, cols = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m.train}{
%%     ~~Describe \code{m.train} here~~
}
  \item{labels.train}{
%%     ~~Describe \code{labels.train} here~~
}
  \item{phen.train}{
%%     ~~Describe \code{phen.train} here~~
}
  \item{m.test}{
%%     ~~Describe \code{m.test} here~~
}
  \item{labels.test}{
%%     ~~Describe \code{labels.test} here~~
}
  \item{phen.test}{
%%     ~~Describe \code{phen.test} here~~
}
  \item{row.names}{
%%     ~~Describe \code{row.names} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
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
function (m.train, labels.train, phen.train, m.test = NULL, labels.test = NULL, 
    phen.test = NULL, row.names = NULL, cols = NULL) 
{
    if (is.null(cols)) 
        cols <- sample(colors(), 20)
    pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)
    S1 <- pca$x[, 1]
    S2 <- pca$x[, 2]
    X1 <- pca$rotation[, 1]
    X2 <- pca$rotation[, 2]
    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1/max.X
    X2 <- max.S * X2/max.X
    max.A <- max(max.S, max.X)
    color <- cols[class.labels]
    num.samples <- length(S1)
    num.rows <- length(m.train[, 1])
    class.labels <- labels.train
    if (min(class.labels) == 0) 
        class.labels <- class.labels + 1
    class.phen <- phen.train
    if (!is.null(m.test)) {
        test.scores <- predict(pca, t(m.test))
        S1 <- c(pca$x[, 1], test.scores[, 1])
        S2 <- c(pca$x[, 2], test.scores[, 2])
        max.S <- max(sqrt(S1 * S1 + S2 * S2))
        max.X <- max(sqrt(X1 * X1 + X2 * X2))
        X1 <- max.S * X1/max.X
        X2 <- max.S * X2/max.X
        num.samples <- length(S1)
        if (min(labels.test) == 0) 
            labels.test <- labels.test + 1
        class.labels <- c(labels.train, labels.test + max(labels.train))
        class.phen <- c(class.phen, phen.test)
    }
    x11(height = 9, width = 14)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 
        1), heights = 1, respect = FALSE)
    plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), 
        type = "n", main = "Biplot", sub = input.ds)
    for (j in 1:num.samples) {
        if (!is.null(m.test)) {
            if (j <= length(labels.train)) 
                pch <- 22
            else pch <- 21
        }
        else {
            pch <- 22
        }
        points(S1[j], S2[j], pch = pch, type = "p", cex = 2, 
            bg = cols[class.labels[j]], col = "black")
    }
    for (j in 1:num.rows) {
        x.coor <- X1[j] * 0.925
        y.coor <- X2[j] * 0.925
        arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, 
            angle = 20, col = "black")
        leg.txt <- ifelse(!is.null(row.names), row.names[j], 
            paste("F", j, sep = ""))
        text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, 
            offset = 0.5, vfont = NULL, cex = 1, col = "grey")
    }
    ang <- vector(length = k.proj, mode = "numeric")
    for (j in 1:k.proj) {
        ang[j] <- ifelse(atan2(X2[j], X1[j]) > 0, atan2(X2[j], 
            X1[j]), 2 * pi + atan2(X2[j], X1[j]))
    }
    ang.index <- order(ang, decreasing = F)
    ang2 <- ang[ang.index]
    for (j in 1:num.rows) {
        if (j == num.rows) {
            angle.in.between <- (ang2[1] - ang2[j] - 2 * pi)/2 + 
                ang2[j] - 2 * pi
        }
        else {
            angle.in.between <- (ang2[j + 1] - ang2[j])/2 + ang2[j]
        }
        x <- max.S * cos(angle.in.between)
        y <- max.S * sin(angle.in.between)
        arrows(0, 0, x, y, lwd = 4, length = 0, lty = 3, col = "grey50")
    }
    if (!is.null(m.test)) {
        leg.txt <- c("Train:", phen.train, "Test:", phen.test)
        p.vec <- c(rep(22, length(phen.train) + 1), rep(21, length(phen.test) + 
            1))
        c.vec <- c("white", cols[1:length(phen.train)], "white", 
            cols[seq(length(phen.train) + 1, length(class.phen))])
        b.vec <- c("white", rep("black", length(phen.train)), 
            "white", rep("black", length(phen.test)))
    }
    else {
        leg.txt <- class.phen
        n.phen <- length(class.phen)
        p.vec <- rep(21, n.phen)
        c.vec <- cols[1:n.phen]
        b.vec <- rep("black", n.phen)
    }
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = b.vec, cex = 1.5, 
        pt.cex = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
