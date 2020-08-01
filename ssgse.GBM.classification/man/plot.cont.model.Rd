\name{plot.cont.model}
\alias{plot.cont.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.cont.model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.cont.model(feature, feature.mean, feature.sd, num.points.extrapolation = 500, x, prob.i, beta0, beta1, col.vec, color.map = c("red", "green"), target.feature)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{feature}{
%%     ~~Describe \code{feature} here~~
}
  \item{feature.mean}{
%%     ~~Describe \code{feature.mean} here~~
}
  \item{feature.sd}{
%%     ~~Describe \code{feature.sd} here~~
}
  \item{num.points.extrapolation}{
%%     ~~Describe \code{num.points.extrapolation} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{prob.i}{
%%     ~~Describe \code{prob.i} here~~
}
  \item{beta0}{
%%     ~~Describe \code{beta0} here~~
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
}
  \item{col.vec}{
%%     ~~Describe \code{col.vec} here~~
}
  \item{color.map}{
%%     ~~Describe \code{color.map} here~~
}
  \item{target.feature}{
%%     ~~Describe \code{target.feature} here~~
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
function (feature, feature.mean, feature.sd, num.points.extrapolation = 500, 
    x, prob.i, beta0, beta1, col.vec, color.map = c("red", "green"), 
    target.feature) 
{
    range.x <- max(x) - min(x)
    prob.m <- matrix(0, nrow = num.points.extrapolation, ncol = 3)
    x.m <- vector(length = num.points.extrapolation, mode = "numeric")
    for (j in 1:num.points.extrapolation) {
        x.m[j] <- min(x) + j * (range.x/500)
        p.vec <- exp(beta0 + beta1 * x.m[j])/(1 + exp(beta0 + 
            beta1 * x.m[j]))
        prob.m[j, ] <- quantile(p.vec, probs = c(0.05, 0.5, 0.95))
    }
    istar <- min(x) + which.min(abs(0.5 - prob.m[, 2])) * (range.x/num.points.extrapolation)
    x.index <- order(x, decreasing = F)
    x.order <- x[x.index]
    prob.i.order <- prob.i[x.index, ]
    col.vec.order <- col.vec[x.index]
    col.vec.order <- ifelse(col.vec.order == 1, color.map[1], 
        color.map[2])
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.5, 
        1), heights = 1, respect = FALSE)
    plot(x.order, prob.i.order[, 2], main = feature, sub = paste("Boundary: ", 
        signif(istar, 3)), pch = 20, ylim = c(-0.3, 1.07), col = 0, 
        cex = 2, xlab = "Activation score", ylab = "Probability")
    for (h in 1:length(x.m)) points(c(x.m[h], x.m[h]), c(prob.m[h, 
        1], prob.m[h, 3]), type = "l", col = "gray90", lty = 1, 
        lwd = 2)
    points(x.m, prob.m[, 2], type = "l", lwd = 2, col = 1, lty = 1, 
        cex = 1)
    points(x.m, prob.m[, 1], type = "l", col = 1, lty = 1, cex = 1)
    points(x.m, prob.m[, 3], type = "l", col = 1, lty = 1, cex = 1)
    range.x <- range(x.order)
    points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = "gray", 
        lwd = 2)
    points(range.x, c(0, 0), type = "l", lty = 1, col = 1, lwd = 1)
    points(c(istar, istar), c(-0.3, 1.07), type = "l", lty = 3, 
        col = "gray", lwd = 2)
    points(x.order, prob.i.order[, 2], pch = 21, bg = col.vec.order, 
        col = 1, cex = 1.5)
    red.points <- x.order[col.vec.order == color.map[1]]
    green.points <- x.order[col.vec.order == color.map[2]]
    points(range.x, c(-0.1, -0.1), type = "l", lty = 1, col = 1, 
        lwd = 1)
    points(range.x, c(-0.2, -0.2), type = "l", lty = 1, col = 1, 
        lwd = 1)
    points(red.points, rep(-0.1, length(red.points)), pch = 21, 
        bg = color.map[1], col = 1, cex = 1.5)
    points(green.points, rep(-0.2, length(green.points)), pch = 21, 
        bg = color.map[2], col = 1, cex = 1.5)
    leg.txt <- c(target.class, other.class)
    p.vec <- rep(21, 21)
    c.vec <- color.map
    old.par <- par(no.readonly = TRUE)
    par(mar = c(1, 1, 1, 1))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    text(x = 0.5, y = 0.8, labels = target.feature, cex = 1.3)
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.2, pt.cex = 1.2)
    par(old.par)
    return(istar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
