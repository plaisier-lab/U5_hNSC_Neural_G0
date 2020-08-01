\name{MSIG.Score.Plot}
\alias{MSIG.Score.Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Score.Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Score.Plot(z, main = "", phen.cmap, char.rescale = 1, col.classes, col.labels, create.legend = T, create.window = T, xlab = " ", ylab = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{phen.cmap}{
%%     ~~Describe \code{phen.cmap} here~~
}
  \item{char.rescale}{
%%     ~~Describe \code{char.rescale} here~~
}
  \item{col.classes}{
%%     ~~Describe \code{col.classes} here~~
}
  \item{col.labels}{
%%     ~~Describe \code{col.labels} here~~
}
  \item{create.legend}{
%%     ~~Describe \code{create.legend} here~~
}
  \item{create.window}{
%%     ~~Describe \code{create.window} here~~
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
function (z, main = "", phen.cmap, char.rescale = 1, col.classes, 
    col.labels, create.legend = T, create.window = T, xlab = " ", 
    ylab = " ") 
{
    size <- length(z)
    if (create.window == T) {
        x11(width = 19, height = 11)
    }
    if ((create.window == T && create.legend == T)) {
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(6, 
            1), heights = 1, respect = FALSE)
    }
    barplot(z, xlab = xlab, ylab = ylab, main = main, font.axis = 1, 
        cex.lab = 1, cex.axis = 1, cex.names = 1, width = 1, 
        space = 0, col = phen.cmap[col.labels])
    if (create.legend == T) {
        leg.txt <- col.classes
        p.vec <- rep(21, 21)
        c.vec <- phen.cmap
        par(mar = c(0, 0, 0, 0))
        plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), 
            axes = F, type = "n", xlab = "", ylab = "")
        legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
            yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
            cex = char.rescale, pt.cex = char.rescale * 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
