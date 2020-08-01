\name{MSIG.Preprocess.Dataset}
\alias{MSIG.Preprocess.Dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Preprocess.Dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Preprocess.Dataset(input.ds, output.ds, thres = NULL, ceil = NULL, shift = NULL, fold = NULL, delta = NULL, normalization = NULL, cntrl.genes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
}
  \item{ceil}{
%%     ~~Describe \code{ceil} here~~
}
  \item{shift}{
%%     ~~Describe \code{shift} here~~
}
  \item{fold}{
%%     ~~Describe \code{fold} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{normalization}{
%%     ~~Describe \code{normalization} here~~
}
  \item{cntrl.genes}{
%%     ~~Describe \code{cntrl.genes} here~~
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
function (input.ds, output.ds, thres = NULL, ceil = NULL, shift = NULL, 
    fold = NULL, delta = NULL, normalization = NULL, cntrl.genes = NULL) 
{
    print(c("Running MSIG.Preprocess.Dataset... on:", input.ds))
    print(c("output file:", output.ds))
    print(c("normalization =", normalization))
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    if (!is.null(thres)) {
        m[m < thres] <- thres
    }
    if (!is.null(ceil)) {
        m[m > ceil] <- ceil
    }
    if (!is.null(shift)) {
        m <- m + shift
    }
    if (!is.null(cntrl.genes)) {
        gene.names2 <- intersect(cntrl.genes, gs.names)
        locs <- match(gene.names2, gs.names, nomatch = 0)
        msig.cntrl <- m[locs, ]
        msig.cntrl.genes <- gs.names[locs]
        msig.cntrl.descs <- gs.descs[locs]
        m <- m[-locs, ]
        gs.names <- gs.names[-locs]
        gs.descs <- gs.descs[-locs]
    }
    if ((!is.null(fold)) && (!is.null(delta))) {
        temp <- MSIG.VarFilter(V = m, fold = fold, delta = delta, 
            gene.names = gs.names, gene.descs = gs.descs)
        m <- temp$V
        gs.names <- temp$new.gene.names
        gs.descs <- temp$new.gene.descs
        dim(m)
    }
    if (!is.null(cntrl.genes)) {
        m <- rbind(m, msig.cntrl)
        gs.names <- c(gs.names, msig.cntrl.genes)
        gs.descs <- c(gs.descs, msig.cntrl.descs)
    }
    if (!is.null(normalization)) {
        if (normalization == 1) {
            m <- MSIG.NormalizeCols.Rank(m)
        }
        else if (normalization == 2) {
            m <- MSIG.NormalizeCols.Rank(m)/length(m[, 1])
        }
        else if (normalization == 3) {
            m <- GSEA.NormalizeCols(m) + 3
            m <- GSEA.Threshold(m, 0.001, 1e+05)
        }
        else if (normalization == 4) {
            m <- MSIG.NormalizeCols.Rank(m)/length(m[, 1])
        }
        else if (normalization == 5) {
            m <- MSIG.NormalizeCols.Rescale(m)
        }
        else if (normalization == 6) {
            cols <- length(m[1, ])
            for (j in 1:cols) {
                m[, j] <- rank(m[, j], ties.method = "average") - 
                  1
            }
            m <- 10000 * m/(length(m[, 1]) - 1)
        }
        else if (normalization == 7) {
            m <- ((100 * MSIG.NormalizeCols.Rank(m))\%/\%length(m[, 
                1]) + 1)
        }
        else if (normalization == 8) {
            row.mean <- apply(m, MARGIN = 1, FUN = mean)
            for (i in 1:length(m[, 1])) {
                m[i, ] <- m[i, ]/row.mean[i]
            }
        }
    }
    V <- data.frame(m)
    names(V) <- sample.names
    row.names(V) <- gs.names
    write.gct(gct.data.frame = V, descs = gs.descs, filename = output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
