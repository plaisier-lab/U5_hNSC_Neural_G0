\name{MSIG.Extract.Factors}
\alias{MSIG.Extract.Factors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Extract.Factors
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Extract.Factors(input.ds, input.cls = "", output.W.file, output.H.file, k.proj = 2, alg = "NMF.div", niter = 1000, seed = 1234, theta = 0, sort.factors = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.W.file}{
%%     ~~Describe \code{output.W.file} here~~
}
  \item{output.H.file}{
%%     ~~Describe \code{output.H.file} here~~
}
  \item{k.proj}{
%%     ~~Describe \code{k.proj} here~~
}
  \item{alg}{
%%     ~~Describe \code{alg} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{sort.factors}{
%%     ~~Describe \code{sort.factors} here~~
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
function (input.ds, input.cls = "", output.W.file, output.H.file, 
    k.proj = 2, alg = "NMF.div", niter = 1000, seed = 1234, theta = 0, 
    sort.factors = F) 
{
    print(c("Running MSIG.Extract.Factors... on", input.ds))
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dim(m)
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    if (input.cls != "") {
        CLS <- ReadClsFile(file = input.cls)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
    }
    else {
        class.labels <- rep(1, Ns)
        class.phen <- "Samples"
    }
    if (alg == "PCA") {
        svd.proj <- svd(m, nv = Ns, nu = Ns)
        H.full <- diag(x = svd.proj$d, nrow = Ns, ncol = Ns) \%*\% 
            t(svd.proj$v)
        H <- H.full[1:k.proj, ]
        W <- svd.proj$u[, 1:k.proj]
    }
    else if (alg == "NMF.div") {
        NMF.out <- NMF.div(V = m, k = k.proj, maxniter = niter, 
            seed = seed, stopconv = 40, stopfreq = 10)
        H <- NMF.out$H
        W <- NMF.out$W
    }
    else if (alg == "NMF") {
        NMF.out <- NMF(V = m, k = k.proj, maxniter = niter, seed = seed, 
            stopconv = 40, stopfreq = 10)
        H <- NMF.out$H
        W <- NMF.out$W
    }
    else if (alg == "SNMF") {
        NMF.out <- SNMF(V = m, k = k.proj, maxniter = niter, 
            seed = seed, lambda = 1, stopconv = 40, stopfreq = 10)
        H <- NMF.out$H
        W <- NMF.out$W
    }
    else if (alg == "NSNMF.div") {
        NMF.out <- NSNMF.div(V = m, k = k.proj, theta = theta, 
            maxniter = niter, seed = seed, stopconv = 40, stopfreq = 10)
        H <- NMF.out$H
        W <- NMF.out$W
    }
    else {
        stop(c("unknown algorithm:", alg))
    }
    if (sort.factors == T) {
        dist.matrix <- dist(t(W))
        HC <- hclust(dist.matrix, method = "complete")
        W <- W[, HC$order]
        H <- H[HC$order, ]
    }
    factor.names <- paste("F", seq(1, k.proj), sep = "")
    factor.descs <- paste("NMF Extracted Factor Number ", seq(1, 
        k.proj), sep = "")
    V <- data.frame(W)
    names(V) <- factor.names
    row.names(V) <- gs.names
    write.gct(gct.data.frame = V, descs = gs.descs, filename = output.W.file)
    V <- data.frame(H)
    names(V) <- sample.names
    row.names(V) <- factor.names
    write.gct(gct.data.frame = V, descs = factor.descs, filename = output.H.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
