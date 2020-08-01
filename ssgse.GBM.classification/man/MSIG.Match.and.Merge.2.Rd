\name{MSIG.Match.and.Merge.2}
\alias{MSIG.Match.and.Merge.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Match.and.Merge.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Match.and.Merge.2(input1.ds, input1.cls = "", input2.ds, input2.cls = "", output.ds, output.cls, mode = "intersection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input1.ds}{
%%     ~~Describe \code{input1.ds} here~~
}
  \item{input1.cls}{
%%     ~~Describe \code{input1.cls} here~~
}
  \item{input2.ds}{
%%     ~~Describe \code{input2.ds} here~~
}
  \item{input2.cls}{
%%     ~~Describe \code{input2.cls} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{output.cls}{
%%     ~~Describe \code{output.cls} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
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
function (input1.ds, input1.cls = "", input2.ds, input2.cls = "", 
    output.ds, output.cls, mode = "intersection") 
{
    print(c("Running MSIG.Match.and.Merge... on: ", input1.ds, 
        " ", input2.ds))
    dataset1 <- MSIG.Gct2Frame(filename = input1.ds)
    m1 <- data.matrix(dataset1$ds)
    gs.names1 <- dataset1$row.names
    gs.descs1 <- dataset1$descs
    sample.names1 <- dataset1$names
    N1 <- length(m1[1, ])
    Ng1 <- length(m1[, 1])
    print(c("dataset 1:", N1, " samples"))
    print(c("dataset 1:", Ng1, " genes"))
    dataset2 <- MSIG.Gct2Frame(filename = input2.ds)
    m2 <- data.matrix(dataset2$ds)
    gs.names2 <- dataset2$row.names
    gs.descs2 <- dataset2$descs
    sample.names2 <- dataset2$names
    N2 <- length(m2[1, ])
    Ng2 <- length(m2[, 1])
    print(c("dataset 2:", N2, " samples"))
    print(c("dataset 2:", Ng2, " genes"))
    if ((input1.cls != "") & (input2.cls != "")) {
        CLS1 <- ReadClsFile(file = input1.cls)
        class.labels1 <- CLS1$class.v
        class.phen1 <- CLS1$phen
        CLS2 <- ReadClsFile(file = input2.cls)
        class.labels2 <- CLS2$class.v
        class.phen2 <- CLS2$phen
    }
    sample.names3 <- c(sample.names1, sample.names2)
    N3 <- N1 + N2
    if (mode == "intersection") {
        gs.names3 <- intersect(gs.names1, gs.names2)
        locations1 <- match(gs.names3, gs.names1, nomatch = 0)
        m1 <- m1[locations1, ]
        gs.descs1 <- gs.descs1[locations1]
        locations2 <- match(gs.names3, gs.names2, nomatch = 0)
        m2 <- m2[locations2, ]
        gs.descs2 <- gs.descs2[locations2]
        m3 <- cbind(m1, m2)
    }
    else if (mode == "union") {
        gs.names3 <- union(gs.names1, gs.names2)
        M3 <- length(gs.names3)
        m3 <- matrix(0, nrow = M3, ncol = N3)
        locations1 <- match(gs.names3, gs.names1, nomatch = 0)
        locations2 <- match(gs.names3, gs.names2, nomatch = 0)
        for (i in 1:M3) {
            if (locations1[i] != 0) {
                m3[i, 1:N1] <- m1[locations1[i], ]
            }
            if (locations2[i] != 0) {
                m3[i, (N1 + 1):N3] <- m2[locations2[i], ]
            }
        }
    }
    else {
        stop(c("unknown mode", mode))
    }
    if ((input1.cls != "") & (input2.cls != "")) {
        class.labels3 <- c(class.labels1, class.labels2 + length(class.phen1))
        class.phen3 <- c(class.phen1, class.phen2)
    }
    N3 <- length(m3[1, ])
    Ng3 <- length(m3[, 1])
    print(c("dataset 3:", N3, " samples"))
    print(c("dataset 3:", Ng3, " genes"))
    V <- data.frame(m3)
    names(V) <- sample.names3
    row.names(V) <- gs.names3
    gs.descs1 <- gs.names3
    write.gct(gct.data.frame = V, descs = gs.descs1, filename = output.ds)
    if ((input1.cls != "") & (input2.cls != "")) {
        write.cls(class.v = class.labels3, phen = class.phen3, 
            filename = output.cls)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
