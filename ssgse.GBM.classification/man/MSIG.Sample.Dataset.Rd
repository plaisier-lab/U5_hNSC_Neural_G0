\name{MSIG.Sample.Dataset}
\alias{MSIG.Sample.Dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Sample.Dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Sample.Dataset(input.ds, input.cls = "", column.subset.fraction = "ALL", row.subset.fraction = "ALL", output.ds, output.cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{column.subset.fraction}{
%%     ~~Describe \code{column.subset.fraction} here~~
}
  \item{row.subset.fraction}{
%%     ~~Describe \code{row.subset.fraction} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
}
  \item{output.cls}{
%%     ~~Describe \code{output.cls} here~~
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
function (input.ds, input.cls = "", column.subset.fraction = "ALL", 
    row.subset.fraction = "ALL", output.ds, output.cls) 
{
    print("Running MSIG.Sample.Dataset...")
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    if (input.cls != "") {
        CLS <- ReadClsFile(file = input.cls)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
    }
    if (column.subset.fraction[1] == "ALL") {
        m2 <- m
        sample.names2 <- sample.names
        if (input.cls != "") {
            class.labels2 <- class.labels
        }
    }
    else {
        ncol <- length(m[1, ])
        column.subset <- sample(x = seq(1, ncol), size = ceiling(column.subset.fraction * 
            ncol), replace = FALSE)
        m2 <- m[, column.subset]
        sample.names2 <- sample.names[column.subset]
        if (input.cls != "") {
            class.labels2 <- class.labels[column.subset]
        }
    }
    if (row.subset.fraction[1] == "ALL") {
        m3 <- m2
        gs.names2 <- gs.names
        gs.descs2 <- gs.descs
    }
    else {
        nrow <- length(m[, 1])
        row.subset <- sample(x = seq(1, nrow), size = ceiling(row.subset.fraction * 
            nrow), replace = FALSE)
        m3 <- m2[row.subset, ]
        gs.names2 <- gs.names[row.subset]
        gs.descs2 <- gs.descs[row.subset]
    }
    V <- data.frame(m3)
    names(V) <- sample.names2
    row.names(V) <- gs.names2
    write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)
    if (input.cls != "") {
        write.cls(class.v = class.labels2, phen = class.phen, 
            filename = output.cls)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
