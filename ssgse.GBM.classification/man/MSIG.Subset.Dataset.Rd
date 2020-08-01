\name{MSIG.Subset.Dataset}
\alias{MSIG.Subset.Dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Subset.Dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Subset.Dataset(input.ds, input.cls = NULL, column.subset = "ALL", column.sel.type = "samples", row.subset = "ALL", output.ds, output.cls = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{column.subset}{
%%     ~~Describe \code{column.subset} here~~
}
  \item{column.sel.type}{
%%     ~~Describe \code{column.sel.type} here~~
}
  \item{row.subset}{
%%     ~~Describe \code{row.subset} here~~
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
function (input.ds, input.cls = NULL, column.subset = "ALL", 
    column.sel.type = "samples", row.subset = "ALL", output.ds, 
    output.cls = NULL) 
{
    print(c("Running MSIG.Subset.Dataset... on GCT file:", input.ds))
    print(c("Running MSIG.Subset.Dataset... on CLS file:", input.cls))
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    if (!is.null(input.cls)) {
        CLS <- MSIG.ReadClsFile(file = input.cls)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
        class.list <- CLS$class.list
    }
    if (column.sel.type == "samples") {
        if (column.subset[1] == "ALL") {
            m2 <- m
            sample.names2 <- sample.names
            if (!is.null(input.cls)) {
                class.labels2 <- class.labels
            }
        }
        else {
            m2 <- m[, column.subset]
            if (is.numeric(column.subset[1])) {
                sample.names2 <- sample.names[column.subset]
                if (!is.null(input.cls)) {
                  class.labels2 <- class.labels[column.subset]
                }
            }
            else {
                locations <- match(column.subset, sample.names)
                sample.names2 <- sample.names[locations]
                if (!is.null(input.cls)) {
                  class.labels2 <- class.labels[locations]
                }
            }
        }
    }
    else if (column.sel.type == "phenotype") {
        locations <- !is.na(match(class.list, column.subset))
        sample.names2 <- sample.names[locations]
        m2 <- m[, locations]
        if (!is.null(input.cls)) {
            class.labels2 <- class.labels[locations]
        }
    }
    if (row.subset[1] == "ALL") {
        m3 <- m2
        gs.names2 <- gs.names
        gs.descs2 <- gs.descs
    }
    else {
        m3 <- m2[row.subset, ]
        if (is.numeric(row.subset[1])) {
            gs.names2 <- gs.names[row.subset]
            gs.descs2 <- gs.descs[row.subset]
        }
        else {
            locations <- match(row.subset, gs.names)
            gs.names2 <- gs.names[locations]
            gs.descs2 <- gs.descs[locations]
        }
    }
    V <- data.frame(m3)
    names(V) <- sample.names2
    row.names(V) <- gs.names2
    write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)
    if (!is.null(input.cls)) {
        write.cls(class.v = class.labels2, phen = class.phen, 
            filename = output.cls)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
