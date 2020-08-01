\name{MSIG.StJude_to_GCT}
\alias{MSIG.StJude_to_GCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.StJude_to_GCT
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.StJude_to_GCT(input.ds, output.ds, output.cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
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
function (input.ds, output.ds, output.cls) 
{
    print(c("Running MSIG.StJude_to_GCT... on:", input.ds))
    filename <- ds <- read.delim(filename, header = T, sep = "\t", 
        skip = 1, row.names = 1, blank.lines.skip = T, comment.char = "", 
        as.is = T, na.strings = "")
    descs <- ds[, 1]
    ds <- ds[-1]
    row.names <- row.names(ds)
    names <- names(ds)
    return(list(ds = ds, row.names = row.names, descs = descs, 
        names = names))
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
    if (column.subset[1] == "ALL") {
        m2 <- m
        sample.names2 <- sample.names
        if (input.cls != "") {
            class.labels2 <- class.labels
        }
    }
    else {
        m2 <- m[, column.subset]
        if (is.numeric(column.subset[1])) {
            sample.names2 <- sample.names[column.subset]
            if (input.cls != "") {
                class.labels2 <- class.labels[column.subset]
            }
        }
        else {
            locations <- match(column.subset, sample.names)
            sample.names2 <- sample.names[locations]
            if (input.cls != "") {
                class.labels2 <- class.labels[locations]
            }
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
