\name{MSIG.Match.and.Select}
\alias{MSIG.Match.and.Select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Match.and.Select
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Match.and.Select(input1.ds, input2.ds, output.ds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input1.ds}{
%%     ~~Describe \code{input1.ds} here~~
}
  \item{input2.ds}{
%%     ~~Describe \code{input2.ds} here~~
}
  \item{output.ds}{
%%     ~~Describe \code{output.ds} here~~
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
function (input1.ds, input2.ds, output.ds) 
{
    print(c("Running MSIG.Match.and.Select... on: ", input1.ds, 
        " ", input2.ds))
    dataset1 <- MSIG.Gct2Frame(filename = input1.ds)
    m1 <- data.matrix(dataset1$ds)
    gs.names1 <- dataset1$row.names
    gs.descs1 <- dataset1$descs
    sample.names1 <- dataset1$names
    dataset2 <- MSIG.Gct2Frame(filename = input2.ds)
    m2 <- data.matrix(dataset2$ds)
    gs.names2 <- dataset2$row.names
    gs.descs2 <- dataset2$descs
    sample.names2 <- dataset2$names
    gs.names3 <- intersect(gs.names1, gs.names2)
    locations2 <- match(gs.names3, gs.names2, nomatch = 0)
    gs.names2 <- gs.names2[locations2]
    gs.descs2 <- gs.descs2[locations2]
    m2 <- m2[locations2, ]
    V <- data.frame(m2)
    names(V) <- sample.names2
    row.names(V) <- gs.names2
    write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
