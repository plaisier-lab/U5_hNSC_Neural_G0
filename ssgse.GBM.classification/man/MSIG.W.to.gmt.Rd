\name{MSIG.W.to.gmt}
\alias{MSIG.W.to.gmt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.W.to.gmt
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.W.to.gmt(input.ds, threshold = 1, output.gmt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{output.gmt}{
%%     ~~Describe \code{output.gmt} here~~
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
function (input.ds, threshold = 1, output.gmt) 
{
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    N <- length(m[, 1])
    M <- length(m[1, ])
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    for (i in 1:M) {
        gene.set <- gs.names[m[, i] >= threshold]
        row.header <- paste("Factor_", i, sep = "")
        output.line <- paste(gene.set, sep = "\t", collapse = "\t")
        output.line <- paste(row.header, row.header, output.line, 
            sep = "\t", collapse = "")
        print(paste("factor =", i, " length=", length(gene.set), 
            sep = ""))
        if (i == 1) {
            write(noquote(output.line), file = output.gmt, append = F, 
                ncolumns = length(gene.set) + 2)
        }
        else {
            write(noquote(output.line), file = output.gmt, append = T, 
                ncolumns = length(gene.set) + 2)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
