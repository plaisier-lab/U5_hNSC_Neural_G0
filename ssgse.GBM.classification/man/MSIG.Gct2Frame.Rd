\name{MSIG.Gct2Frame}
\alias{MSIG.Gct2Frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Gct2Frame
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Gct2Frame(filename = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
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
function (filename = "NULL") 
{
    ds <- read.delim(filename, header = T, sep = "\t", skip = 2, 
        row.names = 1, blank.lines.skip = T, comment.char = "", 
        as.is = T, na.strings = "")
    descs <- ds[, 1]
    ds <- ds[-1]
    row.names <- row.names(ds)
    names <- names(ds)
    return(list(ds = ds, row.names = row.names, descs = descs, 
        names = names))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
