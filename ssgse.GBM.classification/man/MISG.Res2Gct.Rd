\name{MISG.Res2Gct}
\alias{MISG.Res2Gct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MISG.Res2Gct
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MISG.Res2Gct(res.file, gct.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.file}{
%%     ~~Describe \code{res.file} here~~
}
  \item{gct.file}{
%%     ~~Describe \code{gct.file} here~~
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
function (res.file, gct.file) 
{
    dataset <- MSIG.Res2Frame(filename = res.file)
    A <- dataset$ds
    row.names(A) <- dataset$row.names
    colnames(A) <- dataset$names
    descs <- dataset$descs
    write.gct(gct.data.frame = A, descs = descs, filename = gct.file)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
