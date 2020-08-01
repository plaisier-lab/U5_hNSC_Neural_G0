\name{OPAM.write.param.line}
\alias{OPAM.write.param.line}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.write.param.line
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.write.param.line(param, param.name, file, append = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{param.name}{
%%     ~~Describe \code{param.name} here~~
}
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{append}{
%%     ~~Describe \code{append} here~~
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
function (param, param.name, file, append = T) 
{
    if (typeof(param) == "character") {
        output.line <- paste("'", param, "'", sep = "", collapse = ",")
    }
    else {
        output.line <- paste(param, collapse = ",")
    }
    output.line
    if (length(param) > 1) {
        output.line <- paste(param.name, paste("c(", noquote(output.line), 
            ")", sep = ""), sep = "\t")
    }
    else {
        output.line <- paste(param.name, noquote(output.line), 
            sep = "\t")
    }
    output.line
    write(output.line, file = file, append = append, ncolumns = length(param) + 
        1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
