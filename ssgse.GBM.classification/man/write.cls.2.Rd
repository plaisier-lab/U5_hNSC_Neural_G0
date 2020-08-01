\name{write.cls.2}
\alias{write.cls.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
write.cls.2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.cls.2(class.v, phen, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class.v}{
%%     ~~Describe \code{class.v} here~~
}
  \item{phen}{
%%     ~~Describe \code{phen} here~~
}
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
function (class.v, phen, filename) 
{
    f <- file(filename, "w")
    n <- length(phen)
    l <- length(class.v)
    cat(l, n, "1", "\n", file = f, append = TRUE, sep = " ")
    cat("#", unlist(phen), "\n", file = f, append = TRUE, sep = " ")
    if (is.vector(class.v)) {
        class.v <- phen[class.v]
        cat(class.v, "\n", file = f, append = TRUE, sep = " ")
    }
    else {
        class.list <- matrix(0, nrow = length(class.v[, 1]), 
            ncol = length(class.v[1, ]))
        for (i in 1:length(class.v[, 1])) {
            class.list[i, ] <- unlist(phen[[i]])[class.v[i, ]]
            cat(class.list[i, ], "\n", file = f, append = TRUE, 
                sep = " ")
        }
    }
    close(f)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
