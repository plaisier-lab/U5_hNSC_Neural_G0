\name{CNMF.write.gct}
\alias{CNMF.write.gct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CNMF.write.gct
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNMF.write.gct(gct, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gct}{
%%     ~~Describe \code{gct} here~~
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
function (gct, filename) 
{
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct)[1], "     ", dim(gct)[2], "\n", file = f, append = TRUE, 
        sep = "")
    cat("Name", "       ", file = f, append = TRUE, sep = "")
    cat("Description", file = f, append = TRUE, sep = "")
    names <- names(gct)
    cat("       ", names[1], file = f, append = TRUE, sep = "")
    for (j in 2:length(names)) {
        cat("   ", names[j], file = f, append = TRUE, sep = "")
    }
    cat("\n", file = f, append = TRUE, sep = "")
    oldWarn <- options(warn = -1)
    m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] + 2)
    m[, 1] <- row.names(gct)
    m[, 2] <- row.names(gct)
    index <- 3
    for (i in 1:dim(gct)[2]) {
        m[, index] <- gct[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "      ", 
        eol = "\n", col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)
    return(gct)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
