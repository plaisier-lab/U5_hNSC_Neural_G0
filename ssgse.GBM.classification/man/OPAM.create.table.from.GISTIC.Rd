\name{OPAM.create.table.from.GISTIC}
\alias{OPAM.create.table.from.GISTIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.create.table.from.GISTIC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.create.table.from.GISTIC(gistic.ds, mapping.ds = NULL, output.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gistic.ds}{
%%     ~~Describe \code{gistic.ds} here~~
}
  \item{mapping.ds}{
%%     ~~Describe \code{mapping.ds} here~~
}
  \item{output.file}{
%%     ~~Describe \code{output.file} here~~
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
function (gistic.ds, mapping.ds = NULL, output.file) 
{
    gistic <- read.delim(gistic.ds, header = T, sep = "\t", blank.lines.skip = T, 
        comment.char = "", as.is = T)
    temp <- strsplit(gistic.ds, split = "/")
    s <- length(temp[[1]])
    gistic.file.name <- temp[[1]][s]
    temp <- strsplit(gistic.file.name, split = ".txt")
    gistic.file.prefix <- temp[[1]][1]
    if (!is.null(mapping.ds)) {
        map <- read.delim(mapping.ds, header = T, sep = "\t", 
            blank.lines.skip = T, comment.char = "", as.is = T)
        common.gistic.map <- intersect(colnames(gistic), map[, 
            1])
        loc.gistic <- match(common.gistic.map, colnames(gistic))
        loc.map <- match(common.gistic.map, map[, 1])
        colnames(gistic)[loc.gistic] <- map[loc.map, 2]
    }
    gistic2 <- NULL
    for (i in 1:length(row.names(gistic))) {
        if (regexpr(pattern = "log2", gistic[i, 1]) != -1) 
            next
        str <- paste(substr(gistic[i, 1], 1, 3), gistic[i, 2], 
            sep = "_")
        gistic2 <- rbind(cbind(str, gistic[i, ]), gistic2)
    }
    gistic2[gistic2 == 2] <- 1
    gistic3 <- cbind(gistic2[, 1], gistic2[, 11:length(gistic[1, 
        ])])
    colnames(gistic3) <- c("Sample", colnames(gistic3)[2:length(gistic3[1, 
        ])])
    gistic3 <- t(gistic3)
    write.table(gistic3, file = output.file, quote = F, col.names = F, 
        row.names = T, append = F, sep = "\t")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
