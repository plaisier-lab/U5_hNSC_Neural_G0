\name{MSIG.Relabel.Phenotype}
\alias{MSIG.Relabel.Phenotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Relabel.Phenotype
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Relabel.Phenotype(input.cls, renaming.list, output.cls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{renaming.list}{
%%     ~~Describe \code{renaming.list} here~~
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
function (input.cls, renaming.list, output.cls) 
{
    CLS <- MSIG.ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
    for (k in 1:length(class.phen)) {
        value <- class.phen[k]
        for (i in 1:length(renaming.list)) {
            if (names(renaming.list[i]) == value) {
                class.phen[k] <- renaming.list[i][[1]]
            }
        }
    }
    write.cls(class.v = class.labels, phen = class.phen, filename = output.cls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
