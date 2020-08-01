\name{MSIG.Match.Tables}
\alias{MSIG.Match.Tables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Match.Tables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Match.Tables(input1.txt, input2.txt, output.txt, match.rows = T, mode = "intersection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input1.txt}{
%%     ~~Describe \code{input1.txt} here~~
}
  \item{input2.txt}{
%%     ~~Describe \code{input2.txt} here~~
}
  \item{output.txt}{
%%     ~~Describe \code{output.txt} here~~
}
  \item{match.rows}{
%%     ~~Describe \code{match.rows} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
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
function (input1.txt, input2.txt, output.txt, match.rows = T, 
    mode = "intersection") 
{
    tab1 <- read.delim(input1.txt, header = T, row.names = 1, 
        sep = "\t", skip = 0, blank.lines.skip = T, comment.char = "", 
        as.is = T)
    tab2 <- read.delim(input2.txt, header = T, row.names = 1, 
        sep = "\t", skip = 0, blank.lines.skip = T, comment.char = "", 
        as.is = T)
    if (match.rows == T) {
        entries.tab1 <- row.names(tab1)
        N1 <- length(entries.tab1)
        M1 <- length(names(tab1))
        entries.tab2 <- row.names(tab2)
        N2 <- length(entries.tab2)
        M2 <- length(names(tab2))
        print(c("table 1:", N1, " entries"))
        print(c("table 1:", M1, " columns"))
        print(c("table 2:", N2, " entries"))
        print(c("table 2:", M2, " columns"))
        if (mode == "intersection") {
            entries.tab3 <- intersect(entries.tab1, entries.tab2)
            N3 <- length(entries.tab3)
            M3 <- M1 + M2
            locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
            tab1 <- tab1[locations1, ]
            locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)
            tab2 <- tab2[locations2, ]
            tab3 <- cbind(tab1, tab2)
        }
        else if (mode == "union") {
            entries.tab3 <- union(entries.tab1, entries.tab2)
            N3 <- length(entries.tab3)
            M3 <- M1 + M2
            m3 <- matrix("UNK", nrow = N3, ncol = M3)
            tab3 <- matrix("NA", nrow = N3, ncol = M3)
            locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
            locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)
            for (i in 1:N3) {
                if (locations1[i] != 0) {
                  tab3[i, 1:M1] <- as.matrix(tab1[locations1[i], 
                    ])
                }
                if (locations2[i] != 0) {
                  tab3[i, (M1 + 1):M3] <- as.matrix(tab2[locations2[i], 
                    ])
                }
            }
            colnames(tab3) <- c(colnames(tab1), colnames(tab2))
        }
        else if (mode == "match.to.first") {
            N3 <- length(entries.tab1)
            M3 <- M1 + M2
            m3 <- matrix("UNK", nrow = N3, ncol = M3)
            tab3 <- matrix("NA", nrow = N3, ncol = M3)
            locations2 <- match(entries.tab1, entries.tab2, nomatch = 0)
            for (i in 1:N3) {
                tab3[i, 1:M1] <- as.matrix(tab1[i, ])
                if (locations2[i] != 0) {
                  tab3[i, (M1 + 1):M3] <- as.matrix(tab2[locations2[i], 
                    ])
                }
            }
            colnames(tab3) <- c(colnames(tab1), colnames(tab2))
            row.names(tab3) <- entries.tab1
        }
        else {
            stop(c("unknown mode", mode))
        }
        print(c("table 3:", N3, " entries"))
        print(c("table 3:", M3, " columns"))
        col.names <- paste(colnames(tab3), collapse = "\t")
        col.names <- paste("SAMPLE", col.names, sep = "\t")
        write(noquote(col.names), file = output.txt, append = F, 
            ncolumns = length(col.names))
        write.table(tab3, file = output.txt, quote = F, col.names = F, 
            row.names = T, append = T, sep = "\t")
    }
    else {
        entries.tab1 <- colnames(tab1)
        M1 <- length(entries.tab1)
        N1 <- length(row.names(tab1))
        entries.tab2 <- colnames(tab2)
        M2 <- length(entries.tab2)
        N2 <- length(row.names(tab2))
        print(c("table 1:", N1, " rows"))
        print(c("table 1:", M1, " columns"))
        print(c("table 2:", N2, " rows"))
        print(c("table 2:", M2, " columns"))
        if (mode == "intersection") {
            entries.tab3 <- intersect(entries.tab1, entries.tab2)
            M3 <- length(entries.tab3)
            N3 <- N1 + N2
            locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
            tab1 <- tab1[, locations1]
            locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)
            tab2 <- tab2[, locations2]
            tab3 <- rbind(tab1, tab2)
        }
        else if (mode == "union") {
            entries.tab3 <- union(entries.tab1, entries.tab2)
            M3 <- length(entries.tab3)
            N3 <- N1 + N2
            m3 <- matrix("UNK", nrow = N3, ncol = M3)
            tab3 <- matrix("NA", nrow = N3, ncol = M3)
            locations1 <- match(entries.tab3, entries.tab1, nomatch = 0)
            locations2 <- match(entries.tab3, entries.tab2, nomatch = 0)
            for (i in 1:M3) {
                if (locations1[i] != 0) {
                  tab3[1:N1, i] <- as.matrix(tab1[, locations1[i]])
                }
                if (locations2[i] != 0) {
                  tab3[(N1 + 1):N3, i] <- as.matrix(tab2[, locations2[i]])
                }
            }
            row.names(tab3) <- c(row.names(tab1), row.names(tab2))
        }
        else if (mode == "match.to.first") {
            M3 <- length(entries.tab1)
            N3 <- N1 + N2
            m3 <- matrix("UNK", nrow = N3, ncol = M3)
            tab3 <- matrix("NA", nrow = N3, ncol = M3)
            locations2 <- match(entries.tab1, entries.tab2, nomatch = 0)
            for (i in 1:M3) {
                tab3[1:N1, i] <- as.matrix(tab1[, i])
                if (locations2[i] != 0) {
                  tab3[(N1 + 1):N3, i] <- as.matrix(tab2[, locations2[i]])
                }
            }
            row.names(tab3) <- c(row.names(tab1), row.names(tab2))
            colnames(tab3) <- entries.tab1
        }
        else {
            stop(c("unknown mode", mode))
        }
        print(c("table 3:", N3, " entries"))
        print(c("table 3:", M3, " columns"))
        col.names <- paste(colnames(tab3), collapse = "\t")
        col.names <- paste("SAMPLE", col.names, sep = "\t")
        write(noquote(col.names), file = output.txt, append = F, 
            ncolumns = length(col.names))
        write.table(tab3, file = output.txt, quote = F, col.names = F, 
            row.names = T, append = T, sep = "\t")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
