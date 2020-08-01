\name{MSIG.Define.Dataset.from.Table}
\alias{MSIG.Define.Dataset.from.Table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Define.Dataset.from.Table
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Define.Dataset.from.Table(input.gct, table.txt, output.gct, output.cls, prefix_entries = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.gct}{
%%     ~~Describe \code{input.gct} here~~
}
  \item{table.txt}{
%%     ~~Describe \code{table.txt} here~~
}
  \item{output.gct}{
%%     ~~Describe \code{output.gct} here~~
}
  \item{output.cls}{
%%     ~~Describe \code{output.cls} here~~
}
  \item{prefix_entries}{
%%     ~~Describe \code{prefix_entries} here~~
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
function (input.gct, table.txt, output.gct, output.cls, prefix_entries = F) 
{
    library(RColorBrewer)
    dataset1 <- MSIG.Gct2Frame(filename = input.gct)
    m <- data.matrix(dataset1$ds)
    gene.names <- dataset1$row.names
    gene.decs <- dataset1$descs
    sample.names.gct <- dataset1$names
    Ns <- length(sample.names.gct)
    tab <- read.delim(table.txt, header = T, row.names = 1, sep = "\t", 
        skip = 0, blank.lines.skip = T, comment.char = "", as.is = T)
    sample.names.tab <- row.names(tab)
    phen.names <- names(tab)
    overlap <- intersect(sample.names.tab, sample.names.gct)
    print("sample names GCT")
    print(sample.names.gct)
    print("sample names TAB")
    print(sample.names.tab)
    locs.gct <- match(overlap, sample.names.gct)
    print(match(sample.names.tab, sample.names.gct))
    print(match(sample.names.gct, sample.names.tab))
    locs.tab <- match(overlap, sample.names.tab)
    print(locs.tab)
    print(c("GCT matching set (", length(locs.gct), " samples):", 
        sample.names.gct[locs.gct]))
    print(c("TAB matching set (", length(overlap), " samples):", 
        sample.names.tab[locs.tab]))
    print(c("overlap set (", length(overlap), " samples):", overlap))
    m2 <- m[, locs.gct]
    sample.names.gct <- sample.names.gct[locs.gct]
    sample.names.tab <- sample.names.tab[locs.tab]
    cls.table <- t(tab[locs.tab, ])
    if (prefix_entries == TRUE) {
        for (i in 1:length(cls.table[, 1])) {
            cls.table[i, ] <- paste(colnames(tab)[i], tab[, i], 
                sep = ".")
        }
    }
    if (!is.null(output.gct)) {
        V <- data.frame(m2)
        names(V) <- sample.names.gct
        row.names(V) <- gene.names
        write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
    }
    class.phen <- unique(cls.table)
    n <- length(class.phen)
    l <- length(cls.table[1, ])
    col.list <- c(brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, 
        name = "Dark2"), brewer.pal(n = 7, name = "Set1"), brewer.pal(n = 8, 
        name = "Accent"), brewer.pal(n = 10, name = "Spectral"), 
        brewer.pal(n = 8, name = "Set3"), brewer.pal(n = 8, name = "BuGn"), 
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), 
        brewer.pal(n = 7, name = "Set1"), brewer.pal(n = 8, name = "Accent"), 
        brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, 
            name = "Set3"), brewer.pal(n = 8, name = "BuGn"), 
        brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, name = "Dark2"), 
        brewer.pal(n = 7, name = "Set1"), brewer.pal(n = 8, name = "Accent"), 
        brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, 
            name = "Set3"), brewer.pal(n = 8, name = "BuGn"))
    num <- 0
    class.order.list <- NULL
    for (i in 1:length(cls.table[, 1])) {
        num <- num + length(unique(cls.table[i, ]))
        class.order.list <- c(class.order.list, unique(cls.table[i, 
            ]))
    }
    phen.names.string <- paste("phen.names:", paste(phen.names, 
        collapse = " "), sep = " ")
    sig.col <- col.list[1:num]
    col.phen.string <- paste("col.phen:", paste(sig.col, collapse = " "), 
        sep = " ")
    cat(paste(l, num, length(cls.table[, 1]), phen.names.string, 
        col.phen.string, sep = " "), "\n", file = output.cls, 
        append = FALSE, sep = "")
    cat("# ", paste(class.order.list, collapse = " "), "\n", 
        file = output.cls, append = TRUE, sep = "")
    for (i in 1:length(cls.table[, 1])) {
        cat(paste(cls.table[i, ], collapse = " "), "\n", file = output.cls, 
            append = TRUE, sep = "")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
