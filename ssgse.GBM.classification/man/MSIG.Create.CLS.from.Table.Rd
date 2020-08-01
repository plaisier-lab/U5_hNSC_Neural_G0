\name{MSIG.Create.CLS.from.Table}
\alias{MSIG.Create.CLS.from.Table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Create.CLS.from.Table
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Create.CLS.from.Table(file.gct, table.txt, output.gct = NULL, output.cls, sort.by = NULL, then.sort.by = NULL, only.matches = T, rename = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.gct}{
%%     ~~Describe \code{file.gct} here~~
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
  \item{sort.by}{
%%     ~~Describe \code{sort.by} here~~
}
  \item{then.sort.by}{
%%     ~~Describe \code{then.sort.by} here~~
}
  \item{only.matches}{
%%     ~~Describe \code{only.matches} here~~
}
  \item{rename}{
%%     ~~Describe \code{rename} here~~
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
function (file.gct, table.txt, output.gct = NULL, output.cls, 
    sort.by = NULL, then.sort.by = NULL, only.matches = T, rename = T) 
{
    dataset1 <- MSIG.Gct2Frame(filename = file.gct)
    m <- data.matrix(dataset1$ds)
    gene.names <- dataset1$row.names
    gene.decs <- dataset1$descs
    sample.names.gct <- dataset1$names
    Ns <- length(sample.names.gct)
    tab <- read.delim(table.txt, header = T, row.names = 1, sep = "\t", 
        skip = 0, blank.lines.skip = T, comment.char = "", as.is = T)
    sample.names.tab <- row.names(tab)
    phen.names <- names(tab)
    print(c("phen names:", phen.names))
    locs <- match(sample.names.gct, sample.names.tab, nomatch = NULL)
    print(c("locs:", locs))
    tab2 <- matrix(0, nrow = length(locs), ncol = length(tab[1, 
        ]))
    tab2 <- tab[locs, ]
    row.names.tab2 <- row.names(tab)[locs]
    cls.table <- as.matrix(t(tab2))
    print(c("tab2", tab2))
    print(c("sample.names.gct:", sample.names.gct))
    print(c("sample.names.tab:", sample.names.tab))
    print(c("Total samples GCT file:", length(sample.names.gct)))
    print(c("Total samples table file:", length(sample.names.tab)))
    if (only.matches) {
        temp <- substr(row.names.tab2, 1, 2)
        locs <- !is.na(temp)
        print(c("locs:", locs))
        if (is.vector(cls.table)) {
            cls.table <- cls.table[locs]
            print(c("Number of matches:", length(cls.table)))
        }
        else {
            cls.table <- cls.table[, locs]
            print(c("Number of matches:", length(cls.table[1, 
                ])))
        }
        if (!is.null(output.gct)) {
            m <- m[, locs]
            sample.names.gct <- sample.names.gct[locs]
        }
    }
    print(c("matching.names: (after)", sample.names.gct))
    if (is.vector(cls.table)) {
        name <- phen.names
        for (j in 1:length(cls.table)) {
            if (rename == T) {
                if (is.na(cls.table[j])) {
                  cls.table[j] <- "UNK"
                }
                else if (cls.table[j] == 1) {
                  cls.table[j] <- name
                }
                else {
                  cls.table[j] <- "WT"
                }
            }
        }
        class.phen <- unique(cls.table)
        n <- length(class.phen)
        l <- length(cls.table)
        cat(l, n, "1", "\n", file = output.cls, append = FALSE, 
            sep = " ")
        cat("#", class.phen, "\n", file = output.cls, append = TRUE, 
            sep = " ")
        cat(cls.table, "\n", file = output.cls, append = TRUE, 
            sep = " ")
    }
    else {
        class.phen <- NULL
        for (i in 1:length(cls.table[, 1])) {
            name <- row.names(cls.table)[i]
            for (j in 1:length(cls.table[1, ])) {
                if (rename == T) {
                  if (is.na(cls.table[i, j])) {
                    cls.table[i, j] <- "UNK"
                  }
                  else if (cls.table[i, j] == 1) {
                    cls.table[i, j] <- name
                  }
                  else {
                    cls.table[i, j] <- "WT"
                  }
                }
            }
            class.phen[[i]] <- unique(cls.table[i, ])
        }
        n <- length(class.phen)
        l <- length(cls.table[1, ])
        cat(l, n, "1", "\n", file = output.cls, append = FALSE, 
            sep = " ")
        cat("#", unlist(class.phen), "\n", file = output.cls, 
            append = TRUE, sep = " ")
        for (i in 1:length(cls.table[, 1])) {
            cat(cls.table[i, ], "\n", file = output.cls, append = TRUE, 
                sep = " ")
        }
    }
    if (!is.null(sort.by)) {
        print(cls.table)
        if (!is.null(then.sort.by)) {
            new.order <- order(cls.table[sort.by, ], cls.table[then.sort.by, 
                ], decreasing = F)
        }
        else {
            new.order <- order(cls.table[sort.by, ], decreasing = F)
        }
        if (is.vector(cls.table[, new.order])) {
            cls.table <- cls.table[new.order]
        }
        else {
            cls.table <- cls.table[, new.order]
        }
        print(cls.table)
        print(dim(cls.table))
        sample.names.gct <- sample.names.gct[new.order]
        if (!is.null(output.gct)) {
            m <- m[, new.order]
        }
    }
    if (!is.null(output.gct)) {
        V <- data.frame(m)
        names(V) <- sample.names.gct
        row.names(V) <- gene.names
        write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
