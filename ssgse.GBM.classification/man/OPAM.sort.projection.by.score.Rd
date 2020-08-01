\name{OPAM.sort.projection.by.score}
\alias{OPAM.sort.projection.by.score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.sort.projection.by.score
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.sort.projection.by.score(input.ds, input.cls, results.dir, normalize.score = T, normalization.type = "zero.one", model, user.colors = NA, decreasing.order = T, output.dataset = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{results.dir}{
%%     ~~Describe \code{results.dir} here~~
}
  \item{normalize.score}{
%%     ~~Describe \code{normalize.score} here~~
}
  \item{normalization.type}{
%%     ~~Describe \code{normalization.type} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{user.colors}{
%%     ~~Describe \code{user.colors} here~~
}
  \item{decreasing.order}{
%%     ~~Describe \code{decreasing.order} here~~
}
  \item{output.dataset}{
%%     ~~Describe \code{output.dataset} here~~
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
function (input.ds, input.cls, results.dir, normalize.score = T, 
    normalization.type = "zero.one", model, user.colors = NA, 
    decreasing.order = T, output.dataset = NA) 
{
    library(gtools)
    library(verification)
    library(ROCR)
    library(MASS)
    library(RColorBrewer)
    library(heatmap.plus)
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    model.names <- dataset$row.names
    model.descs <- dataset$descs
    Ns <- length(m[1, ])
    dim(m)
    sample.names <- dataset$names
    n.models <- length(m[, 1])
    temp <- strsplit(input.ds, split = "/")
    s <- length(temp[[1]])
    test.file.name <- temp[[1]][s]
    temp <- strsplit(test.file.name, split = ".gct")
    test.file.prefix <- temp[[1]][1]
    char.res <- 0.013 * n.models + 0.65
    if (normalize.score == T) {
        if (normalization.type == "zero.one") {
            for (i in 1:n.models) {
                m[i, ] <- (m[i, ] - min(m[i, ]))/(max(m[i, ]) - 
                  min(m[i, ]))
            }
        }
        else if (normalization.type == "z.score") {
            for (i in 1:n.models) {
                m[i, ] <- (m[i, ] - mean(m[i, ]))/sd(m[i, ])
            }
        }
        else if (normalization.type == "r.z.score") {
            for (i in 1:n.models) {
                m[i, ] <- (m[i, ] - median(m[i, ]))/mad(m[i, 
                  ])
            }
        }
    }
    CLS <- MSIG.ReadPhenFile(file = input.cls)
    cls.labels <- CLS$class.v
    cls.phen <- CLS$phen
    cls.list <- CLS$class.list
    if (is.vector(cls.labels)) {
        n.phen <- 1
    }
    else {
        n.phen <- length(cls.labels[, 1])
    }
    if (!is.na(user.colors)) {
        c.test <- user.colors
    }
    else {
        if (!is.null(CLS$col.phen)) {
            c.test <- CLS$col.phen
        }
        else {
            c.test <- c(brewer.pal(n = 7, name = "Set2"), brewer.pal(n = 7, 
                name = "Dark2"), brewer.pal(n = 7, name = "Set1"), 
                brewer.pal(n = 8, name = "Accent"), brewer.pal(n = 9, 
                  name = "Spectral"), brewer.pal(n = 8, name = "Set3"), 
                brewer.pal(n = 8, name = "BuGn"), brewer.pal(n = 7, 
                  name = "Set2"), brewer.pal(n = 7, name = "Dark2"), 
                brewer.pal(n = 7, name = "Set1"), brewer.pal(n = 8, 
                  name = "Accent"), brewer.pal(n = 10, name = "Spectral"), 
                brewer.pal(n = 8, name = "Set3"), brewer.pal(n = 8, 
                  name = "BuGn"), brewer.pal(n = 7, name = "Set2"), 
                brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, 
                  name = "Set1"), brewer.pal(n = 8, name = "Accent"), 
                brewer.pal(n = 10, name = "Spectral"), brewer.pal(n = 8, 
                  name = "Set3"), brewer.pal(n = 8, name = "BuGn"))
        }
    }
    if (!is.null(CLS$phen.names)) {
        phen.names <- CLS$phen.names
    }
    else {
        phen.names <- "NA"
    }
    cls.phen.index <- unlist(cls.phen)
    cls.phen.colors <- c.test[1:length(cls.phen.index)]
    n.classes <- vector(length = n.phen, mode = "numeric")
    if (n.phen == 1) {
        max.classes <- length(cls.phen)
        n.classes[1] <- max.classes
    }
    else {
        max.classes <- max(unlist(lapply(cls.phen, FUN = length)))
        for (i in 1:n.phen) {
            n.classes[i] <- length(cls.phen[[i]])
        }
    }
    filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    loc <- match(model, model.names)
    print(c("loc:", loc))
    s.order <- order(m[loc, ], decreasing = decreasing.order)
    m2 <- m[, s.order]
    sample.names2 <- sample.names[s.order]
    if (is.vector(cls.labels)) {
        cls.labels2 <- cls.labels[s.order]
        cls.list2 <- cls.list[s.order]
    }
    else {
        cls.labels2 <- cls.labels[, s.order]
        cls.list2 <- cls.list[, s.order]
    }
    cls.phen2 <- NULL
    if (is.vector(cls.labels)) {
        classes <- unique(cls.list2)
        cls.phen2 <- classes
        cls.labels2 <- match(cls.list2, cls.phen2)
    }
    else {
        for (kk in 1:length(cls.list2[, 1])) {
            classes <- unique(cls.list2[kk, ])
            cls.phen2[[kk]] <- classes
            cls.labels2[kk, ] <- match(cls.list2[kk, ], cls.phen2[[kk]])
        }
    }
    correl <- cor(t(m2))[, loc]
    m.order <- order(correl, decreasing = T)
    correl2 <- correl[m.order]
    m2 <- m2[m.order, ]
    model.names2 <- model.names[m.order]
    model.descs2 <- paste(model.descs[m.order], signif(correl2, 
        digits = 3))
    phen.list <- unlist(cls.phen2)
    colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
    MSIG.HeatMapPlot.7(V = m2, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = cls.labels2, col.classes = cls.phen2, phen.cmap = colors.list, 
        phen.names = phen.names, col.names = sample.names2, main = " ", 
        xlab = "  ", ylab = "  ", row.norm = T, cmap.type = 3, 
        char.rescale = 1, legend = T)
    dev.off()
    if (!is.na(output.dataset)) {
        V.GCT <- m2
        colnames(V.GCT) <- sample.names2
        row.names(V.GCT) <- model.names2
        write.gct(gct.data.frame = V.GCT, descs = model.descs2, 
            filename = output.dataset)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
