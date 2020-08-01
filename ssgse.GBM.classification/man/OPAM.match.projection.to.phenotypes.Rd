\name{OPAM.match.projection.to.phenotypes}
\alias{OPAM.match.projection.to.phenotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.match.projection.to.phenotypes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.match.projection.to.phenotypes(input.ds, input.cls, results.dir, normalize.score = T, normalization.type = "zero.one", markers.num = 5, user.colors = NA, markers.metric = "ROC", markers.file = NULL, sort.phenotypes = T, sort.decreasing = T, sort.expression = T, sort.decreasing.genes = T, legend = T, char.res = 1, only.up = F, cmap.type = 3, row.norm = T)
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
  \item{markers.num}{
%%     ~~Describe \code{markers.num} here~~
}
  \item{user.colors}{
%%     ~~Describe \code{user.colors} here~~
}
  \item{markers.metric}{
%%     ~~Describe \code{markers.metric} here~~
}
  \item{markers.file}{
%%     ~~Describe \code{markers.file} here~~
}
  \item{sort.phenotypes}{
%%     ~~Describe \code{sort.phenotypes} here~~
}
  \item{sort.decreasing}{
%%     ~~Describe \code{sort.decreasing} here~~
}
  \item{sort.expression}{
%%     ~~Describe \code{sort.expression} here~~
}
  \item{sort.decreasing.genes}{
%%     ~~Describe \code{sort.decreasing.genes} here~~
}
  \item{legend}{
%%     ~~Describe \code{legend} here~~
}
  \item{char.res}{
%%     ~~Describe \code{char.res} here~~
}
  \item{only.up}{
%%     ~~Describe \code{only.up} here~~
}
  \item{cmap.type}{
%%     ~~Describe \code{cmap.type} here~~
}
  \item{row.norm}{
%%     ~~Describe \code{row.norm} here~~
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
    normalization.type = "zero.one", markers.num = 5, user.colors = NA, 
    markers.metric = "ROC", markers.file = NULL, sort.phenotypes = T, 
    sort.decreasing = T, sort.expression = T, sort.decreasing.genes = T, 
    legend = T, char.res = 1, only.up = F, cmap.type = 3, row.norm = T) 
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
    for (i in 1:length(m[, 1])) {
        if (sd(m[i, ]) == 0) {
            val <- m[i, 1]
            m[i, ] <- m[i, ] + runif(n = Ns, min = val - 0.001, 
                max = val + 0.001)
        }
    }
    dim(m)
    sample.names <- dataset$names
    n.models <- length(m[, 1])
    temp <- strsplit(input.ds, split = "/")
    s <- length(temp[[1]])
    test.file.name <- temp[[1]][s]
    temp <- strsplit(test.file.name, split = ".gct")
    test.file.prefix <- temp[[1]][1]
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
            c.test <- c(brewer.pal(n = 7, name = "Dark2"), brewer.pal(n = 7, 
                name = "Set2"), brewer.pal(n = 7, name = "Set1"), 
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
    x <- rbind(sample.names, cls.list, cls.labels)
    print("before loop")
    print(x)
    print(cls.phen)
    print(phen.names)
    filename <- paste(results.dir, test.file.prefix, ".PHEN.MARKERS.", 
        markers.metric, sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    for (k.phen in 1:n.phen) {
        if (is.vector(cls.labels)) {
            k.phen.labels <- cls.labels
            k.phen.list <- cls.list
        }
        else {
            k.phen.labels <- as.vector(cls.labels[k.phen, ])
            k.phen.list <- as.vector(cls.list[k.phen, ])
        }
        if (sort.expression == T) {
            phen.index <- order(k.phen.labels, decreasing = sort.decreasing)
        }
        else {
            phen.index <- seq(1, length(k.phen.labels))
        }
        if (is.vector(cls.labels)) {
            cls.labels2 <- cls.labels[phen.index]
            cls.list2 <- cls.list[phen.index]
        }
        else {
            cls.labels2 <- cls.labels[, phen.index]
            cls.list2 <- cls.list[, phen.index]
        }
        k.phen.labels <- k.phen.labels[phen.index]
        k.phen.list <- k.phen.list[phen.index]
        sample.names2 <- sample.names[phen.index]
        m2 <- m[, phen.index]
        x <- rbind(sample.names2, cls.list2, cls.labels2)
        print(paste("inside loop phen=", k.phen))
        print(x)
        print(cls.phen)
        print(phen.names)
        if (is.vector(cls.labels2)) {
            classes <- unique(cls.list2)
        }
        else {
            classes <- unique(cls.list2[k.phen, ])
        }
        if (length(classes) > 2) {
            k.only.up <- T
        }
        else {
            k.only.up <- only.up
        }
        if (length(classes) == 2) 
            classes <- classes[1]
        markers <- NULL
        markers.descs <- NULL
        metric.list <- NULL
        p.val.list <- NULL
        k.class <- NULL
        for (k in classes) {
            if (is.vector(cls.labels2)) {
                bin.class <- ifelse(cls.list2 == k, 0, 1)
            }
            else {
                bin.class <- ifelse(cls.list2[k.phen, ] == k, 
                  0, 1)
            }
            if (markers.metric == "T.TEST") {
                metric <- vector(length = n.models, mode = "numeric")
                p.val <- vector(length = n.models, mode = "numeric")
                for (i in 1:n.models) {
                  temp <- split(m2[i, ], bin.class)
                  x <- temp[[1]]
                  y <- temp[[2]]
                  metric[i] <- signif(t.test(x = x, y = y)$statistic, 
                    digits = 3)
                  p.val[i] <- signif(t.test(x = x, y = y)$p.value, 
                    digits = 3)
                }
            }
            else if (markers.metric == "ROC") {
                bin.class <- ifelse(bin.class == 1, 0, 1)
                metric <- vector(length = n.models, mode = "numeric")
                p.val <- vector(length = n.models, mode = "numeric")
                for (i in 1:n.models) {
                  m.score <- m2[i, ]
                  m.score.norm <- (m.score - min(m.score))/(max(m.score) - 
                    min(m.score))
                  perf.auc <- roc.area(bin.class, m.score.norm)
                  metric[i] <- signif(perf.auc$A, digits = 3)
                  p.val[i] <- signif(perf.auc$p.value, digits = 3)
                }
            }
            else if (markers.metric == "MEAN.DIFF") {
                bin.class <- ifelse(bin.class == 1, 0, 1)
                metric <- vector(length = n.models, mode = "numeric")
                p.val <- vector(length = n.models, mode = "numeric")
                for (i in 1:n.models) {
                  temp <- split(m2[i, ], bin.class)
                  x <- temp[[1]]
                  y <- temp[[2]]
                  metric[i] <- signif(mean(x) - mean(y), digits = 3)
                  p.val[i] <- signif(t.test(x = x, y = y)$p.value, 
                    digits = 3)
                }
            }
            if (is.na(sort.decreasing.genes)) {
                metric.order <- seq(1, length(metric))
            }
            else {
                metric.order <- order(metric, decreasing = sort.decreasing.genes)
            }
            if (only.up == TRUE) {
                if (length(classes) == 2) {
                  k.markers.num <- ifelse(markers.num > n.models, 
                    n.models, markers.num)
                }
                else {
                  k.markers.num <- ifelse(length(classes) * markers.num > 
                    n.models, floor(n.models/length(classes)), 
                    markers.num)
                }
                markers <- c(markers, model.names[metric.order][1:k.markers.num])
                markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num])
                metric.list <- c(metric.list, metric[metric.order][1:k.markers.num])
                p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num])
                k.class <- c(k.class, rep(k, k.markers.num))
            }
            else {
                k.markers.num <- ifelse(length(classes) * markers.num > 
                  n.models, floor(n.models/length(classes)), 
                  markers.num)
                markers <- c(markers, model.names[metric.order][1:k.markers.num], 
                  model.names[metric.order][(length(model.names) - 
                    k.markers.num + 1):length(model.names)])
                markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num], 
                  model.descs[metric.order][(length(model.names) - 
                    k.markers.num + 1):length(model.names)])
                metric.list <- c(metric.list, metric[metric.order][1:k.markers.num], 
                  metric[metric.order][(length(model.names) - 
                    k.markers.num + 1):length(model.names)])
                p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num], 
                  p.val[metric.order][(length(model.names) - 
                    k.markers.num + 1):length(model.names)])
                k.class <- c(k.class, rep(k, k.markers.num), 
                  rep(paste("not", k), k.markers.num))
            }
        }
        V3 <- m2[markers, ]
        model.descs2 <- paste(metric.list, p.val.list, k.class, 
            markers.descs)
        height <- ifelse(length(markers) + n.phen >= 9, 10, (length(markers) + 
            n.phen) * 0.44 + 5)
        if (sort.expression == T) {
            for (j in unique(k.phen.labels)) {
                V4 <- V3[, k.phen.labels == j]
                sn <- sample.names2[k.phen.labels == j]
                if (is.vector(cls.labels)) {
                  clab <- cls.labels2[k.phen.labels == j]
                  clis <- cls.list2[k.phen.labels == j]
                }
                else {
                  clab <- cls.labels2[, k.phen.labels == j]
                  clis <- cls.list2[, k.phen.labels == j]
                }
                dist.matrix <- dist(t(V4))
                HC <- hclust(dist.matrix, method = "complete")
                V4 <- V4[, HC$order]
                sn <- sn[HC$order]
                if (is.vector(cls.labels2)) {
                  clab <- clab[HC$order]
                  clis <- clis[HC$order]
                }
                else {
                  clab <- clab[, HC$order]
                  clis <- clis[, HC$order]
                }
                V3[, k.phen.labels == j] <- V4
                sample.names2[k.phen.labels == j] <- sn
                if (is.vector(cls.labels2)) {
                  cls.labels2[k.phen.labels == j] <- clab
                  cls.list2[k.phen.labels == j] <- clis
                }
                else {
                  cls.labels2[, k.phen.labels == j] <- clab
                  cls.list2[, k.phen.labels == j] <- clis
                }
            }
        }
        x <- rbind(sample.names2, cls.list2, cls.labels2)
        print(paste("inside loop after in-class sort phen=", 
            k.phen))
        print(x)
        print(cls.phen)
        print(phen.names)
        cls.phen2 <- NULL
        if (is.vector(cls.labels2)) {
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
        x <- rbind(sample.names2, cls.list2, cls.labels2)
        print(paste("inside loop after cls.phen renorm phen=", 
            k.phen))
        print(cls.phen2)
        print(phen.names)
        library(gmodels)
        if (!is.vector(cls.labels2)) {
            if (sort.phenotypes == T) {
                phen.score <- vector(length = n.phen, mode = "numeric")
                for (k.lab in 1:n.phen) {
                  tab <- table(as.vector(cls.list2[k.lab, ]), 
                    k.phen.list)
                  print(tab)
                  CT <- CrossTable(tab, chisq = T)
                  phen.score[k.lab] <- CT$chisq$p.value
                  print(phen.score[k.lab])
                }
                phen.order <- order(phen.score, decreasing = T)
                print(phen.order)
                cls.labels2 <- cls.labels2[phen.order, ]
                cls.phen2 <- cls.phen2[phen.order]
                phen.names2 <- phen.names[phen.order]
                main.string <- paste(test.file.prefix, " - ", 
                  phen.names2[n.phen], markers.metric, " order")
            }
            else {
                phen.names2 <- phen.names
                main.string <- paste(test.file.prefix, " - ", 
                  phen.names2[k.phen], markers.metric, " order")
            }
        }
        else {
            phen.names2 <- phen.names[1]
            main.string <- paste(test.file.prefix, " - ", phen.names2, 
                markers.metric, " order")
        }
        x <- rbind(sample.names2, cls.list2, cls.labels2)
        print(paste("inside loop after phen sort before figure phen=", 
            k.phen))
        print(x)
        print(cls.phen2)
        print(phen.names2)
        phen.list <- unlist(cls.phen2)
        colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
        print(rbind(phen.list, colors.list))
        markers <- paste(markers, seq(1, length(markers)), sep = "_")
        MSIG.HeatMapPlot.7(V = V3, row.names = markers, row.names2 = model.descs2, 
            col.labels = cls.labels2, col.classes = cls.phen2, 
            phen.cmap = colors.list, phen.names = phen.names2, 
            col.names = sample.names2, main = main.string, xlab = "  ", 
            ylab = "  ", row.norm = row.norm, cmap.type = cmap.type, 
            char.rescale = char.res, legend = legend)
        V3 <- data.frame(V3)
        colnames(V3) <- sample.names2
        row.names(V3) <- markers
        if (!is.null(markers.file)) {
            write.gct(gct.data.frame = V3, descs = model.descs2, 
                filename = markers.file)
        }
    }
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
