\name{OPAM.analyze.projection}
\alias{OPAM.analyze.projection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.analyze.projection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.analyze.projection(input.ds, input.cls, results.dir, top.class, top.phen, normalize.score = T, normalization.type = zero.one, feature.sel.thres = 0.05, markers.num = 5, user.colors = NA, k.proj = NA, markers.metric = "S2N", markers.file = NULL)
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
  \item{top.class}{
%%     ~~Describe \code{top.class} here~~
}
  \item{top.phen}{
%%     ~~Describe \code{top.phen} here~~
}
  \item{normalize.score}{
%%     ~~Describe \code{normalize.score} here~~
}
  \item{normalization.type}{
%%     ~~Describe \code{normalization.type} here~~
}
  \item{feature.sel.thres}{
%%     ~~Describe \code{feature.sel.thres} here~~
}
  \item{markers.num}{
%%     ~~Describe \code{markers.num} here~~
}
  \item{user.colors}{
%%     ~~Describe \code{user.colors} here~~
}
  \item{k.proj}{
%%     ~~Describe \code{k.proj} here~~
}
  \item{markers.metric}{
%%     ~~Describe \code{markers.metric} here~~
}
  \item{markers.file}{
%%     ~~Describe \code{markers.file} here~~
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
function (input.ds, input.cls, results.dir, top.class, top.phen, 
    normalize.score = T, normalization.type = zero.one, feature.sel.thres = 0.05, 
    markers.num = 5, user.colors = NA, k.proj = NA, markers.metric = "S2N", 
    markers.file = NULL) 
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
    sample.names <- dataset$names
    Ns <- length(m[1, ])
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
        cls.labels2 <- cls.labels
        n.phen <- 1
    }
    else {
        cls.labels2 <- as.vector(cls.labels[top.phen, ])
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
                name = "Dark2"), brewer.pal(n = 7, name = "Set1"))
        }
    }
    if (!is.null(CLS$phen.names)) {
        phen.names <- CLS$phen.names
    }
    else {
        phen.names = "NA"
    }
    print("CLS:")
    print(CLS)
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
    print("computing ROC...")
    roc.array <- array(data = 0, dim = c(n.models, n.phen, max.classes), 
        dimnames = NULL)
    p.val.array <- array(data = 0, dim = c(n.models, n.phen, 
        max.classes), dimnames = NULL)
    for (i in 1:n.models) {
        for (j in 1:n.phen) {
            for (k in 1:n.classes[[j]]) {
                if (is.vector(cls.labels)) {
                  bin.class <- ifelse(cls.labels == k, 1, 0)
                }
                else {
                  bin.class <- ifelse(cls.labels[j, ] == k, 1, 
                    0)
                }
                m.score <- m[i, ]
                m.score.norm <- (m.score - min(m.score))/(max(m.score) - 
                  min(m.score))
                perf.auc <- roc.area(bin.class, m.score.norm)
                roc.array[i, j, k] <- signif(perf.auc$A, digits = 3)
                p.val.array[i, j, k] <- signif(perf.auc$p.value, 
                  digits = 3)
            }
        }
    }
    print("computing ROC for top class")
    top.roc.vector <- vector(length = n.models, mode = "numeric")
    top.p.val.vector <- vector(length = n.models, mode = "numeric")
    for (i in 1:n.models) {
        if (is.vector(cls.labels)) {
            bin.class <- ifelse(cls.list == top.class, 1, 0)
        }
        else {
            bin.class <- ifelse(cls.list[top.phen, ] == top.class, 
                1, 0)
        }
        if (markers.metric == "ROC") {
            m.score <- m[i, ]
            m.score.norm <- (m.score - min(m.score))/(max(m.score) - 
                min(m.score))
            perf.auc <- roc.area(bin.class, m.score.norm)
            top.roc.vector[i] <- signif(perf.auc$A, digits = 3)
            top.p.val.vector[i] <- signif(perf.auc$p.value, digits = 3)
        }
        else if (markers.metric == "S2N") {
            top.roc.vector[i] <- signif(S2N(m[i, ], bin.class), 
                digits = 3)
            temp <- split(m[i, ], bin.class)
            x <- temp[[1]]
            y <- temp[[2]]
            t.test.out <- t.test(x = x, y = y)
            top.p.val.vector[i] <- signif(t.test.out$p.value, 
                digits = 3)
        }
    }
    roc.table <- data.frame(cbind(model.names, model.descs, top.roc.vector, 
        top.p.val.vector))
    roc.order <- order(top.roc.vector, decreasing = T)
    roc.table <- roc.table[roc.order, ]
    names(roc.table) <- c("Feature:", "Description:", markers.metric, 
        "p-value:")
    top.pos.features.m <- m[top.p.val.vector <= feature.sel.thres, 
        ]
    top.pos.features <- model.names[top.p.val.vector <= feature.sel.thres]
    top.pos.features.descs <- model.descs[top.p.val.vector <= 
        feature.sel.thres]
    top.pos.features.roc <- top.roc.vector[top.p.val.vector <= 
        feature.sel.thres]
    top.pos.features.p.val <- top.p.val.vector[top.p.val.vector <= 
        feature.sel.thres]
    top.neg.features.m <- m[top.p.val.vector >= 1 - feature.sel.thres, 
        ]
    top.neg.features <- model.names[top.p.val.vector >= 1 - feature.sel.thres]
    top.neg.features.descs <- model.descs[top.p.val.vector >= 
        1 - feature.sel.thres]
    top.neg.features.roc <- top.roc.vector[top.p.val.vector >= 
        1 - feature.sel.thres]
    top.neg.features.p.val <- top.p.val.vector[top.p.val.vector >= 
        1 - feature.sel.thres]
    top.features <- data.frame(rbind(cbind(top.pos.features, 
        top.pos.features.descs, top.pos.features.roc, top.pos.features.p.val), 
        cbind(top.neg.features, top.neg.features.descs, top.neg.features.roc, 
            top.neg.features.p.val)))
    names(top.features) <- c("Feature:", "Description:", markers.metric, 
        "p-value:")
    top.m <- rbind(top.pos.features.m, top.neg.features.m)
    top.features.order <- order(top.features[, 3], decreasing = T)
    top.features <- top.features[top.features.order, ]
    top.features
    print("plot in original order")
    height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 
        0.44 + 5)
    x11(width = 14, height = height)
    MSIG.HeatMapPlot.6(V = m, row.names = model.names, row.names2 = model.descs, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names, 
        main = paste(test.file.prefix, "-- Original Order"), 
        xlab = "  ", ylab = "  ", sub = "Original order", row.norm = T, 
        cmap.type = 3, char.rescale = char.res, legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    MSIG.HeatMapPlot.6(V = m, row.names = model.names, row.names2 = model.descs, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names, 
        main = paste(test.file.prefix, "-- Original Order"), 
        xlab = "  ", ylab = "  ", sub = "Original order", row.norm = T, 
        cmap.type = 3, char.rescale = char.res, legend = T)
    dev.off()
    cls.labels.renam <- cls.labels
    if (!is.vector(cls.labels)) {
        max.lab <- max(cls.labels[, 1])
        for (k in 2:length(cls.labels[, 1])) {
            cls.labels.renam[k, ] <- cls.labels[k, ] + max.lab
            max.label <- max(cls.labels[k, ])
        }
    }
    row.order <- order(top.roc.vector, decreasing = TRUE)
    V1 <- m[row.order, ]
    model.names2 <- model.names[row.order]
    model.descs2 <- paste(model.descs[row.order], top.roc.vector[row.order], 
        top.p.val.vector[row.order])
    score.index <- order(V1[1, ], decreasing = TRUE)
    cls.labels2.sorted <- cls.labels2[score.index]
    if (is.vector(cls.labels)) {
        cls.labels.sorted <- cls.labels[score.index]
    }
    else {
        cls.labels.sorted <- cls.labels[, score.index]
    }
    sample.names.sorted <- sample.names[score.index]
    print("plot in top model order. Best match to target class")
    height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 
        0.44 + 5)
    x11(width = 14, height = height)
    MSIG.HeatMapPlot.6(V = V1[, score.index], row.names = model.names2, 
        row.names2 = model.descs2, col.labels = cls.labels.sorted, 
        col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names, 
        col.names = sample.names.sorted, main = paste(test.file.prefix, 
            "-- Model Order", markers.metric), xlab = "  ", ylab = "  ", 
        sub = "Top Model order", row.norm = T, cmap.type = 3, 
        char.rescale = char.res, legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP.MODEL", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    MSIG.HeatMapPlot.6(V = V1[, score.index], row.names = model.names2, 
        row.names2 = model.descs2, col.labels = cls.labels.sorted, 
        col.classes = cls.phen, phen.cmap = c.test, phen.names = phen.names, 
        col.names = sample.names.sorted, main = paste(test.file.prefix, 
            "-- Top Model Order"), xlab = "  ", ylab = "  ", 
        sub = "Model order", row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    dev.off()
    dist.matrix <- dist(m)
    HC <- hclust(dist.matrix, method = "complete")
    V2 <- m[HC$order, ]
    model.names2 <- model.names[HC$order]
    model.descs2 <- model.descs[HC$order]
    top.roc.vector2 <- top.roc.vector[HC$order]
    top.p.val.vector2 <- top.p.val.vector[HC$order]
    sample.names2 <- sample.names
    num.phen <- length(unique(cls.labels2))
    for (k in 1:num.phen) {
        V3 <- V2[, cls.labels2 == k]
        sn <- sample.names2[cls.labels2 == k]
        if (is.vector(cls.labels)) {
            cl <- cls.labels[cls.labels2 == k]
        }
        else {
            cl <- cls.labels[, cls.labels2 == k]
        }
        dist.matrix <- dist(t(V3))
        HC <- hclust(dist.matrix, method = "complete")
        V3 <- V3[, HC$order]
        sn <- sn[HC$order]
        if (is.vector(cls.labels)) {
            cl <- cl[HC$order]
            cls.labels[cls.labels2 == k] <- cl
        }
        else {
            cl <- cl[, HC$order]
            cls.labels[, cls.labels2 == k] <- cl
        }
        V2[, cls.labels2 == k] <- V3
        sample.names2[cls.labels2 == k] <- sn
    }
    print("plot sorted inside each phenotype")
    model.descs2 <- paste(model.descs2, top.roc.vector2, top.p.val.vector2)
    height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 
        0.44 + 5)
    x11(width = 14, height = height)
    MSIG.HeatMapPlot.6(V = V2, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names2, 
        main = paste(test.file.prefix, "-- Sorted Inside Class"), 
        xlab = "  ", ylab = "  ", sub = "Sorted Inside Class", 
        row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP.SORT.PHEN", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    MSIG.HeatMapPlot.6(V = V2, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names2, 
        main = paste(test.file.prefix, "-- Sorted Inside Class"), 
        xlab = "  ", ylab = "  ", sub = "Sorted Inside Class", 
        row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    dev.off()
    print("plot in marker order")
    if (is.vector(cls.labels)) {
        classes <- unique(cls.list)
    }
    else {
        classes <- unique(cls.list[top.phen, ])
    }
    if (length(classes) > 2) {
        only.up <- T
    }
    else {
        only.up <- F
    }
    if (markers.num == 0) {
        if (is.vector(cls.labels)) {
            bin.class <- ifelse(cls.list == top.class, 0, 1)
        }
        else {
            bin.class <- ifelse(cls.list[top.phen, ] == top.class, 
                0, 1)
        }
        if (markers.metric == "S2N") {
            metric <- Gene.ranking(m, bin.class, method = "RS2N")
        }
        else if (markers.metric == "ROC") {
            bin.class <- ifelse(bin.class == 1, 0, 1)
            metric <- vector(length = n.models, mode = "numeric")
            for (i in 1:n.models) {
                m.score <- m[i, ]
                m.score.norm <- (m.score - min(m.score))/(max(m.score) - 
                  min(m.score))
                perf.auc <- roc.area(bin.class, m.score.norm)
                metric[i] <- perf.auc$A
            }
        }
        metric <- signif(metric, digits = 3)
        metric.order <- order(metric, decreasing = T)
        markers <- model.names[metric.order]
        markers.descs <- model.descs[metric.order]
        metric.list <- metric[metric.order]
        markers.num <- length(m[, 1])/2
        k.class <- rep(top.class, markers.num)
    }
    else {
        if (length(classes) == 2) 
            classes <- classes[1]
        markers <- NULL
        markers.descs <- NULL
        metric.list <- NULL
        k.class <- NULL
        for (k in classes) {
            if (is.vector(cls.labels)) {
                bin.class <- ifelse(cls.list == k, 0, 1)
            }
            else {
                bin.class <- ifelse(cls.list[top.phen, ] == k, 
                  0, 1)
            }
            if (markers.metric == "S2N") {
                metric <- Gene.ranking(m, bin.class, method = "RS2N")
            }
            else if (markers.metric == "ROC") {
                bin.class <- ifelse(bin.class == 1, 0, 1)
                metric <- vector(length = n.models, mode = "numeric")
                for (i in 1:n.models) {
                  m.score <- m[i, ]
                  m.score.norm <- (m.score - min(m.score))/(max(m.score) - 
                    min(m.score))
                  perf.auc <- roc.area(bin.class, m.score.norm)
                  metric[i] <- perf.auc$A
                }
            }
            metric <- signif(metric, digits = 2)
            metric.order <- order(metric, decreasing = T)
            if (only.up == T) {
                markers <- c(markers, model.names[metric.order][1:markers.num])
                markers.descs <- c(markers.descs, model.descs[metric.order][1:markers.num])
                metric.list <- c(metric.list, metric[metric.order][1:markers.num])
                k.class <- c(k.class, rep(k, markers.num))
            }
            else {
                markers <- c(markers, model.names[metric.order][1:markers.num], 
                  model.names[metric.order][(length(model.names) - 
                    markers.num + 1):length(model.names)])
                markers.descs <- c(markers.descs, model.descs[metric.order][1:markers.num], 
                  model.descs[metric.order][(length(model.names) - 
                    markers.num + 1):length(model.names)])
                metric.list <- c(metric.list, metric[metric.order][1:markers.num], 
                  metric[metric.order][(length(model.names) - 
                    markers.num + 1):length(model.names)])
                k.class <- c(k.class, rep(k, markers.num), rep(paste("not", 
                  k), markers.num))
            }
        }
    }
    print("markers")
    V3 <- V2[markers, ]
    model.descs2 <- paste(markers.descs, metric.list, k.class)
    height <- ifelse(length(markers) + n.phen >= 9, 9, (length(markers) + 
        n.phen) * 0.44 + 5)
    char.res <- 0.0125 * length(markers) + 0.7
    x11(width = 15, height = height)
    MSIG.HeatMapPlot.6(V = V3, row.names = markers, row.names2 = model.descs2, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names2, 
        main = paste(test.file.prefix, "-- Markers Order"), xlab = "  ", 
        ylab = "  ", sub = paste("Markers Order (metric: ", metric, 
            ")"), row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP.MARKERS.", 
        markers.metric, sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 14, 
        width = 11)
    MSIG.HeatMapPlot.6(V = V3, row.names = markers, row.names2 = model.descs2, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names2, 
        main = paste(test.file.prefix, "-- Markers Order"), xlab = "  ", 
        ylab = "  ", sub = paste("Markers Order (metric: ", metric, 
            ")"), row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    dev.off()
    names(V3)[seq(1, Ns)] <- sample.names2
    row.names(V3) <- paste(markers, seq(1, length(markers)), 
        sep = "_")
    if (!is.null(markers.file)) {
        write.gct(gct.data.frame = V3, descs = model.descs2, 
            filename = markers.file)
    }
    print("NMF projection")
    if (is.na(k.proj)) {
        k.proj <- ifelse(length(model.names2) < 10, 3, ceiling(length(model.names2)/8))
    }
    NMF.out <- NMF.div(V = V2, k = k.proj, maxniter = 4000, seed = 1234, 
        stopconv = 40, stopfreq = 10)
    H <- NMF.out$H
    W <- NMF.out$W
    k.proj.names <- paste("NMF_", seq(1, k.proj), sep = "")
    height <- ifelse(length(k.proj.names) + n.phen >= 9, 9, (length(markers) + 
        n.phen) * 0.44 + 5)
    char.res <- 0.0125 * length(markers) + 0.7
    print("H matrix")
    x11(width = 15, height = height)
    MSIG.HeatMapPlot.6(V = H, row.names = k.proj.names, row.names2 = k.proj.names, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names2, 
        main = paste(test.file.prefix, "-- H matrix"), xlab = "  ", 
        ylab = "  ", sub = "H Matrix", row.norm = T, cmap.type = 3, 
        char.rescale = char.res, legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP.SORT.PHEN.H.MATRIX", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    MSIG.HeatMapPlot.6(V = H, row.names = k.proj.names, row.names2 = k.proj.names, 
        col.labels = cls.labels, col.labels2 = cls.labels2, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names2, 
        main = paste(test.file.prefix, "-- H matrix"), xlab = "  ", 
        ylab = "  ", sub = "H Matrix", row.norm = T, cmap.type = 3, 
        char.rescale = char.res, legend = T)
    dev.off()
    print("W matrix")
    height <- ifelse(length(k.proj.names) + n.phen >= 9, 9, (length(markers) + 
        n.phen) * 0.44 + 5)
    char.res <- 0.0125 * length(markers) + 0.7
    x11(width = 15, height = height)
    MSIG.HeatMapPlot.6(V = W, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = seq(1, k.proj), col.classes = k.proj.names, 
        phen.cmap = c.test, col.names = k.proj.names, main = paste(test.file.prefix, 
            "-- W matrix"), xlab = "  ", ylab = "  ", sub = "W Matrix", 
        row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP.SORT.PHEN.W.MATRIX", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    MSIG.HeatMapPlot.6(V = W, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = seq(1, k.proj), col.classes = k.proj.names, 
        phen.cmap = c.test, col.names = k.proj.names, main = paste(test.file.prefix, 
            "-- W matrix"), xlab = "  ", ylab = "  ", sub = "W Matrix", 
        row.norm = T, cmap.type = 3, char.rescale = char.res, 
        legend = T)
    dev.off()
    print("sorting entire set")
    dist.matrix <- dist(m)
    HC <- hclust(dist.matrix, method = "complete")
    V2 <- m[HC$order, ]
    model.names2 <- model.names[HC$order]
    model.descs2 <- model.descs[HC$order]
    sample.names2 <- sample.names
    dist.matrix <- dist(t(V2))
    HC <- hclust(dist.matrix, method = "complete")
    V2 <- V2[, HC$order]
    sample.names.sorted <- sample.names2[HC$order]
    if (is.vector(cls.labels)) {
        cls.labels.sorted <- cls.labels[HC$order]
    }
    else {
        cls.labels.sorted <- cls.labels[, HC$order]
    }
    print("plot bisorted")
    height <- ifelse(n.models + n.phen >= 9, 9, (n.models + n.phen) * 
        0.44 + 5)
    x11(width = 14, height = height)
    MSIG.HeatMapPlot.6(V = V2, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = cls.labels.sorted, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names.sorted, 
        main = paste(test.file.prefix, "-- Bisorted"), xlab = "  ", 
        ylab = "  ", sub = "Bisorted", row.norm = F, cmap.type = 3, 
        char.rescale = char.res, legend = T)
    filename <- paste(results.dir, test.file.prefix, ".HEATMAP.BISORT", 
        sep = "")
    pdf(file = paste(filename, ".pdf", sep = ""), height = 8.5, 
        width = 11)
    MSIG.HeatMapPlot.6(V = V2, row.names = model.names2, row.names2 = model.descs2, 
        col.labels = cls.labels.sorted, col.classes = cls.phen, 
        phen.cmap = c.test, phen.names = phen.names, col.names = sample.names.sorted, 
        main = paste(test.file.prefix, "-- Bisorted"), xlab = "  ", 
        ylab = "  ", sub = "Bisorted", row.norm = F, cmap.type = 3, 
        char.rescale = char.res, legend = T)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
