\name{OPAM.create.model}
\alias{OPAM.create.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.create.model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.create.model(input.ds, input.cls, input2.ds, input2.cls, models.dir, target.class, target.class2, model.name, model.description, sample.norm.type = "rank", marker.disc = "MEAN.DIFF", top.markers.up = 20, top.markers.dn = 20, top.markers.up2 = 20, top.markers.dn2 = 20, statistic = "area.under.RES", weight = 0.25, msig.cntrl.genes = NA, random.seed = 12345, nperm = 0, link.function = "logit", burnin.iter = 5000, mcmc.iter = 25000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{input2.ds}{
%%     ~~Describe \code{input2.ds} here~~
}
  \item{input2.cls}{
%%     ~~Describe \code{input2.cls} here~~
}
  \item{models.dir}{
%%     ~~Describe \code{models.dir} here~~
}
  \item{target.class}{
%%     ~~Describe \code{target.class} here~~
}
  \item{target.class2}{
%%     ~~Describe \code{target.class2} here~~
}
  \item{model.name}{
%%     ~~Describe \code{model.name} here~~
}
  \item{model.description}{
%%     ~~Describe \code{model.description} here~~
}
  \item{sample.norm.type}{
%%     ~~Describe \code{sample.norm.type} here~~
}
  \item{marker.disc}{
%%     ~~Describe \code{marker.disc} here~~
}
  \item{top.markers.up}{
%%     ~~Describe \code{top.markers.up} here~~
}
  \item{top.markers.dn}{
%%     ~~Describe \code{top.markers.dn} here~~
}
  \item{top.markers.up2}{
%%     ~~Describe \code{top.markers.up2} here~~
}
  \item{top.markers.dn2}{
%%     ~~Describe \code{top.markers.dn2} here~~
}
  \item{statistic}{
%%     ~~Describe \code{statistic} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{msig.cntrl.genes}{
%%     ~~Describe \code{msig.cntrl.genes} here~~
}
  \item{random.seed}{
%%     ~~Describe \code{random.seed} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{link.function}{
%%     ~~Describe \code{link.function} here~~
}
  \item{burnin.iter}{
%%     ~~Describe \code{burnin.iter} here~~
}
  \item{mcmc.iter}{
%%     ~~Describe \code{mcmc.iter} here~~
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
function (input.ds, input.cls, input2.ds, input2.cls, models.dir, 
    target.class, target.class2, model.name, model.description, 
    sample.norm.type = "rank", marker.disc = "MEAN.DIFF", top.markers.up = 20, 
    top.markers.dn = 20, top.markers.up2 = 20, top.markers.dn2 = 20, 
    statistic = "area.under.RES", weight = 0.25, msig.cntrl.genes = NA, 
    random.seed = 12345, nperm = 0, link.function = "logit", 
    burnin.iter = 5000, mcmc.iter = 25000) 
{
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    library(gtools)
    library(verification)
    library(ROCR)
    library(MASS)
    library(RColorBrewer)
    library(MCMCpack)
    c1 <- c("black", "lightgrey")
    set.seed(random.seed)
    models.file <- paste(models.dir, "/", model.name, ".mod", 
        sep = "")
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gene.names <- dataset$row.names
    gene.descs <- dataset$descs
    sample.names <- dataset$names
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    CLS <- MSIG.ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
    if (is.na(match(target.class, class.phen))) 
        stop(c("target class is not phenotype in:", input.cls))
    for (i in 1:length(class.list)) class.labels[i] <- ifelse(class.list[i] == 
        target.class, 1, 2)
    col.index <- order(class.labels, decreasing = F)
    for (j in 1:Ng) m[j, ] <- m[j, col.index]
    sample.names <- sample.names[col.index]
    class.labels <- class.labels[col.index]
    class.list <- class.list[col.index]
    class.phen <- c(target.class, "CNTL")
    control.class <- "CNTL"
    if (sample.norm.type == "rank") {
        for (j in 1:Ns) {
            m[, j] <- rank(m[, j], ties.method = "average")
        }
        m <- 10000 * m/Ng
    }
    else if (sample.norm.type == "log.rank") {
        for (j in 1:Ns) {
            m[, j] <- rank(m[, j], ties.method = "average")
        }
        m <- log(10000 * m/Ng + exp(1))
    }
    else if (sample.norm.type == "log") {
        m[m < 1] <- 1
        m <- log(m + exp(1))
    }
    if (!is.na(msig.cntrl.genes)) {
        gene.names.int <- intersect(msig.cntrl.genes, gene.names)
        locs <- match(gene.names.int, gene.names, nomatch = 0)
        msig.cntrl <- m[locs, ]
        msig.cntrl.genes <- gene.names[locs]
        msig.cntrl.descs <- gene.descs[locs]
        msig.cntrl.size <- length(locs)
        if (msig.cntrl.size < 1) 
            msig.cntrl.center <- rep(1, Ns)
        else if (msig.cntrl.size == 1) 
            msig.cntrl.center <- msig.cntrl
        else if (msig.cntrl.size > 1) 
            msig.cntrl.center <- apply(msig.cntrl, MARGIN = 2, 
                FUN = mean)
        for (i in 1:Ng) {
            m[i, ] <- m[i, ]/msig.cntrl.center
        }
    }
    temp <- Gene.ranking(m, class.labels, method = marker.disc)
    gene.index <- order(temp, decreasing = T)
    gene.scores <- temp[gene.index]
    msig.up <- m[gene.index[1:top.markers.up], ]
    msig.up.size <- top.markers.up
    msig.up.genes <- gene.names[gene.index[1:top.markers.up]]
    msig.up.descs <- gene.descs[gene.index[1:top.markers.up]]
    msig.dn <- m[gene.index[seq(Ng, Ng - top.markers.dn + 1, 
        -1)], ]
    msig.dn.size <- top.markers.dn
    msig.dn.genes <- gene.names[gene.index[seq(Ng, Ng - top.markers.dn + 
        1, -1)]]
    msig.dn.descs <- gene.descs[gene.index[seq(Ng, Ng - top.markers.dn + 
        1, -1)]]
    print("Signatures (up/dn) created from gene marker selection")
    print(c("msig.up.size:", msig.up.size))
    print(c("msig.up.genes:", msig.up.genes))
    print(c("..."))
    print(c("msig.dn.size:", msig.dn.size))
    print(c("msig.dn.genes:", msig.dn.genes))
    print(c("..."))
    if (!is.na(msig.cntrl.genes)) {
        print(c("msig.cntrl.size:", msig.cntrl.size))
        print(c("msig.cntrl.genes:", msig.cntrl.genes))
    }
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(V = msig.up, row.names = msig.up.genes, 
        col.labels = class.labels, col.classes = class.phen, 
        phen.cmap = c1, col.names = sample.names, main = paste(model.name, 
            " UP signature"), xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1)
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(V = msig.dn, row.names = msig.dn.genes, 
        col.labels = class.labels, col.classes = class.phen, 
        phen.cmap = c1, col.names = sample.names, main = paste(model.name, 
            " DN signature"), xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1)
    if (!is.na(msig.cntrl.genes)) {
        x11(height = 9, width = 12)
        MSIG.HeatMapPlot.3(V = msig.cntrl, row.names = msig.cntrl.genes, 
            col.labels = class.labels, col.classes = class.phen, 
            phen.cmap = c1, col.names = sample.names, main = paste(model.name, 
                " CNTRL signature"), xlab = " ", ylab = " ", 
            sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1)
    }
    dataset <- MSIG.Gct2Frame(filename = input2.ds)
    m2 <- data.matrix(dataset$ds)
    gene.names2 <- dataset$row.names
    gene.descs2 <- dataset$descs
    sample.names2 <- dataset$names
    Ns2 <- length(m2[1, ])
    Ng2 <- length(m2[, 1])
    temp <- strsplit(input2.ds, split = "/")
    s <- length(temp[[1]])
    temp <- temp[[1]][4]
    temp <- strsplit(temp, split = ".gct")
    results.prefix2 <- temp[[1]][1]
    results.prefix2
    CLS <- MSIG.ReadClsFile(file = input2.cls)
    class.labels2 <- CLS$class.v
    class.phen2 <- CLS$phen
    class.list2 <- CLS$class.list
    if (is.na(match(target.class2, class.phen2))) 
        stop(c("target class2 is not phenotype in:", input2.cls))
    for (i in 1:length(class.list2)) class.labels2[i] <- ifelse(class.list2[i] == 
        target.class2, 1, 2)
    col.index <- order(class.labels2, decreasing = F)
    for (j in 1:Ng2) m2[j, ] <- m2[j, col.index]
    sample.names2 <- sample.names2[col.index]
    class.labels2 <- class.labels2[col.index]
    class.list2 <- class.list2[col.index]
    class.phen2 <- c(target.class2, "CNTL")
    control.class2 <- "CNTL"
    if (sample.norm.type == "rank") {
        for (j in 1:Ns2) {
            m2[, j] <- rank(m2[, j], ties.method = "average")
        }
        m2 <- 10000 * m2/Ng2
    }
    else if (sample.norm.type == "log.rank") {
        for (j in 1:Ns2) {
            m2[, j] <- rank(m2[, j], ties.method = "average")
        }
        m2 <- log(10000 * m2/Ng2 + exp(1))
    }
    else if (sample.norm.type == "log") {
        m2[m2 < 1] <- 1
        m2 <- log(m2 + exp(1))
    }
    if (!is.na(msig.cntrl.genes)) {
        gene.names.int2 <- intersect(msig.cntrl.genes, gene.names2)
        locs <- match(gene.names.int2, gene.names2, nomatch = 0)
        msig.cntrl2 <- m2[locs, ]
        msig.cntrl.genes2 <- gene.names2[locs]
        msig.cntrl.descs2 <- gene.descs2[locs]
        msig.cntrl.size2 <- length(locs)
        if (msig.cntrl.size2 < 1) 
            msig.cntrl.center2 <- rep(1, Ns2)
        else if (msig.cntrl.size2 == 1) 
            msig.cntrl.center2 <- msig.cntrl2
        else if (msig.cntrl.size2 > 1) 
            msig.cntrl.center2 <- apply(msig.cntrl2, MARGIN = 2, 
                FUN = mean)
        for (i in 1:Ng2) {
            m2[i, ] <- m2[i, ]/msig.cntrl.center2
        }
    }
    gene.names.int2 <- intersect(msig.up.genes, gene.names2)
    locs <- match(gene.names.int2, gene.names2, nomatch = 0)
    msig.up2 <- m2[locs, ]
    msig.up.genes2 <- gene.names2[locs]
    msig.up.descs2 <- gene.descs2[locs]
    msig.up.size2 <- length(locs)
    gene.names.int2 <- intersect(msig.dn.genes, gene.names2)
    locs <- match(gene.names.int2, gene.names2, nomatch = 0)
    msig.dn2 <- m2[locs, ]
    msig.dn.genes2 <- gene.names2[locs]
    msig.dn.descs2 <- gene.descs2[locs]
    msig.dn.size2 <- length(locs)
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(V = msig.up2, row.names = msig.up.genes2, 
        col.labels = class.labels2, col.classes = class.phen2, 
        phen.cmap = c1, col.names = sample.names2, main = paste(model.name, 
            " UP signature2"), xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1)
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(V = msig.dn2, row.names = msig.dn.genes2, 
        col.labels = class.labels2, col.classes = class.phen2, 
        phen.cmap = c1, col.names = sample.names2, main = paste(model.name, 
            " DN signature2"), xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1)
    if (!is.na(msig.cntrl.genes)) {
        x11(height = 9, width = 12)
        MSIG.HeatMapPlot.3(V = msig.cntrl2, row.names = msig.cntrl.genes2, 
            col.labels = class.labels2, col.classes = class.phen2, 
            phen.cmap = c1, col.names = sample.names2, main = paste(model.name, 
                " CNTRL signature2"), xlab = " ", ylab = " ", 
            sub = " ", row.norm = T, cmap.type = 4, char.rescale = 1)
    }
    m3 <- rbind(msig.up2, msig.dn2)
    Ng3 <- length(m3[, 1])
    gene.names3 <- c(msig.up.genes2, msig.dn.genes2)
    gene.descs3 <- c(msig.up.descs2, msig.dn.descs2)
    temp <- Gene.ranking(m3, class.labels2, method = marker.disc)
    gene.index <- order(temp, decreasing = T)
    gene.scores3 <- temp[gene.index]
    msig.up3 <- m3[gene.index[1:top.markers.up2], ]
    msig.up.size3 <- top.markers.up2
    msig.up.genes3 <- gene.names3[gene.index[1:top.markers.up2]]
    msig.up.descs3 <- gene.descs3[gene.index[1:top.markers.up2]]
    msig.dn3 <- m3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 
        1, -1)], ]
    msig.dn.size3 <- top.markers.dn2
    msig.dn.genes3 <- gene.names3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 
        1, -1)]]
    msig.dn.descs3 <- gene.descs3[gene.index[seq(Ng3, Ng3 - top.markers.dn2 + 
        1, -1)]]
    print("Signatures3 (up/dn) created from gene marker selection")
    print(c("msig.up.size3:", msig.up.size3))
    print(c("msig.up.genes3:", msig.up.genes3))
    print(c("..."))
    print(c("msig.dn.size3:", msig.dn.size3))
    print(c("msig.dn.genes3:", msig.dn.genes3))
    print(c("..."))
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(V = msig.up3, row.names = msig.up.genes3, 
        col.labels = class.labels2, col.classes = class.phen2, 
        phen.cmap = c1, col.names = sample.names2, main = paste(model.name, 
            " UP signature3"), xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1)
    x11(height = 9, width = 12)
    MSIG.HeatMapPlot.3(V = msig.dn3, row.names = msig.dn.genes3, 
        col.labels = class.labels2, col.classes = class.phen2, 
        phen.cmap = c1, col.names = sample.names2, main = paste(model.name, 
            " DN signature3"), xlab = " ", ylab = " ", sub = " ", 
        row.norm = T, cmap.type = 4, char.rescale = 1)
    OPAM <- OPAM.Projection(m2, gene.names2, Ns2, Ng2, weight, 
        statistic, msig.up.genes3, nperm = nperm)
    score.up <- OPAM$ES.vector
    OPAM <- OPAM.Projection(m2, gene.names2, Ns2, Ng2, weight, 
        statistic, msig.dn.genes3, nperm = nperm)
    score.dn <- OPAM$ES.vector
    score <- score.up - score.dn
    x11(width = 14, height = 9)
    nf <- layout(matrix(c(1, 2, 3, 0, 4, 0), 3, 2, byrow = T), 
        widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
    par(mar = c(2, 4, 2, 4))
    barplot(score.up, main = "OPAM Score UP (refinement)", font.axis = 1.25, 
        cex.lab = 1.5, cex.axis = 1.25, cex.names = 1.25, width = 1, 
        space = 0, col = c1[class.labels])
    leg.txt <- class.phen
    p.vec <- rep(22, length(leg.txt))
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c1, col = "black", cex = 1.25, 
        pt.cex = 2.5)
    par(mar = c(2, 4, 2, 4))
    barplot(score.dn, main = "OPAM Score DOWN (refinement)", 
        font.axis = 1.25, cex.lab = 1.5, cex.axis = 1.25, cex.names = 1.25, 
        width = 1, space = 0, col = c1[class.labels])
    par(mar = c(2, 4, 2, 4))
    barplot(score, main = "OPAM Total Score (refinement)", font.axis = 1.25, 
        cex.lab = 1.5, cex.axis = 1.25, cex.names = 1.25, width = 1, 
        space = 0, col = c1[class.labels])
    target.var <- ifelse(class.list2 == target.class2, 1, 0)
    Bayesian.function <- ifelse(link.function == "logit", "MCMClogit(", 
        "MCMCprobit(")
    model.formula <- paste(Bayesian.function, "target.var ~ score,  burnin = burnin.iter, mcmc = mcmc.iter, bayes.resid=T)", 
        sep = "")
    model.formula
    reg.model <- eval(parse(text = model.formula))
    beta.0 <- reg.model[, 1]
    beta.1 <- reg.model[, 2]
    model.formula <- "beta.0 + beta.1 * score[i]"
    model.formula
    prob.i <- matrix(0, nrow = Ns2, ncol = 3)
    model.score <- vector(length = Ns2, mode = "numeric")
    for (i in 1:Ns2) {
        model.score[i] <- eval(parse(text = model.formula))
        if (link.function == "logit") {
            p.vec <- paste("inv.logit(x=", model.formula, ")", 
                sep = "")
        }
        else if (link.function == "probit") {
            p.vec <- paste("(erf(", model.formula, ") + 1)/2", 
                sep = "")
        }
        else {
            stop("Unknown link function")
        }
        val <- eval(parse(text = p.vec))
        prob.i[i, 1] <- quantile(val, probs = 0.5)
        prob.i[i, 2] <- quantile(val, probs = 0.05)
        prob.i[i, 3] <- quantile(val, probs = 0.95)
    }
    xmin <- min(model.score)
    xmax <- max(model.score)
    range.x <- xmax - xmin
    n.points <- 1000
    prob.m <- matrix(0, nrow = n.points, ncol = 3)
    x.m <- vector(length = n.points, mode = "numeric")
    for (k in 1:n.points) {
        x.m[k] <- xmin + k * (range.x/n.points)
        if (link.function == "logit") {
            p.vec <- paste("inv.logit(x=", x.m[k], ")", sep = "")
        }
        else if (link.function == "probit") {
            p.vec <- paste("(erf(", x.m[k], ") + 1)/2", sep = "")
        }
        else {
            stop("Unknown link function")
        }
        val <- eval(parse(text = p.vec))
        prob.m[k, 1] <- quantile(val, probs = 0.5)
        prob.m[k, 2] <- quantile(val, probs = 0.05)
        prob.m[k, 3] <- quantile(val, probs = 0.95)
    }
    istar <- which.min(abs(0.5 - prob.m[, 1]))
    istar <- xmin + istar * (range.x/1000)
    x.index <- order(model.score, decreasing = F)
    x.order <- model.score[x.index]
    prob.i.order <- prob.i[x.index, ]
    target.var.order <- ifelse(target.var[x.index] == 1, c1[1], 
        c1[2])
    class.labels.order <- class.labels[x.index]
    boundary <- istar
    pred.class <- ifelse(prob.i.order[, 1] >= 0.5, 2, 1)
    z.range <- range(x.order)
    x11(height = 7, width = 9.5)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 
        1), heights = 1, respect = FALSE)
    plot(x.order, prob.i.order[, 1], sub = model.name, pch = 20, 
        col = 0, cex = 2, xlab = "Activation Index", ylab = "Probability")
    points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, lty = 1, 
        cex = 1)
    points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, cex = 1)
    points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, cex = 1)
    arrows(x.order, prob.i.order[, 2], x.order, prob.i.order[, 
        3], col = 4, angle = 90, code = 3, length = 0)
    range.x <- range(x.order)
    points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, 
        lwd = 2)
    points(range.x, c(-0.15, -0.15), type = "l", lty = 1, col = 1, 
        lwd = 2)
    points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, 
        col = 1, lwd = 2)
    points(x.order, prob.i.order[, 1], pch = 21, bg = target.var.order, 
        col = 1, cex = 2)
    points(x.order, rep(-0.15, length(x.order)), pch = 21, bg = target.var.order, 
        col = 1, cex = 2)
    leg.txt <- class.phen
    p.vec <- rep(22, length(leg.txt))
    c.vec <- c1
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.25, pt.cex = 2.5)
    x11(width = 14, height = 9)
    nf <- layout(matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T), 
        widths = c(5, 1), heights = c(1, 1, 1), respect = FALSE)
    par(mar = c(4, 7, 4, 7))
    MSIG.Score.Plot(z = score, main = paste(model.name, " Model Score (train)"), 
        phen.cmap = c1, char.rescale = 1, col.classes = class.phen, 
        col.labels = class.labels, xlab = "Samples", ylab = "Score", 
        create.window = F, create.legend = T)
    par(mar = c(4, 7, 4, 7))
    norm.score <- (score - min(score))/(max(score) - min(score))
    MSIG.Score.Plot(z = norm.score, main = paste(model.name, 
        " Normalized Model Score (train)"), phen.cmap = c1, char.rescale = 1, 
        col.classes = class.phen, col.labels = class.labels, 
        xlab = "Samples", ylab = "Normalized Score", create.window = F, 
        create.legend = T)
    par(mar = c(4, 7, 4, 7))
    MSIG.Score.Plot(z = prob.i[, 1], main = paste(model.name, 
        " Probabiliy (train)"), phen.cmap = c1, char.rescale = 1, 
        col.classes = class.phen, col.labels = class.labels, 
        xlab = "Samples", ylab = "Probability", create.window = F, 
        create.legend = T)
    MSIG.HeatMapPlot.6(V = rbind(score, model.score, t(prob.i[, 
        1])), row.names = c("raw.score", "model.score", "probability"), 
        row.names2 = c(model.name, model.name, model.name), col.labels = class.labels2, 
        col.labels2 = class.labels2, col.classes = class.phen2, 
        phen.cmap = c1, phen.names = model.name, col.names = sample.names2, 
        main = model.name, xlab = "  ", ylab = "  ", sub = "   ", 
        row.norm = T, cmap.type = 3, char.rescale = 1, legend = T)
    model.creation.date <- date()
    OPAM.write.param.line(param = model.creation.date, param.name = "model.creation.date", 
        file = models.file, append = F)
    OPAM.write.param.line(param = input.ds, param.name = "input.ds", 
        file = models.file, append = T)
    OPAM.write.param.line(param = input.ds, param.name = "input.ds", 
        file = models.file, append = T)
    OPAM.write.param.line(param = input.cls, param.name = "input.cls", 
        file = models.file, append = T)
    OPAM.write.param.line(param = input2.ds, param.name = "input2.ds", 
        file = models.file, append = T)
    OPAM.write.param.line(param = input2.cls, param.name = "input2.cls", 
        file = models.file, append = T)
    OPAM.write.param.line(param = target.class, param.name = "target.class", 
        file = models.file, append = T)
    OPAM.write.param.line(param = target.class2, param.name = "target.class2", 
        file = models.file, append = T)
    OPAM.write.param.line(param = model.name, param.name = "model.name", 
        file = models.file, append = T)
    OPAM.write.param.line(param = model.description, param.name = "model.description", 
        file = models.file, append = T)
    OPAM.write.param.line(param = sample.norm.type, param.name = "sample.norm.type", 
        file = models.file, append = T)
    OPAM.write.param.line(param = marker.disc, param.name = "marker.disc", 
        file = models.file, append = T)
    OPAM.write.param.line(param = top.markers.up, param.name = "top.markers.up", 
        file = models.file, append = T)
    OPAM.write.param.line(param = top.markers.dn, param.name = "top.markers.dn", 
        file = models.file, append = T)
    OPAM.write.param.line(param = top.markers.up2, param.name = "top.markers.up2", 
        file = models.file, append = T)
    OPAM.write.param.line(param = top.markers.dn2, param.name = "top.markers.dn2", 
        file = models.file, append = T)
    OPAM.write.param.line(param = statistic, param.name = "statistic", 
        file = models.file, append = T)
    OPAM.write.param.line(param = weight, param.name = "weight", 
        file = models.file, append = T)
    OPAM.write.param.line(param = random.seed, param.name = "random.seed", 
        file = models.file, append = T)
    OPAM.write.param.line(param = nperm, param.name = "nperm", 
        file = models.file, append = T)
    OPAM.write.param.line(param = link.function, param.name = "link.function", 
        file = models.file, append = T)
    OPAM.write.param.line(param = c1, param.name = "c1", file = models.file, 
        append = T)
    OPAM.write.param.line(param = msig.cntrl.genes, param.name = "msig.cntrl.genes", 
        file = models.file, append = T)
    OPAM.write.param.line(param = msig.up.genes, param.name = "msig.up.genes", 
        file = models.file, append = T)
    OPAM.write.param.line(param = msig.dn.genes, param.name = "msig.dn.genes", 
        file = models.file, append = T)
    OPAM.write.param.line(param = msig.up.genes2, param.name = "msig.up.genes2", 
        file = models.file, append = T)
    OPAM.write.param.line(param = msig.dn.genes2, param.name = "msig.dn.genes2", 
        file = models.file, append = T)
    OPAM.write.param.line(param = msig.up.genes3, param.name = "msig.up.genes3", 
        file = models.file, append = T)
    OPAM.write.param.line(param = msig.dn.genes3, param.name = "msig.dn.genes3", 
        file = models.file, append = T)
    OPAM.write.param.line(param = beta.0, param.name = "beta.0", 
        file = models.file, append = T)
    OPAM.write.param.line(param = beta.1, param.name = "beta.1", 
        file = models.file, append = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
