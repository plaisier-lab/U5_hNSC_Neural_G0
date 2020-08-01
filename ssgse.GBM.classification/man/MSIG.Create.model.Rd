\name{MSIG.Create.model}
\alias{MSIG.Create.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Create.model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Create.model(gct.file, cls.file, output.dir, database.dir, identifier, target.class, column.subset = "ALL", column.sel.type = "samples", thres = "NULL", ceil = "NULL", shift = "NULL", fold = 1, delta = 0, norm = 6, seed = 1234, topgs = 25, link.function = "logit", model.type = "Classic", burnin.iter = 5000, mcmc.iter = 25000, col.target = "darkgreen", col.control = "yellow", no.call.range.min = 0.3, no.call.range.max = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gct.file}{
%%     ~~Describe \code{gct.file} here~~
}
  \item{cls.file}{
%%     ~~Describe \code{cls.file} here~~
}
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
}
  \item{database.dir}{
%%     ~~Describe \code{database.dir} here~~
}
  \item{identifier}{
%%     ~~Describe \code{identifier} here~~
}
  \item{target.class}{
%%     ~~Describe \code{target.class} here~~
}
  \item{column.subset}{
%%     ~~Describe \code{column.subset} here~~
}
  \item{column.sel.type}{
%%     ~~Describe \code{column.sel.type} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
}
  \item{ceil}{
%%     ~~Describe \code{ceil} here~~
}
  \item{shift}{
%%     ~~Describe \code{shift} here~~
}
  \item{fold}{
%%     ~~Describe \code{fold} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{norm}{
%%     ~~Describe \code{norm} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{topgs}{
%%     ~~Describe \code{topgs} here~~
}
  \item{link.function}{
%%     ~~Describe \code{link.function} here~~
}
  \item{model.type}{
%%     ~~Describe \code{model.type} here~~
}
  \item{burnin.iter}{
%%     ~~Describe \code{burnin.iter} here~~
}
  \item{mcmc.iter}{
%%     ~~Describe \code{mcmc.iter} here~~
}
  \item{col.target}{
%%     ~~Describe \code{col.target} here~~
}
  \item{col.control}{
%%     ~~Describe \code{col.control} here~~
}
  \item{no.call.range.min}{
%%     ~~Describe \code{no.call.range.min} here~~
}
  \item{no.call.range.max}{
%%     ~~Describe \code{no.call.range.max} here~~
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
function (gct.file, cls.file, output.dir, database.dir, identifier, 
    target.class, column.subset = "ALL", column.sel.type = "samples", 
    thres = "NULL", ceil = "NULL", shift = "NULL", fold = 1, 
    delta = 0, norm = 6, seed = 1234, topgs = 25, link.function = "logit", 
    model.type = "Classic", burnin.iter = 5000, mcmc.iter = 25000, 
    col.target = "darkgreen", col.control = "yellow", no.call.range.min = 0.3, 
    no.call.range.max = 0.7) 
{
    c1 <- c(col.target, col.control)
    set.seed(seed)
    O <- MSIG.Subset.Dataset(input.ds = gct.file, input.cls = cls.file, 
        column.subset = column.subset, column.sel.type = column.sel.type, 
        row.subset = "ALL", output.ds = paste(output.dir, "temp1.gct", 
            sep = ""), output.cls = paste(output.dir, "temp1.cls", 
            sep = ""))
    O <- MSIG.Preprocess.Dataset(input.ds = paste(output.dir, 
        "temp1.gct", sep = ""), output.ds = paste(output.dir, 
        "temp2.gct", sep = ""), thres = thres, ceil = ceil, fold = fold, 
        delta = delta, normalization = norm)
    dataset <- MSIG.Gct2Frame(filename = paste(output.dir, "temp2.gct", 
        sep = ""))
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dim(m)
    Ns <- length(m[1, ])
    Ng <- length(m[, 1])
    CLS <- MSIG.ReadClsFile(file = paste(output.dir, "temp1.cls", 
        sep = ""))
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    class.list <- CLS$class.list
    if (is.na(match(target.class, class.phen))) 
        stop(c("target class is not phenotype in:", model.cls))
    print("Renaming phenotypes...")
    for (i in 1:length(class.list)) {
        if (class.list[i] == target.class) {
            class.labels[i] <- 1
        }
        else {
            class.list[i] <- "CNTL"
            class.labels[i] <- 0
        }
    }
    col.index <- order(class.labels, decreasing = T)
    for (j in 1:Ng) {
        m[j, ] <- m[j, col.index]
    }
    sample.names <- sample.names[col.index]
    class.labels <- class.labels[col.index]
    class.list <- class.list[col.index]
    class.phen <- c(target.class, "CNTL")
    control.class <- "CNTL"
    print("Executing marker selection...")
    topgs <- ifelse(topgs > floor(Ng/length(class.phen)), floor(Ng/length(class.phen)), 
        topgs)
    sample.molsig.sorted.subset <- matrix(0, nrow = length(class.phen) * 
        topgs, ncol = Ns)
    sample.molsig.sorted.subset.gs <- vector(length = length(class.phen) * 
        topgs, mode = "character")
    sample.molsig.sorted.s2n <- vector(length = length(class.phen) * 
        topgs, mode = "character")
    sample.molsig.sorted.class <- vector(length = length(class.phen) * 
        topgs, mode = "character")
    sample.molsig.sorted.index <- vector(length = length(class.phen) * 
        topgs, mode = "character")
    num.k <- 1
    for (k in class.phen) {
        print(c("Executing marker selection for class:", k))
        class.k.labels <- ifelse(class.list == k, 0, 1)
        col.index <- order(class.k.labels, decreasing = F)
        m1 <- m
        for (j in 1:Ng) {
            m1[j, ] <- m[j, col.index]
        }
        names(m1) <- sample.names[col.index]
        class.k.labels <- class.k.labels[col.index]
        O <- GSEA.GeneRanking(m1, class.k.labels, gs.names, 1, 
            permutation.type = 0, sigma.correction = "GeneCluster", 
            fraction = 1, replace = F, reverse.sign = F)
        order.matrix <- O$order.matrix
        obs.order.matrix <- O$obs.order.matrix
        correl.matrix <- O$s2n.matrix
        obs.correl.matrix <- O$obs.s2n.matrix
        rm(O)
        obs.s2n.orig <- apply(obs.correl.matrix, 1, median)
        obs.index <- order(obs.s2n.orig, decreasing = T)
        obs.s2n <- sort(obs.s2n.orig, decreasing = T)
        sample.molsig.sorted <- m[obs.index, ]
        gs.names.sorted <- gs.names[obs.index]
        start <- (num.k - 1) * topgs + 1
        end <- num.k * topgs
        sample.molsig.sorted.subset[start:end, ] <- sample.molsig.sorted[1:topgs, 
            ]
        sample.molsig.sorted.subset.gs[start:end] <- gs.names.sorted[1:topgs]
        sample.molsig.sorted.s2n[start:end] <- signif(obs.s2n[1:topgs], 
            digits = 3)
        sample.molsig.sorted.class[start:end] <- class.phen[num.k]
        sample.molsig.sorted.index[start:end] <- obs.index[1:topgs]
        num.k <- num.k + 1
    }
    c1 <- c("grey3", "cadetblue", "darkgreen", "chartreuse2", 
        "red1", "darkred", "orange", "blue2", "lightblue", "pink", 
        "coral")
    print("Making heat maps...")
    MSIG.HeatMapPlot.5(V = t(sample.molsig.sorted.subset), row.names = sample.names, 
        col.labels = c(rep(1, topgs), rep(0, topgs)), col.classes = c("C1", 
            "C2"), col.names = sample.molsig.sorted.subset.gs, 
        main = "Original Expression (Norm. Signature genes)", 
        xlab = " ", ylab = " ", sub = gct.file, row.norm = T, 
        cmap.type = 4, rotated.col.labels = T)
    msize <- 2 * topgs
    sig.matrix <- array(0, dim = c(msize, Ns))
    for (k in 1:Ns) {
        sig.matrix[, k] <- rank(sample.molsig.sorted.subset[, 
            k], ties.method = "average")
    }
    MSIG.HeatMapPlot.5(V = t(sig.matrix), row.names = sample.names, 
        col.labels = c(rep(1, topgs), rep(0, topgs)), col.classes = c("C1", 
            "C2"), col.names = sample.molsig.sorted.subset.gs, 
        main = "Signature -- Training", xlab = " ", ylab = " ", 
        sub = gct.file, row.norm = F, cmap.type = 2, rotated.col.labels = T)
    print("Defining signature...")
    t.class.point <- apply(sig.matrix[, class.list == target.class], 
        MARGIN = 1, FUN = mean)
    c.class.point <- apply(sig.matrix[, class.list == control.class], 
        MARGIN = 1, FUN = mean)
    d.t.class <- vector(length = Ns, mode = "numeric")
    d.c.class <- vector(length = Ns, mode = "numeric")
    d.c.t.class <- sum(abs(t.class.point - c.class.point))
    x <- vector(length = Ns, mode = "numeric")
    y <- vector(length = Ns, mode = "numeric")
    for (i in 1:Ns) {
        d.t.class[i] <- sum(abs(t.class.point - sig.matrix[, 
            i]))/d.c.t.class
        d.c.class[i] <- sum(abs(c.class.point - sig.matrix[, 
            i]))/d.c.t.class
        x[i] <- (d.t.class[i]^2 - d.c.class[i]^2 - 1)/(-2)
        y[i] <- sqrt(d.c.class[i]^2 - x[i]^2)
    }
    x.range <- range(c(x, 0, 1))
    y.range <- range(c(y, 0))
    x11(height = 24, width = 30)
    plot(x, y, xlim = x.range, ylim = y.range, type = "n")
    points(0, 0, cex = 2, pch = 21, col = 1, bg = 3)
    points(1, 0, cex = 2, pch = 21, col = 1, bg = 2)
    points(x[class.list == control.class], y[class.list == control.class], 
        cex = 1, pch = 21, col = 1, bg = 3)
    points(x[class.list == target.class], y[class.list == target.class], 
        cex = 1, pch = 21, col = 1, bg = 2)
    target.var <- ifelse(class.list == target.class, 1, 0)
    if (model.type == "Bayesian") {
        if (link.function == "logit") {
            reg.model <- MCMClogit(target.var ~ x, burnin = burnin.iter, 
                mcmc = mcmc.iter, bayes.resid = T)
        }
        else if (link.function == "probit") {
            reg.model <- MCMCprobit(target.var ~ x, burnin = burnin.iter, 
                mcmc = mcmc.iter, bayes.resid = T)
        }
        else {
            stop("Unknown link function")
        }
    }
    else if (model.type == "Classic") {
        if (link.function == "logit") {
            reg.model <- glm(target.var ~ x, family = binomial("logit"))
        }
        else if (link.function == "probit") {
            reg.model <- glm(target.var ~ x, family = binomial("probit"))
        }
        else {
            stop("Unknown link function")
        }
    }
    else {
        stop("Unknown model type")
    }
    if (model.type == "Bayesian") {
        beta0 <- reg.model[, 1]
        beta1 <- reg.model[, 2]
        print(c("beta0=", beta0, " beta1=", beta1))
        prob.i <- matrix(0, nrow = Ns, ncol = 3)
    }
    else if (model.type == "Classic") {
        beta0 <- reg.model[[1]][1]
        beta1 <- reg.model[[1]][2]
        print(c("beta0=", beta0, " beta1=", beta1))
        prob.i <- matrix(0, nrow = Ns, ncol = 3)
    }
    else {
        stop("Unknown model type")
    }
    for (i in 1:Ns) {
        if (link.function == "logit") {
            p.vec <- (exp(beta0 + beta1 * x[i])/(1 + exp(beta0 + 
                beta1 * x[i])))
        }
        else if (link.function == "probit") {
            p.vec <- (erf(beta0 + beta1 * x[i]) + 1)/2
        }
        else {
            stop("Unknown link function")
        }
        prob.i[i, 1] <- quantile(p.vec, probs = 0.5)
        prob.i[i, 2] <- quantile(p.vec, probs = 0.05)
        prob.i[i, 3] <- quantile(p.vec, probs = 0.95)
    }
    xmin <- min(x)
    xmax <- max(x)
    range.x <- xmax - xmin
    prob.m <- matrix(0, nrow = 1000, ncol = 3)
    x.m <- vector(length = 1000, mode = "numeric")
    for (k in 1:1000) {
        x.m[k] <- xmin + k * (range.x/1000)
        if (link.function == "logit") {
            p.vec <- (exp(beta0 + beta1 * x.m[k])/(1 + exp(beta0 + 
                beta1 * x.m[k])))
        }
        else if (link.function == "probit") {
            p.vec <- (erf(beta0 + beta1 * x.m[k]) + 1)/2
        }
        else {
            stop("Unknown link function")
        }
        prob.m[k, 1] <- quantile(p.vec, probs = 0.5)
        prob.m[k, 2] <- quantile(p.vec, probs = 0.05)
        prob.m[k, 3] <- quantile(p.vec, probs = 0.95)
    }
    istar <- which.min(abs(0.5 - prob.m[, 1]))
    istar <- xmin + istar * (range.x/1000)
    x.index <- order(x, decreasing = F)
    x.order <- x[x.index]
    prob.i.order <- prob.i[x.index, ]
    target.var.order <- target.var[x.index]
    target.var.order <- ifelse(target.var.order == 1, c1[1], 
        c1[2])
    target.var <- ifelse(target.var == 1, c1[1], c1[2])
    x11(height = 7, width = 9.5)
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 
        1), heights = 1, respect = FALSE)
    plot(x.order, prob.i.order[, 1], sub = gct.file, pch = 20, 
        ylim = c(-0.2, 1.07), xlim = c(-0.1, 1.1), col = 0, cex = 2, 
        xlab = "Activation Index", ylab = "Probability")
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
    p.vec <- rep(21, 21)
    c.vec <- c1
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 1.2, pt.cex = 2)
    activation.indicator <- ifelse(prob.i[, 1] >= 0.5, 1, 0)
    activation.indicator <- ifelse((prob.i[, 1] >= no.call.range.max) | 
        (prob.i[, 1] <= no.call.range.min), activation.indicator, 
        0.5)
    z <- rbind(prob.i[, 1], activation.indicator)
    print("z:")
    print(z)
    print("class.labels:")
    print(class.labels)
    print("class.phen:")
    print(class.phen)
    print("sample.names:")
    print(sample.names)
    MSIG.HeatMapPlot.5(V = z, row.names = c("P(A)", "A"), col.labels = class.labels, 
        col.classes = class.phen, col.names = sample.names, main = "Activation Assessment (training set)", 
        xlab = " ", ylab = " ", sub = gct.file, row.norm = F, 
        cmap.type = 2)
    filename <- paste(database.dir, identifier, ".msig.params", 
        sep = "")
    write(paste("seed", seed, sep = "\t"), file = filename, append = F, 
        ncolumns = 3)
    write(paste("topgs", topgs, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("link.function", link.function, sep = "\t"), 
        file = filename, append = T, ncolumns = 3)
    write(paste("model.type", model.type, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("burnin.iter", burnin.iter, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("mcmc.iter", mcmc.iter, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("col.target", col.target, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("col.control", col.control, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("no.call.range.max", no.call.range.max, sep = "\t"), 
        file = filename, append = T, ncolumns = 3)
    write(paste("no.call.range.min", no.call.range.min, sep = "\t"), 
        file = filename, append = T, ncolumns = 3)
    write(paste("beta0", beta0, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("beta1", beta1, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    write(paste("target class", target.class, sep = "\t"), file = filename, 
        append = T, ncolumns = 3)
    lset <- length(sample.molsig.sorted.subset.gs)
    filename <- paste(database.dir, identifier, ".msig.up.gmt", 
        sep = "")
    gene.set <- sample.molsig.sorted.subset.gs[seq(1, lset/2)]
    output.line <- paste(gene.set, sep = "\t", collapse = "\t")
    gene.set.name <- paste(identifier, ".up", sep = "")
    output.line <- paste(gene.set.name, gene.set.name, output.line, 
        sep = "\t", collapse = "")
    write(noquote(output.line), file = filename, append = F, 
        ncolumns = length(gene.set) + 2)
    filename <- paste(database.dir, identifier, ".msig.dn.gmt", 
        sep = "")
    gene.set <- sample.molsig.sorted.subset.gs[seq(lset/2 + 1, 
        lset)]
    output.line <- paste(gene.set, sep = "\t", collapse = "\t")
    gene.set.name <- paste(identifier, ".dn", sep = "")
    output.line <- paste(gene.set.name, gene.set.name, output.line, 
        sep = "\t", collapse = "")
    write(noquote(output.line), file = filename, append = F, 
        ncolumns = length(gene.set) + 2)
    filename <- paste(database.dir, identifier, ".msig.gct", 
        sep = "")
    z <- data.frame(sample.molsig.sorted.subset)
    names(z) <- sample.names
    row.names(z) <- sample.molsig.sorted.subset.gs
    write.gct(gct.data.frame = z, descs = sample.molsig.sorted.subset.gs, 
        filename = filename)
    class.labels <- ifelse(class.labels == 1, 1, 2)
    print("class.labels:")
    print(class.labels)
    filename <- paste(database.dir, identifier, ".msig.cls", 
        sep = "")
    write.cls(class.v = class.labels, phen = class.phen, filename = filename)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
