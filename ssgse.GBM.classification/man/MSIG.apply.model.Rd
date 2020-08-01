\name{MSIG.apply.model}
\alias{MSIG.apply.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.apply.model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.apply.model(gct.file, cls.file, phen.annot.file = NULL, output.dir, database.dir, identifiers, column.subset = "ALL", column.sel.type = "samples", thres = "NULL", ceil = "NULL", shift = "NULL", fold = 1, delta = 0, norm = 6, no.call.range.max = NULL, no.call.range.min = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gct.file}{
%%     ~~Describe \code{gct.file} here~~
}
  \item{cls.file}{
%%     ~~Describe \code{cls.file} here~~
}
  \item{phen.annot.file}{
%%     ~~Describe \code{phen.annot.file} here~~
}
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
}
  \item{database.dir}{
%%     ~~Describe \code{database.dir} here~~
}
  \item{identifiers}{
%%     ~~Describe \code{identifiers} here~~
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
  \item{no.call.range.max}{
%%     ~~Describe \code{no.call.range.max} here~~
}
  \item{no.call.range.min}{
%%     ~~Describe \code{no.call.range.min} here~~
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
function (gct.file, cls.file, phen.annot.file = NULL, output.dir, 
    database.dir, identifiers, column.subset = "ALL", column.sel.type = "samples", 
    thres = "NULL", ceil = "NULL", shift = "NULL", fold = 1, 
    delta = 0, norm = 6, no.call.range.max = NULL, no.call.range.min = NULL) 
{
    print(c("Processing test file: ", gct.file))
    O <- MSIG.Subset.Dataset(input.ds = gct.file, input.cls = cls.file, 
        column.subset = column.subset, column.sel.type = column.sel.type, 
        row.subset = "ALL", output.ds = paste(output.dir, "temp2.gct", 
            sep = ""), output.cls = paste(output.dir, "temp2.cls", 
            sep = ""))
    O <- MSIG.Preprocess.Dataset(input.ds = paste(output.dir, 
        "temp2.gct", sep = ""), output.ds = paste(output.dir, 
        "temp3.gct", sep = ""), thres = thres, ceil = ceil, normalization = "NULL")
    dataset <- MSIG.Gct2Frame(filename = paste(output.dir, "temp3.gct", 
        sep = ""))
    m.test <- data.matrix(dataset$ds)
    gs.names.test <- dataset$row.names
    gs.descs.test <- dataset$descs
    sample.names.test <- dataset$names
    Ns.test <- length(m.test[1, ])
    Ng.test <- length(m.test[, 1])
    CLS <- MSIG.ReadClsFile(file = paste(output.dir, "temp2.cls", 
        sep = ""))
    class.labels.test <- CLS$class.v
    class.phen.test <- CLS$phen
    class.list.test <- CLS$class.list
    for (sig in identifiers) {
        filename <- paste(database.dir, sig, ".msig.params", 
            sep = "")
        temp <- readLines(filename)
        seed <- as.numeric(noquote(unlist(strsplit(temp[[1]], 
            "\t")))[2])
        topgs <- as.numeric(noquote(unlist(strsplit(temp[[2]], 
            "\t")))[2])
        link.function <- unlist(strsplit(temp[[3]], "\t"))[2]
        model.type <- unlist(strsplit(temp[[4]], "\t"))[2]
        burnin.iter <- as.numeric(noquote(unlist(strsplit(temp[[5]], 
            "\t")))[2])
        mcmc.iter <- as.numeric(noquote(unlist(strsplit(temp[[6]], 
            "\t")))[2])
        col.target <- unlist(strsplit(temp[[7]], "\t"))[2]
        col.control <- unlist(strsplit(temp[[8]], "\t"))[2]
        no.call.r.max <- as.numeric(noquote(unlist(strsplit(temp[[9]], 
            "\t")))[2])
        no.call.r.min <- as.numeric(noquote(unlist(strsplit(temp[[10]], 
            "\t")))[2])
        beta0.train <- as.numeric(noquote(unlist(strsplit(temp[[11]], 
            "\t")))[2])
        beta1.train <- as.numeric(noquote(unlist(strsplit(temp[[12]], 
            "\t")))[2])
        target.class <- unlist(strsplit(temp[[13]], "\t"))[2]
        c1 <- c(col.target, col.control)
        if (is.null(no.call.range.max)) {
            no.call.range.max <- no.call.r.max
        }
        if (is.null(no.call.range.min)) {
            no.call.range.min <- no.call.r.min
        }
        filename <- paste(database.dir, sig, ".msig.gct", sep = "")
        dataset <- MSIG.Gct2Frame(filename = filename)
        sample.molsig.sorted.subset <- dataset$ds
        Ns <- length(sample.molsig.sorted.subset[1, ])
        msize.all <- length(sample.molsig.sorted.subset[, 1])
        sample.molsig.sorted.subset.gs <- dataset$row.names
        sample.names <- dataset$names
        filename <- paste(database.dir, sig, ".msig.gct", sep = "")
        dataset <- MSIG.Gct2Frame(filename = filename)
        sample.molsig.sorted.subset <- dataset$ds
        Ns <- length(sample.molsig.sorted.subset[1, ])
        msize.all <- length(sample.molsig.sorted.subset[, 1])
        sample.molsig.sorted.subset.gs <- dataset$row.names
        sample.names <- dataset$names
        filename <- paste(database.dir, sig, ".msig.cls", sep = "")
        CLS <- MSIG.ReadClsFile(file = filename)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
        class.list <- CLS$class.list
        for (i in 1:length(class.list)) {
            if (class.list[i] == target.class) {
                class.labels[i] <- 1
            }
            else {
                class.list[i] <- "CNTL"
                class.labels[i] <- 0
            }
        }
        print(c("Target class:", target.class))
        print(c("Class labels:", class.labels))
        col.index <- order(class.labels, decreasing = T)
        for (j in 1:msize.all) {
            sample.molsig.sorted.subset[j, ] <- sample.molsig.sorted.subset[j, 
                col.index]
        }
        sample.names <- sample.names[col.index]
        class.labels <- class.labels[col.index]
        class.list <- class.list[col.index]
        class.phen <- c(target.class, "CNTL")
        control.class <- "CNTL"
        gs.names2 <- intersect(sample.molsig.sorted.subset.gs, 
            gs.names.test)
        locations <- match(gs.names2, gs.names.test, nomatch = 0)
        m.test2 <- m.test[locations, ]
        locations2 <- match(gs.names2, sample.molsig.sorted.subset.gs)
        m.train <- sample.molsig.sorted.subset[locations2, ]
        print(c("Matched signature and test set: overlap=", length(gs.names2), 
            " Total original signature size= ", length(sample.molsig.sorted.subset.gs)))
        msize <- length(locations)
        sig.matrix <- array(0, dim = c(msize, Ns))
        sig.matrix.test <- array(0, dim = c(msize, Ns.test))
        for (k in 1:Ns) {
            sig.matrix[, k] <- rank(m.train[, k], ties.method = "average")
        }
        for (k in 1:Ns.test) {
            sig.matrix.test[, k] <- rank(m.test2[, k], ties.method = "average")
        }
        sig.matrix.all <- cbind(sig.matrix, sig.matrix.test)
        sample.names.all <- c(sample.names, sample.names.test)
        MSIG.HeatMapPlot.5(V = t(sig.matrix.all), row.names = sample.names.all, 
            col.labels = rep(1, msize), col.classes = "C", col.names = gs.names2, 
            main = paste(sig, gct.file, sep = " / "), xlab = " ", 
            ylab = " ", row.norm = F, cmap.type = 2)
        t.class.point <- apply(sig.matrix[, class.list == target.class], 
            MARGIN = 1, FUN = mean)
        c.class.point <- apply(sig.matrix[, class.list == control.class], 
            MARGIN = 1, FUN = mean)
        d.t.class <- vector(length = Ns, mode = "numeric")
        d.c.class <- vector(length = Ns, mode = "numeric")
        d.c.t.class <- sum(abs(t.class.point - c.class.point))
        x <- vector(length = Ns, mode = "numeric")
        y <- vector(length = Ns, mode = "numeric")
        d.t.class.test <- vector(length = Ns.test, mode = "numeric")
        d.c.class.test <- vector(length = Ns.test, mode = "numeric")
        x.test <- vector(length = Ns.test, mode = "numeric")
        y.test <- vector(length = Ns.test, mode = "numeric")
        for (i in 1:Ns) {
            d.t.class[i] <- sum(abs(t.class.point - sig.matrix[, 
                i]))/d.c.t.class
            d.c.class[i] <- sum(abs(c.class.point - sig.matrix[, 
                i]))/d.c.t.class
            x[i] <- (d.t.class[i]^2 - d.c.class[i]^2 - 1)/(-2)
            y[i] <- sqrt(d.c.class[i]^2 - x[i]^2)
        }
        print(c("Creating regression signature model using overlap..."))
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
        print(c("beta0 train=", beta0.train, " beta0=", beta0))
        print(c("beta1 train=", beta1.train, " beta1=", beta1))
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
                nstop("Unknown link function")
            }
            prob.m[k, 1] <- quantile(p.vec, probs = 0.5)
            prob.m[k, 2] <- quantile(p.vec, probs = 0.05)
            prob.m[k, 3] <- quantile(p.vec, probs = 0.95)
        }
        istar <- which.min(abs(0.5 - prob.m[, 1]))
        istar <- xmin + istar * (range.x/1000)
        for (i in 1:Ns.test) {
            d.t.class.test[i] <- sum(abs(t.class.point - sig.matrix.test[, 
                i]))/d.c.t.class
            d.c.class.test[i] <- sum(abs(c.class.point - sig.matrix.test[, 
                i]))/d.c.t.class
            x.test[i] <- (d.t.class.test[i]^2 - d.c.class.test[i]^2 - 
                1)/(-2)
            y.test[i] <- sqrt(d.c.class.test[i]^2 - x.test[i]^2)
        }
        x.range <- range(c(x, x.test, 0, 1))
        y.range <- range(c(y, y.test, 0))
        x11(height = 24, width = 30)
        plot(x, y, xlim = x.range, ylim = y.range, type = "n", 
            main = sig, sub = gct.file)
        points(0, 0, cex = 2, pch = 21, col = 1, bg = 3)
        points(1, 0, cex = 2, pch = 21, col = 1, bg = 2)
        points(x[class.list == control.class], y[class.list == 
            control.class], cex = 1, pch = 21, col = 1, bg = 3)
        points(x[class.list == target.class], y[class.list == 
            target.class], cex = 1, pch = 21, col = 1, bg = 2)
        k <- 1
        for (i in class.list.test) {
            points(x.test[class.list.test == i], y.test[class.list.test == 
                i], cex = 1, pch = 22, col = 1, bg = k\%\%5)
            k <- k + 1
        }
        prob.i.test <- matrix(0, nrow = Ns.test, ncol = 3)
        for (i in 1:Ns.test) {
            if (link.function == "logit") {
                p.vec.test <- (exp(beta0 + beta1 * x.test[i])/(1 + 
                  exp(beta0 + beta1 * x.test[i])))
            }
            else if (link.function == "probit") {
                p.vec.test <- (erf(beta0 + beta1 * x.test[i]) + 
                  1)/2
            }
            else {
                stop("Unknown link function")
            }
            prob.i.test[i, 1] <- quantile(p.vec.test, probs = 0.5)
            prob.i.test[i, 2] <- quantile(p.vec.test, probs = 0.05)
            prob.i.test[i, 3] <- quantile(p.vec.test, probs = 0.95)
        }
        x.index <- order(x.test, decreasing = F)
        x.order.test <- x.test[x.index]
        prob.i.order.test <- prob.i.test[x.index, ]
        class.list.test.order <- class.list.test[x.index]
        x11(height = 7, width = 9.5)
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3.75, 
            1), heights = 1, respect = FALSE)
        plot(x.order.test, prob.i.order.test[, 1], sub = gct.file, 
            pch = 20, ylim = c(-0.05, 1.07), main = sig, xlim = c(-0.1, 
                1.1), col = 0, cex.axis = 1.35, cex = 3, cex.lab = 1.35, 
            xlab = "Activation Index", ylab = "Probability")
        points(x.m, prob.m[, 1], type = "l", lwd = 2, col = 1, 
            lty = 1, cex = 1)
        points(x.m, prob.m[, 2], type = "l", col = 4, lty = 1, 
            cex = 1)
        points(x.m, prob.m[, 3], type = "l", col = 4, lty = 1, 
            cex = 1)
        arrows(x.order.test, prob.i.order.test[, 2], x.order.test, 
            prob.i.order.test[, 3], col = 4, angle = 90, code = 3, 
            length = 0)
        range.x <- range(x.order.test)
        points(range.x, c(0.5, 0.5), type = "l", lty = 3, col = 1, 
            lwd = 2)
        points(c(istar, istar), c(-0.07, 1.07), type = "l", lty = 3, 
            col = 1, lwd = 2)
        k <- 1
        for (i in class.list.test) {
            points(x.order.test[class.list.test.order == i], 
                prob.i.order.test[class.list.test.order == i, 
                  1], pch = 21, bg = k\%\%5, col = 1, cex = 2)
            k <- k + 1
        }
        leg.txt <- unique(class.list.test.order)
        p.vec <- rep(21, length(unique(class.list.test.order)))
        c.vec <- rep(seq(1, 5), length(unique(class.list.test.order)))
        par(mar = c(0, 0, 0, 0))
        plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), 
            axes = F, type = "n", xlab = "", ylab = "")
        legend(x = 0, y = 0.8, legend = leg.txt, bty = "n", xjust = 0, 
            yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
            cex = 1.2, pt.cex = 2)
        activation.indicator <- ifelse(prob.i.test[, 1] >= 0.5, 
            1, 0)
        activation.indicator <- ifelse((prob.i.test[, 1] >= no.call.range.max) | 
            (prob.i.test[, 1] <= no.call.range.min), activation.indicator, 
            0.5)
        if (!is.null(phen.annot.file)) {
            filename <- phen.annot.file
            dataset <- MSIG.Gct2Frame(filename = filename)
            phen.annot <- data.matrix(dataset$ds)
            phen.annot.gs <- dataset$row.names
            for (i in 1:length(phen.annot[, 1])) {
                phen.annot[i, ] <- (phen.annot[i, ] - min(phen.annot[i, 
                  ]))/(max(phen.annot[i, ]) - min(phen.annot[i, 
                  ]))
            }
            z <- rbind(prob.i.test[, 1], activation.indicator, 
                phen.annot)
            p.lab <- c(paste("P(", sig, ")", sep = ""), paste("A(", 
                sig, ")", sep = ""), phen.annot.gs)
        }
        else {
            z <- rbind(prob.i.test[, 1], activation.indicator)
            p.lab <- c(paste("P(", sig, ")", sep = ""), paste("A(", 
                sig, ")", sep = ""))
        }
        MSIG.HeatMapPlot.5(V = z, row.names = p.lab, col.labels = class.labels.test, 
            col.classes = class.phen.test, col.names = sample.names.test, 
            main = paste(sig, " Activation on Test", sep = ""), 
            xlab = " ", ylab = " ", sub = gct.file, row.norm = F, 
            cmap.type = 3, rotated.col.labels = T)
        if (sig == identifiers[[1]]) {
            z.all <- prob.i.test[, 1]
            z.act.all <- activation.indicator
            if (!is.null(phen.annot.file)) {
                phen.annot.all <- phen.annot
                phen.annot.gs.all <- phen.annot.gs
            }
            p.lab.all <- paste("P(", sig, ")", sep = "")
            p.act.lab.all <- paste("A(", sig, ")", sep = "")
        }
        else {
            z.all <- rbind(z.all, prob.i.test[, 1])
            z.act.all <- rbind(z.act.all, activation.indicator)
            p.lab.all <- c(p.lab.all, paste("P(", sig, ")", sep = ""))
            p.act.lab.all <- c(p.act.lab.all, paste("A(", sig, 
                ")", sep = ""))
        }
    }
    if (!is.null(phen.annot.file)) {
        z.all <- rbind(z.all, phen.annot)
        z.act.all <- rbind(z.act.all, phen.annot)
        p.lab.all <- c(p.lab.all, phen.annot.gs)
        p.act.lab.all <- c(p.act.lab.all, phen.annot.gs)
    }
    print(c("dim z.all=", dim(z.all)))
    MSIG.HeatMapPlot.5(V = z.all, row.names = p.lab.all, col.labels = class.labels.test, 
        col.classes = class.phen.test, col.names = sample.names.test, 
        main = " ", xlab = " ", ylab = " ", sub = gct.file, row.norm = F, 
        cmap.type = 2, rotated.col.labels = T)
    MSIG.HeatMapPlot.5(V = z.act.all, row.names = p.act.lab.all, 
        col.labels = class.labels.test, col.classes = class.phen.test, 
        col.names = sample.names.test, main = " ", xlab = " ", 
        ylab = " ", sub = gct.file, row.norm = F, cmap.type = 2, 
        rotated.col.labels = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
