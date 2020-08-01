\name{MSIG.Evaluate.Projection}
\alias{MSIG.Evaluate.Projection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.Evaluate.Projection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.Evaluate.Projection(input.ds, input.cls, model.set, prediction.results.file, prediction.matrix.file, train.pred.plot, test.pred.plot, pred.2D.plot, col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"), non.interactive.run = F, use.feature.names = F, nchar.phen = 3, high.conf.thres = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{model.set}{
%%     ~~Describe \code{model.set} here~~
}
  \item{prediction.results.file}{
%%     ~~Describe \code{prediction.results.file} here~~
}
  \item{prediction.matrix.file}{
%%     ~~Describe \code{prediction.matrix.file} here~~
}
  \item{train.pred.plot}{
%%     ~~Describe \code{train.pred.plot} here~~
}
  \item{test.pred.plot}{
%%     ~~Describe \code{test.pred.plot} here~~
}
  \item{pred.2D.plot}{
%%     ~~Describe \code{pred.2D.plot} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{use.feature.names}{
%%     ~~Describe \code{use.feature.names} here~~
}
  \item{nchar.phen}{
%%     ~~Describe \code{nchar.phen} here~~
}
  \item{high.conf.thres}{
%%     ~~Describe \code{high.conf.thres} here~~
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
function (input.ds, input.cls, model.set, prediction.results.file, 
    prediction.matrix.file, train.pred.plot, test.pred.plot, 
    pred.2D.plot, col = c("grey3", "cadetblue", "darkgreen", 
        "chartreuse2", "red1", "darkred", "orange", "blue2", 
        "lightblue", "pink", "coral"), non.interactive.run = F, 
    use.feature.names = F, nchar.phen = 3, high.conf.thres = 0.7) 
{
    print(c("Running MSIG.Evaluate.Projection... on:", input.ds))
    library(e1071)
    library(tree)
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dim(m)
    Ns <- length(m[1, ])
    N <- length(m[, 1])
    CLS <- MSIG.ReadClsFile(file = input.cls)
    class.labels <- CLS$class.v
    class.list <- CLS$class.list
    class.phen <- CLS$phen
    num.classes <- length(class.phen)
    print("Reading dataset completed...")
    class.list2 <- vector(length = Ns, mode = "numeric")
    for (i in 1:Ns) {
        class.list2[i] <- substr(class.list[i], 1, nchar.phen)
    }
    class.phen2 <- vector(length = num.classes, mode = "character")
    for (i in 1:num.classes) {
        class.phen2[i] <- substr(class.phen[i], 1, 3)
    }
    true.num.classes <- length(table(class.phen2))
    class.labels2 <- match(class.list2, class.phen2)
    m.train <- m[, model.set]
    n.train <- length(model.set)
    num.samples.train <- n.train
    sample.names.train <- as.factor(sample.names[model.set])
    class.list.train <- class.list2[model.set]
    class.phen.train <- unique(class.list.train)
    class.labels.train <- class.labels2[model.set]
    if (Ns - length(model.set) > 0) {
        m.test <- as.matrix(m[, -model.set])
        n.test <- length(m.test[1, ])
        sample.names.test <- as.factor(sample.names[-model.set])
        class.list.test <- class.list2[-model.set]
        class.phen.test <- unique(class.list.test)
        class.labels.test <- class.labels2[-model.set]
    }
    print("Building SVM model...")
    svm.model <- svm(x = t(m.train), y = class.list.train, scale = F, 
        tolerance = 1e-04, type = "C-classification", kernel = "linear", 
        cost = 3, probability = T)
    print("Computing train set predictions...")
    train.pred <- predict(object = svm.model, newdata = t(m.train), 
        decision.values = T, probability = T)
    dec.vals.train <- attr(train.pred, "decision.values")
    prob.train <- signif(attr(train.pred, "probabilities"), digits = 2)
    confidence.vector <- apply(prob.train, 1, max)
    confidence.call <- ifelse(confidence.vector < high.conf.thres, 
        " L ", " H ")
    error.call <- ifelse(class.list.train == as.character(train.pred), 
        "   ", " * ")
    col.symbols.train <- paste(confidence.call, error.call)
    class.names <- names(data.frame(prob.train))
    bscore <- vector(length = n.train, mode = "numeric")
    for (i in 1:n.train) {
        bscore[i] <- 0
        for (k in 1:length(prob.train[1, ])) {
            if (class.list.train[i] == class.names[k]) {
                bscore[i] <- bscore[i] + (1 - prob.train[i, k])^2
            }
            else {
                bscore[i] <- bscore[i] + prob.train[i, k]^2
            }
        }
        bscore[i] <- signif(bscore[i], digits = 2)
    }
    Brier.train <- signif(mean(bscore), digits = 2)
    train.results <- data.frame(cbind(as.character(sample.names.train), 
        class.list.train, as.character(train.pred), error.call, 
        confidence.call, prob.train, bscore))
    names(train.results)[1] <- "Train Sample Name"
    names(train.results)[2] <- "Actual"
    names(train.results)[3] <- "Predicted"
    names(train.results)[4] <- "Error (*)"
    names(train.results)[5] <- "Conf (H/L)"
    names(train.results)[6 + length(class.phen.train)] <- "Brier score"
    print(train.results)
    print(c("Brier score (Train) = ", Brier.train))
    write("Training Results \n", file = prediction.results.file, 
        append = F)
    write.table(train.results, file = prediction.results.file, 
        append = T, quote = F, row.names = F, sep = "\t")
    write(c("\n\n Brier score (Train) = ", Brier.train), file = prediction.results.file, 
        append = T)
    conf.table.train <- table(class.list.train, train.pred)
    print(conf.table.train)
    write("\n\n Confusion Matrix (Train) \n", file = prediction.results.file, 
        append = T)
    write.table(conf.table.train, file = prediction.results.file, 
        append = T, quote = F, row.names = F, sep = "\t")
    height <- ifelse(length(class.phen.train) > 50, 15, 0.3 * 
        length(class.phen.train) + 4.8)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- train.pred.plot
            x11(height = height, width = 15)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(train.pred.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(train.pred.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(train.pred.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
    }
    MSIG.HeatMapPlot.2(V = t(prob.train), row.names = class.names, 
        col.symbols = col.symbols.train, col.labels = class.list.train, 
        col.names = as.character(sample.names.train), main = "Train Samples Predictions", 
        sub = " ", xlab = " ", ylab = " ", row.norm = F, cmap.type = 2)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = train.pred.plot, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    print("Building SVM model completed. Predicting test data...")
    if (Ns - length(model.set) > 0) {
        test.pred <- predict(object = svm.model, newdata = t(m.test), 
            decision.values = T, probability = T)
        dec.vals.test <- attr(test.pred, "decision.values")
        prob.test <- signif(attr(test.pred, "probabilities"), 
            digits = 2)
        confidence.vector <- apply(prob.test, 1, max)
        confidence.call <- ifelse(confidence.vector < high.conf.thres, 
            " L ", " H ")
        error.call <- ifelse(class.list.test == as.character(test.pred), 
            "   ", " * ")
        col.symbols.test <- paste(confidence.call, error.call)
        class.names <- names(data.frame(prob.test))
        bscore <- vector(length = n.test, mode = "numeric")
        for (i in 1:n.test) {
            bscore[i] <- 0
            for (k in 1:length(prob.test[1, ])) {
                if (class.list.test[i] == class.names[k]) {
                  bscore[i] <- bscore[i] + (1 - prob.test[i, 
                    k])^2
                }
                else {
                  bscore[i] <- bscore[i] + prob.test[i, k]^2
                }
            }
            bscore[i] <- signif(bscore[i], digits = 2)
        }
        Brier.test <- signif(mean(bscore), digits = 2)
        test.results <- data.frame(cbind(as.character(sample.names.test), 
            class.list.test, as.character(test.pred), error.call, 
            confidence.call, prob.test, bscore))
        names(test.results)[1] <- "Test Sample Name"
        names(test.results)[2] <- "Actual"
        names(test.results)[3] <- "Predicted"
        names(test.results)[4] <- "Error (*)"
        names(test.results)[5] <- "Conf (H/L)"
        names(test.results)[6 + length(class.phen.train)] <- "Brier score"
        print(test.results)
        print(c("Brier score (Test) = ", Brier.test))
        write("\n Test Results \n", file = prediction.results.file, 
            append = T)
        write.table(test.results, file = prediction.results.file, 
            append = T, quote = F, row.names = F, sep = "\t")
        write(c("\n\n Brier score (Test) = ", Brier.test), file = prediction.results.file, 
            append = T)
        conf.table.test <- table(class.list.test, test.pred)
        print(conf.table.test)
        write("\n\n Confusion Matrix (Test) \n", file = prediction.results.file, 
            append = T)
        write.table(conf.table.test, file = prediction.results.file, 
            append = T, quote = F, row.names = F, sep = "\t")
        height <- ifelse(length(class.phen.train) > 50, 15, 0.3 * 
            length(class.phen.train) + 4.8)
        if (non.interactive.run == F) {
            if (.Platform$OS.type == "windows") {
                plot.filename <- test.pred.plot
                x11(height = height, width = 15)
            }
            else if (.Platform$OS.type == "unix") {
                plot.filename <- paste(test.pred.plot, ".pdf", 
                  sep = "", collapse = "")
                pdf(file = plot.filename, height = height, width = 15)
            }
        }
        else {
            if (.Platform$OS.type == "unix") {
                plot.filename <- paste(test.pred.plot, ".pdf", 
                  sep = "", collapse = "")
                pdf(file = plot.filename, height = height, width = 15)
            }
            else if (.Platform$OS.type == "windows") {
                plot.filename <- paste(test.pred.plot, ".pdf", 
                  sep = "", collapse = "")
                pdf(file = plot.filename, height = height, width = 15)
            }
        }
        MSIG.HeatMapPlot.2(V = t(prob.test), row.names = names(test.results)[seq(6, 
            6 + length(class.phen.train) - 1)], col.symbols = col.symbols.test, 
            col.labels = class.list.test, col.names = as.character(sample.names.test), 
            main = "Test Samples Predictions", sub = " ", xlab = " ", 
            ylab = " ", row.norm = F, cmap.type = 2)
        if (non.interactive.run == F) {
            if (.Platform$OS.type == "windows") {
                savePlot(filename = test.pred.plot, type = "jpeg", 
                  device = dev.cur())
            }
            else if (.Platform$OS.type == "unix") {
                dev.off()
            }
        }
        else {
            dev.off()
        }
    }
    print(c("dim train=", dim(prob.train)))
    print(c("dim test=", dim(prob.test)))
    V <- cbind(t(prob.train), t(prob.test))
    W <- data.frame(V)
    names(W) <- c(as.character(sample.names.train), as.character(sample.names.test))
    row.names(W) <- names(train.results)[seq(6, 6 + length(class.phen.train) - 
        1)]
    print(W)
    write.gct(gct.data.frame = W, descs = row.names(W), filename = prediction.matrix.file)
    print("Done predicting test data...")
    pca <- prcomp(t(m.train), retx = TRUE, center = TRUE, scale. = TRUE)
    S1 <- pca$x[, 1]
    S2 <- pca$x[, 2]
    S3 <- pca$x[, 3]
    X1 <- pca$rotation[, 1]
    X2 <- pca$rotation[, 2]
    X3 <- pca$rotation[, 3]
    c0 <- col
    c1 <- colors()[match(c0, colors())]
    color <- c1[class.labels]
    height <- 15
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- pred.2D.plot
            x11(height = height, width = 15)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 
        1), heights = 1, respect = FALSE)
    if (Ns - length(model.set) > 0) {
        test.scores <- predict(pca, t(m.test))
        S1 <- c(pca$x[, 1], test.scores[, 1])
        S2 <- c(pca$x[, 2], test.scores[, 2])
        S3 <- c(pca$x[, 3], test.scores[, 3])
    }
    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1/max.X
    X2 <- max.S * X2/max.X
    max.A <- max(max.S, max.X)
    num.samples <- length(S1)
    plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), 
        type = "n", main = "  ", sub = input.ds)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            color.code <- c1[class.labels[j]]
        }
        points(S1[j], S2[j], pch = 21, type = "p", cex = 2, bg = color.code, 
            col = "black")
    }
    for (j in 1:N) {
        x.coor <- X1[j] * 0.925
        y.coor <- X2[j] * 0.925
        arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, 
            angle = 20, col = "grey50")
        if (use.feature.names == FALSE) {
            leg.txt <- paste("F", j, sep = "")
        }
        else {
            leg.txt <- gs.names[j]
        }
        text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, 
            offset = 0.5, vfont = NULL, cex = 1, col = "grey50")
    }
    max.points <- 1000
    max.d <- apply(m, 1, max)
    grid <- matrix(0, nrow = N, ncol = max.points)
    for (i in 1:max.points) {
        for (j in 1:N) {
            grid[j, i] <- runif(1, 0, 1.2 * max.d[j])
        }
    }
    grid.pred <- predict(object = svm.model, newdata = t(grid), 
        decision.values = T, probability = T)
    prob.test <- attr(grid.pred, "probabilities")
    grid.df <- data.frame(t(grid))
    names(grid.df) <- gs.names
    test.scores2 <- predict(pca, grid.df)
    S1 <- test.scores2[, 1]
    S2 <- test.scores2[, 2]
    for (j in 1:max.points) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[match(grid.pred[j], 
                class.list2)] + 1]
        }
        else {
            color.code <- c1[class.labels[match(grid.pred[j], 
                class.list2)]]
        }
        points(S1[j], S2[j], pch = 21, type = "p", cex = 1, bg = color.code, 
            col = color.code)
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- rep(21, n.phen)
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = 2, pt.cex = 3)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = pred.2D.plot, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    MSIG.HeatMapPlot.3(V = t(prob.train), row.names = class.names, 
        col.labels = class.labels.train, col.names = as.character(error.call), 
        col.classes = class.names, phen.cmap = col[1:length(class.names)], 
        main = "Train Samples Predictions", sub = " ", xlab = " ", 
        ylab = " ", row.norm = F, cmap.type = 2)
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- col[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = symbol.scaling * 1.25, pt.cex = symbol.scaling * 
            2.5)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = train.pred.plot, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    print("Building SVM model completed. Predicting test data...")
    if (Ns - length(model.set) > 0) {
        test.pred <- predict(object = svm.model, newdata = t(m.test), 
            decision.values = T, probability = T)
        dec.vals.test <- attr(test.pred, "decision.values")
        prob.test <- signif(attr(test.pred, "probabilities"), 
            digits = 2)
        confidence.vector <- vector(length = n.test, mode = "numeric")
        bscore <- vector(length = n.test, mode = "numeric")
        max.k <- length(prob.train[1, ])
        random.pred.conf <- ((max.k - 1)/max.k)^2 + (max.k - 
            1) * (1/max.k)^2
        for (ii in 1:n.test) {
            probs <- sort(prob.test[ii, ], decreasing = T)
            confidence.vector[ii] <- 1 - ((1 - probs[1])^2 + 
                sum(probs[2:max.k]^2))/random.pred.conf
            confidence.vector[ii] <- signif(confidence.vector[ii], 
                digits = 3)
            if (class.list.test[ii] == as.character(test.pred[ii])) {
                bscore[ii] <- signif((1 - probs[1])^2, digits = 2)
            }
            else {
                bscore[ii] <- signif(probs[1]^2, digits = 2)
            }
        }
        confidence.call <- ifelse(confidence.vector >= high.conf.thres, 
            " H ", " L ")
        error.call <- ifelse(class.list.test == as.character(test.pred), 
            "   ", " * ")
        no.call <- ifelse(confidence.vector >= high.conf.thres, 
            0, 1)
        real.error <- ifelse(((no.call == 0) & (error.call == 
            " * ")), 1, 0)
        correct.call <- ifelse(((no.call == 0) & (error.call == 
            "   ")), 1, 0)
        col.symbols.test <- paste(confidence.call, error.call)
        class.names <- names(data.frame(prob.test))
        Brier.test <- signif(mean(bscore), digits = 2)
        test.results <- data.frame(cbind(as.character(sample.names.test), 
            class.list.test, as.character(test.pred), error.call, 
            confidence.call, confidence.vector, no.call, real.error, 
            correct.call, prob.test, bscore))
        names(test.results)[1] <- "Test Sample Name"
        names(test.results)[2] <- "Actual"
        names(test.results)[3] <- "Predicted"
        names(test.results)[4] <- "Error (*)"
        names(test.results)[5] <- "Conf (H/L)"
        names(test.results)[6] <- "Conf"
        names(test.results)[7] <- "No Call"
        names(test.results)[8] <- "Real Error"
        names(test.results)[9] <- "Correct Call"
        names(test.results)[10 + length(class.phen.train)] <- "Brier score"
        print(c("Brier score (Test) = ", Brier.test))
        write("\n Test Results \n", file = prediction.results.file, 
            append = T)
        write.table(test.results, file = prediction.results.file, 
            append = T, quote = F, row.names = T, sep = "\t")
        write(c("\n\n Brier score (Test) = ", Brier.test), file = prediction.results.file, 
            append = T)
        no.call.list <- split(no.call, class.list.test)
        real.error.list <- split(real.error, class.list.test)
        correct.call.list <- split(correct.call, class.list.test)
        count.class <- c(sapply(no.call.list, length), length(no.call))
        no.call.class <- c(sapply(no.call.list, sum), sum(no.call))
        real.error.class <- c(sapply(real.error.list, sum), sum(real.error))
        correct.call.class <- c(sapply(correct.call.list, sum), 
            sum(correct.call))
        no.call.class.pct <- no.call.class/count.class
        real.error.class.pct <- real.error.class/count.class
        correct.call.class.pct <- correct.call.class/count.class
        perf.table.test <- data.frame(cbind(count.class, no.call.class, 
            no.call.class.pct, real.error.class, real.error.class.pct, 
            correct.call.class, correct.call.class.pct))
        names(perf.table.test) <- c("Count", "No Call", "No Call (\%)", 
            "Real Error", "Real Error (\%)", "Correct Call", "Correct Call (\%)")
        row.names(perf.table.test) <- c(names(no.call.list), 
            "Total")
        write.table(perf.table.test, file = prediction.results.file, 
            append = T, quote = F, row.names = T, sep = "\t")
        print(perf.table.test)
        conf.table.test <- table(class.list.test, test.pred)
        print(conf.table.test)
        write("\n\n Confusion Matrix (Test) \n", file = prediction.results.file, 
            append = T)
        write.table(conf.table.test, file = prediction.results.file, 
            append = T, quote = F, row.names = T, sep = "\t")
        height <- ifelse(length(class.phen.train) > 50, 20, 0.3 * 
            length(class.phen.train) + 10)
        if (non.interactive.run == F) {
            if (.Platform$OS.type == "windows") {
                plot.filename <- test.pred.plot
                x11(height = height, width = 40)
            }
            else if (.Platform$OS.type == "unix") {
                plot.filename <- paste(test.pred.plot, ".pdf", 
                  sep = "", collapse = "")
                pdf(file = plot.filename, height = height, width = 40)
            }
        }
        else {
            if (.Platform$OS.type == "unix") {
                plot.filename <- paste(test.pred.plot, ".pdf", 
                  sep = "", collapse = "")
                pdf(file = plot.filename, height = height, width = 40)
            }
            else if (.Platform$OS.type == "windows") {
                plot.filename <- paste(test.pred.plot, ".pdf", 
                  sep = "", collapse = "")
                pdf(file = plot.filename, height = height, width = 40)
            }
        }
        nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(5, 
            1), heights = 1, respect = FALSE)
        MSIG.HeatMapPlot.3(V = t(prob.test), row.names = names(test.results)[seq(10, 
            10 + length(class.phen.train) - 1)], col.labels = class.labels.test, 
            col.names = as.character(error.call), col.classes = class.names, 
            phen.cmap = col[1:length(class.names)], main = "Test Samples Predictions", 
            sub = " ", xlab = " ", ylab = " ", row.norm = F, 
            cmap.type = 2)
        leg.txt <- class.phen
        n.phen <- length(class.phen)
        p.vec <- symbs[1:n.phen]
        c.vec <- col[1:n.phen]
        par(mar = c(0, 0, 0, 0))
        plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), 
            axes = F, type = "n", xlab = "", ylab = "")
        legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
            yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
            cex = symbol.scaling * 1.25, pt.cex = symbol.scaling * 
                2.5)
        if (non.interactive.run == F) {
            if (.Platform$OS.type == "windows") {
                savePlot(filename = test.pred.plot, type = "jpeg", 
                  device = dev.cur())
            }
            else if (.Platform$OS.type == "unix") {
                dev.off()
            }
        }
        else {
            dev.off()
        }
    }
    if (Ns - length(model.set) > 0) {
        V <- cbind(t(prob.train), t(prob.test))
        W <- data.frame(V)
        names(W) <- c(as.character(sample.names.train), as.character(sample.names.test))
        row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 
            1)]
    }
    else {
        V <- t(prob.train)
        W <- data.frame(V)
        names(W) <- as.character(sample.names.train)
        row.names(W) <- names(train.results)[seq(7, 7 + length(class.phen.train) - 
            1)]
    }
    write.gct(gct.data.frame = W, descs = row.names(W), filename = prediction.matrix.file)
    print("Done predicting test data...")
    pca <- prcomp(t(m.train), retx = TRUE, center = T, scale. = T)
    S1 <- pca$x[, 1]
    S2 <- pca$x[, 2]
    S3 <- pca$x[, 3]
    X1 <- pca$rotation[, 1]
    X2 <- pca$rotation[, 2]
    X3 <- pca$rotation[, 3]
    row.mean <- apply(m.train, MARGIN = 1, FUN = mean)
    row.sd <- apply(m.train, MARGIN = 1, FUN = sd)
    c0 <- col
    c1 <- col
    color <- c1[class.labels]
    height <- 25
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- pred.2D.plot
            x11(height = height, width = 30)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 30)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 30)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(pred.2D.plot, ".pdf", sep = "", 
                collapse = "")
            pdf(file = plot.filename, height = height, width = 30)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 
        1), heights = 1, respect = FALSE)
    if (Ns - length(model.set) > 0) {
        test.scores <- predict(pca, t(m.test))
        S1 <- c(pca$x[, 1], test.scores[, 1])
        S2 <- c(pca$x[, 2], test.scores[, 2])
        S3 <- c(pca$x[, 3], test.scores[, 3])
    }
    max.S <- max(sqrt(S1 * S1 + S2 * S2))
    max.X <- max(sqrt(X1 * X1 + X2 * X2))
    X1 <- max.S * X1/max.X
    X2 <- max.S * X2/max.X
    X3 <- max.S * X3/max.X
    max.A <- max(max.S, max.X)
    num.samples <- length(S1)
    plot(S1, S2, xlim = c(-max.A, max.A), ylim = c(-max.A, max.A), 
        type = "n", main = "  ", sub = input.ds)
    for (j in 1:num.samples) {
        if (min(class.labels) == 0) {
            symb <- symbs[class.labels[j] + 1]
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            symb <- symbs[class.labels[j]]
            color.code <- c1[class.labels[j]]
        }
        points(S1[j], S2[j], pch = symb, type = "p", cex = symbol.scaling * 
            3, bg = color.code, col = "black")
    }
    for (j in 1:N) {
        x.coor <- X1[j] * 0.925
        y.coor <- X2[j] * 0.925
        arrows(0, 0, x.coor, y.coor, lwd = 1, length = 0.15, 
            angle = 20, col = "grey50")
        if (use.feature.names == FALSE) {
            leg.txt <- paste("F", j, sep = "")
        }
        else {
            leg.txt <- gs.names[j]
        }
        text(X1[j], X2[j], labels = leg.txt, adj = NULL, pos = NULL, 
            offset = 0.5, vfont = NULL, cex = 1, col = "grey50")
    }
    points.axis <- 200
    x <- vector(length = points.axis, mode = "numeric")
    pca.x <- matrix(0, nrow = 2, ncol = points.axis * points.axis)
    for (i in 1:points.axis) {
        x[i] <- -max.A + i * 2 * max.A/points.axis
    }
    for (i in 1:points.axis) {
        for (j in 1:points.axis) {
            index.point <- i + (j - 1) * points.axis
            pca.x[1, index.point] <- -max.A + i * 2 * max.A/points.axis
            pca.x[2, index.point] <- -max.A + j * 2 * max.A/points.axis
        }
    }
    grid.H <- pca$rotation[, 1:2] \%*\% pca.x
    for (i in 1:N) {
        grid.H[i, ] <- grid.H[i, ] * row.sd[i] + row.mean[i]
    }
    grid.pred <- predict(object = svm.model, newdata = t(grid.H), 
        decision.values = T, probability = T)
    prob.test <- attr(grid.pred, "probabilities")
    z <- matrix(0, nrow = points.axis, ncol = points.axis)
    z.class <- array(dim = c(length(class.phen.train), points.axis, 
        points.axis))
    max.k <- length(prob.train[1, ])
    random.pred.conf <- ((max.k - 1)/max.k)^2 + (max.k - 1) * 
        (1/max.k)^2
    for (i in 1:points.axis) {
        for (j in 1:points.axis) {
            index.point <- i + (j - 1) * points.axis
            probs <- sort(prob.test[index.point, ], decreasing = T)
            z[i, j] <- 1 - ((1 - probs[1])^2 + sum(probs[2:max.k]^2))/random.pred.conf
            for (k in 1:length(class.phen.train)) {
                if (probs[1] == prob.test[index.point, k]) {
                  z.class[k, i, j] <- z[i, j]
                }
                else {
                  z.class[k, i, j] <- 0
                }
            }
        }
    }
    library("RColorBrewer")
    if (length(levels) > 1) {
        contour(x, x, z, levels = levels, col = brewer.pal(n = 9, 
            name = "Greys")[7], add = T)
    }
    else {
        contour(x, x, z, nlevels = nlevels, col = brewer.pal(n = 9, 
            name = "Greys")[7], add = T)
    }
    contour(x, x, z, levels = 0.01, col = brewer.pal(n = 9, name = "Greys")[7], 
        lwd = 2, add = T)
    for (k in 1:length(class.phen.train)) {
        contour(x, x, z.class[k, , ], levels = high.conf.thres, 
            col = col[k], lwd = 2, add = T)
    }
    leg.txt <- class.phen
    n.phen <- length(class.phen)
    p.vec <- symbs[1:n.phen]
    c.vec <- c1[1:n.phen]
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 0), c(1, 1), xlim = c(0, 1), ylim = c(0, 1), axes = F, 
        type = "n", xlab = "", ylab = "")
    legend(x = 0, y = 1, legend = leg.txt, bty = "n", xjust = 0, 
        yjust = 1, pch = p.vec, pt.bg = c.vec, col = "black", 
        cex = symbol.scaling * 1.15, pt.cex = symbol.scaling * 
            3)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = pred.2D.plot, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
