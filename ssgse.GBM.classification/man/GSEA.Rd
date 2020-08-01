\name{GSEA}
\alias{GSEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GSEA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEA(input.ds, input.cls, gene.ann = "", gs.db, gs.ann = "", output.directory = "", doc.string = "GSEA.analysis", non.interactive.run = F, reshuffling.type = "sample.labels", nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1, fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, topgs = 10, adjust.FDR.q.val = F, gs.size.threshold.min = 25, gs.size.threshold.max = 500, reverse.sign = F, preproc.type = 0, random.seed = 123456, perm.type = 0, fraction = 1, replace = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{gene.ann}{
%%     ~~Describe \code{gene.ann} here~~
}
  \item{gs.db}{
%%     ~~Describe \code{gs.db} here~~
}
  \item{gs.ann}{
%%     ~~Describe \code{gs.ann} here~~
}
  \item{output.directory}{
%%     ~~Describe \code{output.directory} here~~
}
  \item{doc.string}{
%%     ~~Describe \code{doc.string} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{reshuffling.type}{
%%     ~~Describe \code{reshuffling.type} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{weighted.score.type}{
%%     ~~Describe \code{weighted.score.type} here~~
}
  \item{nom.p.val.threshold}{
%%     ~~Describe \code{nom.p.val.threshold} here~~
}
  \item{fwer.p.val.threshold}{
%%     ~~Describe \code{fwer.p.val.threshold} here~~
}
  \item{fdr.q.val.threshold}{
%%     ~~Describe \code{fdr.q.val.threshold} here~~
}
  \item{topgs}{
%%     ~~Describe \code{topgs} here~~
}
  \item{adjust.FDR.q.val}{
%%     ~~Describe \code{adjust.FDR.q.val} here~~
}
  \item{gs.size.threshold.min}{
%%     ~~Describe \code{gs.size.threshold.min} here~~
}
  \item{gs.size.threshold.max}{
%%     ~~Describe \code{gs.size.threshold.max} here~~
}
  \item{reverse.sign}{
%%     ~~Describe \code{reverse.sign} here~~
}
  \item{preproc.type}{
%%     ~~Describe \code{preproc.type} here~~
}
  \item{random.seed}{
%%     ~~Describe \code{random.seed} here~~
}
  \item{perm.type}{
%%     ~~Describe \code{perm.type} here~~
}
  \item{fraction}{
%%     ~~Describe \code{fraction} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
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
function (input.ds, input.cls, gene.ann = "", gs.db, gs.ann = "", 
    output.directory = "", doc.string = "GSEA.analysis", non.interactive.run = F, 
    reshuffling.type = "sample.labels", nperm = 1000, weighted.score.type = 1, 
    nom.p.val.threshold = -1, fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, 
    topgs = 10, adjust.FDR.q.val = F, gs.size.threshold.min = 25, 
    gs.size.threshold.max = 500, reverse.sign = F, preproc.type = 0, 
    random.seed = 123456, perm.type = 0, fraction = 1, replace = F) 
{
    print(" *** Running GSEA Analysis...")
    if (output.directory != "") {
        filename <- paste(output.directory, doc.string, "_params.txt", 
            sep = "", collapse = "")
        time.string <- as.character(as.POSIXlt(Sys.time(), "GMT"))
        write(paste("Run of GSEA on ", time.string), file = filename)
        if (is.data.frame(input.ds)) {
        }
        else {
            write(paste("input.ds=", input.ds, sep = " "), file = filename, 
                append = T)
        }
        if (is.list(input.cls)) {
        }
        else {
            write(paste("input.cls=", input.cls, sep = " "), 
                file = filename, append = T)
        }
        if (is.data.frame(gene.ann)) {
        }
        else {
            write(paste("gene.ann =", gene.ann, sep = " "), file = filename, 
                append = T)
        }
        if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
        }
        else {
            write(paste("gs.db=", gs.db, sep = " "), file = filename, 
                append = T)
        }
        if (is.data.frame(gs.ann)) {
        }
        else {
            write(paste("gs.ann =", gs.ann, sep = " "), file = filename, 
                append = T)
        }
        write(paste("output.directory =", output.directory, sep = " "), 
            file = filename, append = T)
        write(paste("doc.string = ", doc.string, sep = " "), 
            file = filename, append = T)
        write(paste("non.interactive.run =", non.interactive.run, 
            sep = " "), file = filename, append = T)
        write(paste("reshuffling.type =", reshuffling.type, sep = " "), 
            file = filename, append = T)
        write(paste("nperm =", nperm, sep = " "), file = filename, 
            append = T)
        write(paste("weighted.score.type =", weighted.score.type, 
            sep = " "), file = filename, append = T)
        write(paste("nom.p.val.threshold =", nom.p.val.threshold, 
            sep = " "), file = filename, append = T)
        write(paste("fwer.p.val.threshold =", fwer.p.val.threshold, 
            sep = " "), file = filename, append = T)
        write(paste("fdr.q.val.threshold =", fdr.q.val.threshold, 
            sep = " "), file = filename, append = T)
        write(paste("topgs =", topgs, sep = " "), file = filename, 
            append = T)
        write(paste("adjust.FDR.q.val =", adjust.FDR.q.val, sep = " "), 
            file = filename, append = T)
        write(paste("gs.size.threshold.min =", gs.size.threshold.min, 
            sep = " "), file = filename, append = T)
        write(paste("gs.size.threshold.max =", gs.size.threshold.max, 
            sep = " "), file = filename, append = T)
        write(paste("reverse.sign =", reverse.sign, sep = " "), 
            file = filename, append = T)
        write(paste("preproc.type =", preproc.type, sep = " "), 
            file = filename, append = T)
        write(paste("random.seed =", random.seed, sep = " "), 
            file = filename, append = T)
        write(paste("perm.type =", perm.type, sep = " "), file = filename, 
            append = T)
        write(paste("fraction =", fraction, sep = " "), file = filename, 
            append = T)
        write(paste("replace =", replace, sep = " "), file = filename, 
            append = T)
    }
    if (.Platform$OS.type == "windows") {
        memory.limit(6e+09)
        memory.limit()
    }
    set.seed(seed = random.seed, kind = NULL)
    adjust.param <- 0.5
    gc()
    time1 <- proc.time()
    if (is.data.frame(input.ds)) {
        dataset <- input.ds
    }
    else {
        if (regexpr(pattern = ".gct", input.ds) == -1) {
            dataset <- GSEA.Res2Frame(filename = input.ds)
        }
        else {
            dataset <- GSEA.Gct2Frame(filename = input.ds)
        }
    }
    gene.labels <- row.names(dataset)
    sample.names <- names(dataset)
    A <- data.matrix(dataset)
    dim(A)
    cols <- length(A[1, ])
    rows <- length(A[, 1])
    if (preproc.type == 1) {
        A <- GSEA.NormalizeCols(A)
    }
    else if (preproc.type == 2) {
        for (j in 1:cols) {
            A[, j] <- rank(A[, j])
        }
        A <- GSEA.NormalizeRows(A)
    }
    else if (preproc.type == 3) {
        for (j in 1:cols) {
            A[, j] <- rank(A[, j])
        }
    }
    if (is.list(input.cls)) {
        CLS <- input.cls
    }
    else {
        CLS <- GSEA.ReadClsFile(file = input.cls)
    }
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    if (reverse.sign == T) {
        phen1 <- class.phen[2]
        phen2 <- class.phen[1]
    }
    else {
        phen1 <- class.phen[1]
        phen2 <- class.phen[2]
    }
    col.index <- order(class.labels, decreasing = F)
    class.labels <- class.labels[col.index]
    sample.names <- sample.names[col.index]
    for (j in 1:rows) {
        A[j, ] <- A[j, col.index]
    }
    names(A) <- sample.names
    if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
        temp <- gs.db
    }
    else {
        temp <- readLines(gs.db)
    }
    max.Ng <- length(temp)
    temp.size.G <- vector(length = max.Ng, mode = "numeric")
    for (i in 1:max.Ng) {
        temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 
            2
    }
    max.size.G <- max(temp.size.G)
    gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, 
        ncol = max.size.G)
    temp.names <- vector(length = max.Ng, mode = "character")
    temp.desc <- vector(length = max.Ng, mode = "character")
    gs.count <- 1
    for (i in 1:max.Ng) {
        gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 
            2
        gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
        gene.set.name <- gs.line[1]
        gene.set.desc <- gs.line[1]
        gene.set.tags <- vector(length = gene.set.size, mode = "character")
        for (j in 1:gene.set.size) {
            gene.set.tags[j] <- gs.line[j + 2]
        }
        existing.set <- is.element(gene.set.tags, gene.labels)
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < gs.size.threshold.min) || (set.size > 
            gs.size.threshold.max)) 
            next
        temp.size.G[gs.count] <- set.size
        gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", 
            max.size.G - temp.size.G[gs.count]))
        temp.names[gs.count] <- gene.set.name
        temp.desc[gs.count] <- gene.set.desc
        gs.count <- gs.count + 1
    }
    Ng <- gs.count - 1
    gs.names <- vector(length = Ng, mode = "character")
    gs.desc <- vector(length = Ng, mode = "character")
    size.G <- vector(length = Ng, mode = "numeric")
    gs.names <- temp.names[1:Ng]
    gs.desc <- temp.desc[1:Ng]
    size.G <- temp.size.G[1:Ng]
    N <- length(A[, 1])
    Ns <- length(A[1, ])
    print(c("Number of genes:", N))
    print(c("Number of Gene Sets:", Ng))
    print(c("Number of samples:", Ns))
    print(c("Original number of Gene Sets:", max.Ng))
    print(c("Maximum gene set size:", max.size.G))
    all.gene.descs <- vector(length = N, mode = "character")
    all.gene.symbols <- vector(length = N, mode = "character")
    all.gs.descs <- vector(length = Ng, mode = "character")
    if (is.data.frame(gene.ann)) {
        temp <- gene.ann
        a.size <- length(temp[, 1])
        print(c("Number of gene annotation file entries:", a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gene.labels, accs)
        all.gene.descs <- as.character(temp[locs, "Gene.Title"])
        all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
        rm(temp)
    }
    else if (gene.ann == "") {
        for (i in 1:N) {
            all.gene.descs[i] <- " -- "
            all.gene.symbols[i] <- " -- "
        }
    }
    else {
        temp <- read.delim(gene.ann, header = T, sep = ",", comment.char = "", 
            as.is = T)
        a.size <- length(temp[, 1])
        print(c("Number of gene annotation file entries:", a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gene.labels, accs)
        all.gene.descs <- as.character(temp[locs, "Gene.Title"])
        all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
        rm(temp)
    }
    if (is.data.frame(gs.ann)) {
        temp <- gs.ann
        a.size <- length(temp[, 1])
        print(c("Number of gene set annotation file entries:", 
            a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gs.names, accs)
        all.gs.descs <- as.character(temp[locs, "SOURCE"])
        rm(temp)
    }
    else if (gs.ann == "") {
        for (i in 1:N) {
            all.gene.descs[i] <- " -- "
            all.gene.symbols[i] <- " -- "
        }
    }
    else {
        temp <- read.delim(gs.ann, header = T, sep = "\t", comment.char = "", 
            as.is = T)
        a.size <- length(temp[, 1])
        print(c("Number of gene set annotation file entries:", 
            a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gs.names, accs)
        all.gs.descs <- as.character(temp[locs, "SOURCE"])
        rm(temp)
    }
    Obs.indicator <- matrix(nrow = Ng, ncol = N)
    Obs.RES <- matrix(nrow = Ng, ncol = N)
    Obs.ES <- vector(length = Ng, mode = "numeric")
    Obs.arg.ES <- vector(length = Ng, mode = "numeric")
    Obs.ES.norm <- vector(length = Ng, mode = "numeric")
    time2 <- proc.time()
    obs.s2n <- vector(length = N, mode = "numeric")
    signal.strength <- vector(length = Ng, mode = "numeric")
    tag.frac <- vector(length = Ng, mode = "numeric")
    gene.frac <- vector(length = Ng, mode = "numeric")
    coherence.ratio <- vector(length = Ng, mode = "numeric")
    obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
    correl.matrix <- matrix(nrow = N, ncol = nperm)
    obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
    order.matrix <- matrix(nrow = N, ncol = nperm)
    obs.order.matrix <- matrix(nrow = N, ncol = nperm)
    nperm.per.call <- 200
    n.groups <- nperm\%/\%nperm.per.call
    n.rem <- nperm\%\%nperm.per.call
    n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
    n.ends <- cumsum(n.perms)
    n.starts <- n.ends - n.perms + 1
    if (n.rem == 0) {
        n.tot <- n.groups
    }
    else {
        n.tot <- n.groups + 1
    }
    for (nk in 1:n.tot) {
        call.nperm <- n.perms[nk]
        print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", 
            n.starts[nk], "--", n.ends[nk], sep = " "))
        O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, 
            permutation.type = perm.type, sigma.correction = "GeneCluster", 
            fraction = fraction, replace = replace, reverse.sign = reverse.sign)
        gc()
        order.matrix[, n.starts[nk]:n.ends[nk]] <- O$order.matrix
        obs.order.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
        correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
        obs.correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
        rm(O)
    }
    obs.s2n <- apply(obs.correl.matrix, 1, median)
    obs.index <- order(obs.s2n, decreasing = T)
    obs.s2n <- sort(obs.s2n, decreasing = T)
    obs.gene.labels <- gene.labels[obs.index]
    obs.gene.descs <- all.gene.descs[obs.index]
    obs.gene.symbols <- all.gene.symbols[obs.index]
    for (r in 1:nperm) {
        correl.matrix[, r] <- correl.matrix[order.matrix[, r], 
            r]
    }
    for (r in 1:nperm) {
        obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[, 
            r], r]
    }
    gene.list2 <- obs.index
    for (i in 1:Ng) {
        print(paste("Computing observed enrichment for gene set:", 
            i, gs.names[i], sep = " "))
        gene.set <- gs[i, gs[i, ] != "null"]
        gene.set2 <- vector(length = length(gene.set), mode = "numeric")
        gene.set2 <- match(gene.set, gene.labels)
        GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, 
            gene.set = gene.set2, weighted.score.type = weighted.score.type, 
            correl.vector = obs.s2n)
        Obs.ES[i] <- GSEA.results$ES
        Obs.arg.ES[i] <- GSEA.results$arg.ES
        Obs.RES[i, ] <- GSEA.results$RES
        Obs.indicator[i, ] <- GSEA.results$indicator
        if (Obs.ES[i] >= 0) {
            tag.frac[i] <- sum(Obs.indicator[i, 1:Obs.arg.ES[i]])/size.G[i]
            gene.frac[i] <- Obs.arg.ES[i]/N
        }
        else {
            tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
            gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
        }
        signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * 
            (N/(N - size.G[i]))
    }
    phi <- matrix(nrow = Ng, ncol = nperm)
    phi.norm <- matrix(nrow = Ng, ncol = nperm)
    obs.phi <- matrix(nrow = Ng, ncol = nperm)
    if (reshuffling.type == "sample.labels") {
        for (i in 1:Ng) {
            print(paste("Computing random permutations' enrichment for gene set:", 
                i, gs.names[i], sep = " "))
            gene.set <- gs[i, gs[i, ] != "null"]
            gene.set2 <- vector(length = length(gene.set), mode = "numeric")
            gene.set2 <- match(gene.set, gene.labels)
            for (r in 1:nperm) {
                gene.list2 <- order.matrix[, r]
                GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, 
                  gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                  correl.vector = correl.matrix[, r])
                phi[i, r] <- GSEA.results$ES
                obs.gene.list2 <- obs.order.matrix[, r]
                GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
                  gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                  correl.vector = obs.correl.matrix[, r])
                obs.phi[i, r] <- GSEA.results$ES
            }
            gc()
        }
    }
    else if (reshuffling.type == "gene.labels") {
        for (i in 1:Ng) {
            gene.set <- gs[i, gs[i, ] != "null"]
            gene.set2 <- vector(length = length(gene.set), mode = "numeric")
            gene.set2 <- match(gene.set, gene.labels)
            for (r in 1:nperm) {
                reshuffled.gene.labels <- sample(1:rows)
                GSEA.results <- GSEA.EnrichmentScore(gene.list = reshuffled.gene.labels, 
                  gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                  correl.vector = obs.s2n)
                phi[i, r] <- GSEA.results$ES
                obs.gene.list2 <- obs.order.matrix[, r]
                GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
                  gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                  correl.vector = obs.correl.matrix[, r])
                obs.phi[i, r] <- GSEA.results$ES
            }
            gc()
        }
    }
    print("Computing nominal p-values...")
    p.vals <- matrix(0, nrow = Ng, ncol = 2)
    for (i in 1:Ng) {
        pos.phi <- phi[i, phi[i, ] >= 0]
        neg.phi <- phi[i, phi[i, ] < 0]
        if (Obs.ES[i] >= 0) {
            p.vals[i, 1] <- sum(pos.phi >= Obs.ES[i])/length(pos.phi)
            p.vals[i, 1] <- signif(p.vals[i, 1], digits = 5)
        }
        else {
            p.vals[i, 1] <- sum(neg.phi < Obs.ES[i])/length(neg.phi)
            p.vals[i, 1] <- signif(p.vals[i, 1], digits = 5)
        }
    }
    erf <- function(x) {
        2 * pnorm(sqrt(2) * x)
    }
    KS.mean <- function(N) {
        S <- 0
        for (k in -100:100) {
            if (k == 0) 
                next
            S <- S + 4 * (-1)^(k + 1) * (0.25 * exp(-2 * k * 
                k * N) - sqrt(2 * pi) * erf(sqrt(2 * N) * k)/(16 * 
                k * sqrt(N)))
        }
        return(abs(S))
    }
    KS.mean.table <- vector(length = 5000, mode = "numeric")
    for (i in 1:5000) {
        KS.mean.table[i] <- KS.mean(i)
    }
    KS.size <- vector(length = Ng, mode = "numeric")
    print("Computing rescaling normalization for each gene set null...")
    for (i in 1:Ng) {
        pos.phi <- phi[i, phi[i, ] >= 0]
        neg.phi <- phi[i, phi[i, ] < 0]
        pos.m <- mean(pos.phi)
        neg.m <- mean(abs(neg.phi))
        if (Obs.ES[i] >= 0) {
            KS.size[i] <- which.min(abs(KS.mean.table - pos.m))
        }
        else {
            KS.size[i] <- which.min(abs(KS.mean.table - neg.m))
        }
        pos.phi <- pos.phi/pos.m
        neg.phi <- neg.phi/neg.m
        for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                phi.norm[i, j] <- phi[i, j]/pos.m
            }
            else {
                phi.norm[i, j] <- phi[i, j]/neg.m
            }
        }
        for (j in 1:nperm) {
            if (obs.phi[i, j] >= 0) {
                obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
            }
            else {
                obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
            }
        }
        if (Obs.ES[i] >= 0) {
            Obs.ES.norm[i] <- Obs.ES[i]/pos.m
        }
        else {
            Obs.ES.norm[i] <- Obs.ES[i]/neg.m
        }
    }
    print("Computing FWER p-values...")
    max.ES.vals.p <- NULL
    max.ES.vals.n <- NULL
    for (j in 1:nperm) {
        pos.phi <- phi.norm[phi.norm[, j] >= 0, j]
        neg.phi <- phi.norm[phi.norm[, j] < 0, j]
        if (length(pos.phi) > 0) {
            max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
        }
        if (length(neg.phi) > 0) {
            max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
        }
    }
    for (i in 1:Ng) {
        if (Obs.ES.norm[i] >= 0) {
            p.vals[i, 2] <- sum(max.ES.vals.p >= Obs.ES.norm[i])/length(max.ES.vals.p)
        }
        else {
            p.vals[i, 2] <- sum(max.ES.vals.n < Obs.ES.norm[i])/length(max.ES.vals.n)
        }
        p.vals[i, 2] <- signif(p.vals[i, 2], digits = 4)
    }
    print("Computing FDR q-values...")
    NES <- vector(length = Ng, mode = "numeric")
    phi.norm.mean <- vector(length = Ng, mode = "numeric")
    obs.phi.norm.mean <- vector(length = Ng, mode = "numeric")
    phi.norm.median <- vector(length = Ng, mode = "numeric")
    obs.phi.norm.median <- vector(length = Ng, mode = "numeric")
    phi.norm.mean <- vector(length = Ng, mode = "numeric")
    obs.phi.mean <- vector(length = Ng, mode = "numeric")
    FDR.mean <- vector(length = Ng, mode = "numeric")
    FDR.median <- vector(length = Ng, mode = "numeric")
    phi.norm.median.d <- vector(length = Ng, mode = "numeric")
    obs.phi.norm.median.d <- vector(length = Ng, mode = "numeric")
    Obs.ES.index <- order(Obs.ES.norm, decreasing = T)
    Orig.index <- seq(1, Ng)
    Orig.index <- Orig.index[Obs.ES.index]
    Orig.index <- order(Orig.index, decreasing = F)
    Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
    gs.names.sorted <- gs.names[Obs.ES.index]
    for (k in 1:Ng) {
        NES[k] <- Obs.ES.norm.sorted[k]
        count.col <- vector(length = nperm, mode = "numeric")
        obs.count.col <- vector(length = nperm, mode = "numeric")
        for (i in 1:nperm) {
            if (NES[k] >= 0) {
                count.col.norm <- sum(phi.norm[, i] >= 0)
                obs.count.col.norm <- sum(obs.phi.norm[, i] >= 
                  0)
                count.col[i] <- sum(phi.norm[, i] >= NES[k])/count.col.norm
                obs.count.col[i] <- sum(obs.phi.norm[, i] >= 
                  NES[k])/obs.count.col.norm
            }
            else {
                count.col.norm <- sum(phi.norm[, i] < 0)
                obs.count.col.norm <- sum(obs.phi.norm[, i] < 
                  0)
                count.col[i] <- sum(phi.norm[, i] <= NES[k])/count.col.norm
                obs.count.col[i] <- sum(obs.phi.norm[, i] <= 
                  NES[k])/obs.count.col.norm
            }
        }
        phi.norm.mean[k] <- mean(count.col)
        obs.phi.norm.mean[k] <- mean(obs.count.col)
        phi.norm.median[k] <- median(count.col)
        obs.phi.norm.median[k] <- median(obs.count.col)
        FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 
            1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
        FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 
            1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
    }
    if (adjust.FDR.q.val == T) {
        pos.nes <- length(NES[NES >= 0])
        min.FDR.mean <- FDR.mean[pos.nes]
        min.FDR.median <- FDR.median[pos.nes]
        for (k in seq(pos.nes - 1, 1, -1)) {
            if (FDR.mean[k] < min.FDR.mean) {
                min.FDR.mean <- FDR.mean[k]
            }
            if (min.FDR.mean < FDR.mean[k]) {
                FDR.mean[k] <- min.FDR.mean
            }
        }
        neg.nes <- pos.nes + 1
        min.FDR.mean <- FDR.mean[neg.nes]
        min.FDR.median <- FDR.median[neg.nes]
        for (k in seq(neg.nes + 1, Ng)) {
            if (FDR.mean[k] < min.FDR.mean) {
                min.FDR.mean <- FDR.mean[k]
            }
            if (min.FDR.mean < FDR.mean[k]) {
                FDR.mean[k] <- min.FDR.mean
            }
        }
    }
    obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
    phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
    FDR.mean.sorted <- FDR.mean[Orig.index]
    FDR.median.sorted <- FDR.median[Orig.index]
    glob.p.vals <- vector(length = Ng, mode = "numeric")
    NULL.pass <- vector(length = nperm, mode = "numeric")
    OBS.pass <- vector(length = nperm, mode = "numeric")
    for (k in 1:Ng) {
        NES[k] <- Obs.ES.norm.sorted[k]
        if (NES[k] >= 0) {
            for (i in 1:nperm) {
                NULL.pass[i] <- sum(phi.norm[, i] >= NES[k])/sum(phi.norm[, 
                  i] >= 0)
                OBS.pass[i] <- sum(obs.phi.norm[, i] >= NES[k])/sum(obs.phi.norm[, 
                  i] >= 0)
            }
        }
        else {
            for (i in 1:nperm) {
                NULL.pass[i] <- sum(phi.norm[, i] <= NES[k])/sum(phi.norm[, 
                  i] < 0)
                OBS.pass[i] <- sum(obs.phi.norm[, i] <= NES[k])/sum(obs.phi.norm[, 
                  i] < 0)
            }
        }
        glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
    }
    glob.p.vals.sorted <- glob.p.vals[Orig.index]
    print("Producing result tables and plots...")
    Obs.ES <- signif(Obs.ES, digits = 5)
    Obs.ES.norm <- signif(Obs.ES.norm, digits = 5)
    p.vals <- signif(p.vals, digits = 4)
    signal.strength <- signif(signal.strength, digits = 3)
    tag.frac <- signif(tag.frac, digits = 3)
    gene.frac <- signif(gene.frac, digits = 3)
    FDR.mean.sorted <- signif(FDR.mean.sorted, digits = 5)
    FDR.median.sorted <- signif(FDR.median.sorted, digits = 5)
    glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits = 5)
    report <- data.frame(cbind(gs.names, size.G, all.gs.descs, 
        Obs.ES, Obs.ES.norm, p.vals[, 1], FDR.mean.sorted, p.vals[, 
            2], tag.frac, gene.frac, signal.strength, FDR.median.sorted, 
        glob.p.vals.sorted))
    names(report) <- c("GS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", 
        "FDR q-val", "FWER p-val", "Tag \%", "Gene \%", "Signal", 
        "FDR (median)", "glob.p.val")
    report2 <- report
    report.index2 <- order(Obs.ES.norm, decreasing = T)
    for (i in 1:Ng) {
        report2[i, ] <- report[report.index2[i], ]
    }
    report3 <- report
    report.index3 <- order(Obs.ES.norm, decreasing = F)
    for (i in 1:Ng) {
        report3[i, ] <- report[report.index3[i], ]
    }
    phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
    phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
    report.phen1 <- report2[1:phen1.rows, ]
    report.phen2 <- report3[1:phen2.rows, ]
    if (output.directory != "") {
        if (phen1.rows > 0) {
            filename <- paste(output.directory, doc.string, ".results.report.", 
                phen1, ".txt", sep = "", collapse = "")
            write.table(report.phen1, file = filename, quote = F, 
                row.names = F, sep = "\t")
        }
        if (phen2.rows > 0) {
            filename <- paste(output.directory, doc.string, ".results.report.", 
                phen2, ".txt", sep = "", collapse = "")
            write.table(report.phen2, file = filename, quote = F, 
                row.names = F, sep = "\t")
        }
    }
    if (output.directory != "") {
        if (non.interactive.run == F) {
            if (.Platform$OS.type == "windows") {
                glob.filename <- paste(output.directory, doc.string, 
                  ".global.plots", sep = "", collapse = "")
                x11(width = 10, height = 10)
            }
            else if (.Platform$OS.type == "unix") {
                glob.filename <- paste(output.directory, doc.string, 
                  ".global.plots.pdf", sep = "", collapse = "")
                pdf(file = glob.filename, height = 10, width = 10)
            }
        }
        else {
            if (.Platform$OS.type == "unix") {
                glob.filename <- paste(output.directory, doc.string, 
                  ".global.plots.pdf", sep = "", collapse = "")
                pdf(file = glob.filename, height = 10, width = 10)
            }
            else if (.Platform$OS.type == "windows") {
                glob.filename <- paste(output.directory, doc.string, 
                  ".global.plots.pdf", sep = "", collapse = "")
                pdf(file = glob.filename, height = 10, width = 10)
            }
        }
    }
    nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = T), c(1, 
        1), c(1, 1), TRUE)
    location <- 1:N
    max.corr <- max(obs.s2n)
    min.corr <- min(obs.s2n)
    x <- plot(location, obs.s2n, ylab = "Signal to Noise Ratio (S2N)", 
        xlab = "Gene List Location", main = "Gene List Correlation (S2N) Profile", 
        type = "l", lwd = 2, cex = 0.9, col = 1)
    for (i in seq(1, N, 20)) {
        lines(c(i, i), c(0, obs.s2n[i]), lwd = 3, cex = 0.9, 
            col = colors()[12])
    }
    x <- points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, 
        col = 1)
    lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1)
    temp <- order(abs(obs.s2n), decreasing = T)
    arg.correl <- temp[N]
    lines(c(arg.correl, arg.correl), c(min.corr, 0.7 * max.corr), 
        lwd = 2, lty = 3, cex = 0.9, col = 1)
    area.bias <- signif(100 * (sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N]))/sum(abs(obs.s2n[1:N])), 
        digits = 3)
    area.phen <- ifelse(area.bias >= 0, phen1, phen2)
    delta.string <- paste("Corr. Area Bias to \"", area.phen, 
        "\" =", abs(area.bias), "\%", sep = "", collapse = "")
    zero.crossing.string <- paste("Zero Crossing at location ", 
        arg.correl, " (", signif(100 * arg.correl/N, digits = 3), 
        " \%)")
    leg.txt <- c(delta.string, zero.crossing.string)
    legend(x = N/10, y = max.corr, bty = "n", bg = "white", legend = leg.txt, 
        cex = 0.9)
    leg.txt <- paste("\"", phen1, "\" ", sep = "", collapse = "")
    text(x = 1, y = -0.05 * max.corr, adj = c(0, 1), labels = leg.txt, 
        cex = 0.9)
    leg.txt <- paste("\"", phen2, "\" ", sep = "", collapse = "")
    text(x = N, y = 0.05 * max.corr, adj = c(1, 0), labels = leg.txt, 
        cex = 0.9)
    if (Ng > 1) {
        phi.densities.pos <- matrix(0, nrow = 512, ncol = nperm)
        phi.densities.neg <- matrix(0, nrow = 512, ncol = nperm)
        obs.phi.densities.pos <- matrix(0, nrow = 512, ncol = nperm)
        obs.phi.densities.neg <- matrix(0, nrow = 512, ncol = nperm)
        phi.density.mean.pos <- vector(length = 512, mode = "numeric")
        phi.density.mean.neg <- vector(length = 512, mode = "numeric")
        obs.phi.density.mean.pos <- vector(length = 512, mode = "numeric")
        obs.phi.density.mean.neg <- vector(length = 512, mode = "numeric")
        phi.density.median.pos <- vector(length = 512, mode = "numeric")
        phi.density.median.neg <- vector(length = 512, mode = "numeric")
        obs.phi.density.median.pos <- vector(length = 512, mode = "numeric")
        obs.phi.density.median.neg <- vector(length = 512, mode = "numeric")
        x.coor.pos <- vector(length = 512, mode = "numeric")
        x.coor.neg <- vector(length = 512, mode = "numeric")
        for (i in 1:nperm) {
            pos.phi <- phi.norm[phi.norm[, i] >= 0, i]
            if (length(pos.phi) > 2) {
                temp <- density(pos.phi, adjust = adjust.param, 
                  n = 512, from = 0, to = 3.5)
            }
            else {
                temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, 
                  y = rep(0.001, 512))
            }
            phi.densities.pos[, i] <- temp$y
            norm.factor <- sum(phi.densities.pos[, i])
            phi.densities.pos[, i] <- phi.densities.pos[, i]/norm.factor
            if (i == 1) {
                x.coor.pos <- temp$x
            }
            neg.phi <- phi.norm[phi.norm[, i] < 0, i]
            if (length(neg.phi) > 2) {
                temp <- density(neg.phi, adjust = adjust.param, 
                  n = 512, from = -3.5, to = 0)
            }
            else {
                temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, 
                  y = rep(0.001, 512))
            }
            phi.densities.neg[, i] <- temp$y
            norm.factor <- sum(phi.densities.neg[, i])
            phi.densities.neg[, i] <- phi.densities.neg[, i]/norm.factor
            if (i == 1) {
                x.coor.neg <- temp$x
            }
            pos.phi <- obs.phi.norm[obs.phi.norm[, i] >= 0, i]
            if (length(pos.phi) > 2) {
                temp <- density(pos.phi, adjust = adjust.param, 
                  n = 512, from = 0, to = 3.5)
            }
            else {
                temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, 
                  y = rep(0.001, 512))
            }
            obs.phi.densities.pos[, i] <- temp$y
            norm.factor <- sum(obs.phi.densities.pos[, i])
            obs.phi.densities.pos[, i] <- obs.phi.densities.pos[, 
                i]/norm.factor
            neg.phi <- obs.phi.norm[obs.phi.norm[, i] < 0, i]
            if (length(neg.phi) > 2) {
                temp <- density(neg.phi, adjust = adjust.param, 
                  n = 512, from = -3.5, to = 0)
            }
            else {
                temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, 
                  y = rep(0.001, 512))
            }
            obs.phi.densities.neg[, i] <- temp$y
            norm.factor <- sum(obs.phi.densities.neg[, i])
            obs.phi.densities.neg[, i] <- obs.phi.densities.neg[, 
                i]/norm.factor
        }
        phi.density.mean.pos <- apply(phi.densities.pos, 1, mean)
        phi.density.mean.neg <- apply(phi.densities.neg, 1, mean)
        obs.phi.density.mean.pos <- apply(obs.phi.densities.pos, 
            1, mean)
        obs.phi.density.mean.neg <- apply(obs.phi.densities.neg, 
            1, mean)
        phi.density.median.pos <- apply(phi.densities.pos, 1, 
            median)
        phi.density.median.neg <- apply(phi.densities.neg, 1, 
            median)
        obs.phi.density.median.pos <- apply(obs.phi.densities.pos, 
            1, median)
        obs.phi.density.median.neg <- apply(obs.phi.densities.neg, 
            1, median)
        x <- c(x.coor.neg, x.coor.pos)
        x.plot.range <- range(x)
        y1 <- c(phi.density.mean.neg, phi.density.mean.pos)
        y2 <- c(obs.phi.density.mean.neg, obs.phi.density.mean.pos)
        y.plot.range <- c(-0.3 * max(c(y1, y2)), max(c(y1, y2)))
        plot(x, y1, xlim = x.plot.range, ylim = 1.5 * y.plot.range, 
            type = "l", lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", 
            main = "Global Observed and Null Densities (Area Normalized)")
        y1.point <- y1[seq(1, length(x), 2)]
        y2.point <- y2[seq(2, length(x), 2)]
        x1.point <- x[seq(1, length(x), 2)]
        x2.point <- x[seq(2, length(x), 2)]
        points(x, y1, type = "l", lwd = 2, col = colors()[555])
        points(x, y2, type = "l", lwd = 2, col = colors()[29])
        for (i in 1:Ng) {
            col <- ifelse(Obs.ES.norm[i] > 0, 2, 3)
            lines(c(Obs.ES.norm[i], Obs.ES.norm[i]), c(-0.2 * 
                max(c(y1, y2)), 0), lwd = 1, lty = 1, col = 1)
        }
        leg.txt <- paste("Neg. ES: \"", phen2, " \" ", sep = "", 
            collapse = "")
        text(x = x.plot.range[1], y = -0.25 * max(c(y1, y2)), 
            adj = c(0, 1), labels = leg.txt, cex = 0.9)
        leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep = "", 
            collapse = "")
        text(x = x.plot.range[2], y = -0.25 * max(c(y1, y2)), 
            adj = c(1, 1), labels = leg.txt, cex = 0.9)
        leg.txt <- c("Null Density", "Observed Density", "Observed NES values")
        c.vec <- c(colors()[555], colors()[29], 1)
        lty.vec <- c(1, 1, 1)
        lwd.vec <- c(2, 2, 2)
        legend(x = 0, y = 1.5 * y.plot.range[2], bty = "n", bg = "white", 
            legend = leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, 
            cex = 0.9)
        B <- A[obs.index, ]
        GSEA.HeatMapPlot(V = B, col.labels = class.labels, col.classes = class.phen, 
            main = "Heat Map for Genes in Dataset")
        nom.p.vals <- p.vals[Obs.ES.index, 1]
        FWER.p.vals <- p.vals[Obs.ES.index, 2]
        plot.range <- 1.25 * range(NES)
        plot(NES, FDR.mean, ylim = c(0, 1), xlim = plot.range, 
            col = 1, bg = 1, type = "p", pch = 22, cex = 0.75, 
            xlab = "NES", main = "p-values vs. NES", ylab = "p-val/q-val")
        points(NES, nom.p.vals, type = "p", col = 2, bg = 2, 
            pch = 22, cex = 0.75)
        points(NES, FWER.p.vals, type = "p", col = colors()[577], 
            bg = colors()[577], pch = 22, cex = 0.75)
        leg.txt <- c("Nominal p-value", "FWER p-value", "FDR q-value")
        c.vec <- c(2, colors()[577], 1)
        pch.vec <- c(22, 22, 22)
        legend(x = -0.5, y = 0.5, bty = "n", bg = "white", legend = leg.txt, 
            pch = pch.vec, col = c.vec, pt.bg = c.vec, cex = 0.9)
        lines(c(min(NES), max(NES)), c(nom.p.val.threshold, nom.p.val.threshold), 
            lwd = 1, lty = 2, col = 2)
        lines(c(min(NES), max(NES)), c(fwer.p.val.threshold, 
            fwer.p.val.threshold), lwd = 1, lty = 2, col = colors()[577])
        lines(c(min(NES), max(NES)), c(fdr.q.val.threshold, fdr.q.val.threshold), 
            lwd = 1, lty = 2, col = 1)
        if (non.interactive.run == F) {
            if (.Platform$OS.type == "windows") {
                savePlot(filename = glob.filename, type = "jpeg", 
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
    if (topgs > floor(Ng/2)) {
        topgs <- floor(Ng/2)
    }
    for (i in 1:Ng) {
        if ((p.vals[i, 1] <= nom.p.val.threshold) || (p.vals[i, 
            2] <= fwer.p.val.threshold) || (FDR.mean.sorted[i] <= 
            fdr.q.val.threshold) || (is.element(i, c(Obs.ES.index[1:topgs], 
            Obs.ES.index[(Ng - topgs + 1):Ng])))) {
            kk <- 1
            gene.number <- vector(length = size.G[i], mode = "character")
            gene.names <- vector(length = size.G[i], mode = "character")
            gene.symbols <- vector(length = size.G[i], mode = "character")
            gene.descs <- vector(length = size.G[i], mode = "character")
            gene.list.loc <- vector(length = size.G[i], mode = "numeric")
            core.enrichment <- vector(length = size.G[i], mode = "character")
            gene.s2n <- vector(length = size.G[i], mode = "numeric")
            gene.RES <- vector(length = size.G[i], mode = "numeric")
            rank.list <- seq(1, N)
            if (Obs.ES[i] >= 0) {
                set.k <- seq(1, N, 1)
            }
            else {
                set.k <- seq(N, 1, -1)
            }
            for (k in set.k) {
                if (Obs.indicator[i, k] == 1) {
                  gene.number[kk] <- kk
                  gene.names[kk] <- obs.gene.labels[k]
                  gene.symbols[kk] <- substr(obs.gene.symbols[k], 
                    1, 15)
                  gene.descs[kk] <- substr(obs.gene.descs[k], 
                    1, 40)
                  gene.list.loc[kk] <- k
                  gene.s2n[kk] <- signif(obs.s2n[k], digits = 3)
                  gene.RES[kk] <- signif(Obs.RES[i, k], digits = 3)
                  if (Obs.ES[i] >= 0) {
                    core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= 
                      Obs.arg.ES[i], "YES", "NO")
                  }
                  else {
                    core.enrichment[kk] <- ifelse(gene.list.loc[kk] > 
                      Obs.arg.ES[i], "YES", "NO")
                  }
                  kk <- kk + 1
                }
            }
            gene.report <- data.frame(cbind(gene.number, gene.names, 
                gene.symbols, gene.descs, gene.list.loc, gene.s2n, 
                gene.RES, core.enrichment))
            names(gene.report) <- c("#", "PROBE_ID", "SYMBOL", 
                "DESC", "LIST LOC", "S2N", "RES", "CORE_ENRICHMENT")
            if (output.directory != "") {
                filename <- paste(output.directory, doc.string, 
                  ".", gs.names[i], ".report.txt", sep = "", 
                  collapse = "")
                write.table(gene.report, file = filename, quote = F, 
                  row.names = F, sep = "\t")
                if (non.interactive.run == F) {
                  if (.Platform$OS.type == "windows") {
                    gs.filename <- paste(output.directory, doc.string, 
                      ".", gs.names[i], sep = "", collapse = "")
                    x11(width = 18, height = 6)
                  }
                  else if (.Platform$OS.type == "unix") {
                    gs.filename <- paste(output.directory, doc.string, 
                      ".", gs.names[i], ".pdf", sep = "", collapse = "")
                    pdf(file = gs.filename, height = 6, width = 18)
                  }
                }
                else {
                  if (.Platform$OS.type == "unix") {
                    gs.filename <- paste(output.directory, doc.string, 
                      ".", gs.names[i], ".pdf", sep = "", collapse = "")
                    pdf(file = gs.filename, height = 6, width = 18)
                  }
                  else if (.Platform$OS.type == "windows") {
                    gs.filename <- paste(output.directory, doc.string, 
                      ".", gs.names[i], ".pdf", sep = "", collapse = "")
                    pdf(file = gs.filename, height = 6, width = 18)
                  }
                }
            }
            nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), 
                1, c(1, 1, 1), TRUE)
            ind <- 1:N
            min.RES <- min(Obs.RES[i, ])
            max.RES <- max(Obs.RES[i, ])
            if (max.RES < 0.3) 
                max.RES <- 0.3
            if (min.RES > -0.3) 
                min.RES <- -0.3
            delta <- (max.RES - min.RES) * 0.5
            min.plot <- min.RES - 2 * delta
            max.plot <- max.RES
            max.corr <- max(obs.s2n)
            min.corr <- min(obs.s2n)
            Obs.correl.vector.norm <- (obs.s2n - min.corr)/(max.corr - 
                min.corr) * 1.25 * delta + min.plot
            zero.corr.line <- (-min.corr/(max.corr - min.corr)) * 
                1.25 * delta + min.plot
            col <- ifelse(Obs.ES[i] > 0, 2, 4)
            sub.string <- paste("Number of genes: ", N, " (in list), ", 
                size.G[i], " (in gene set)", sep = "", collapse = "")
            main.string <- paste("Gene Set ", i, ":", gs.names[i])
            plot(ind, Obs.RES[i, ], main = main.string, sub = sub.string, 
                xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", 
                xlim = c(1, N), ylim = c(min.plot, max.plot), 
                type = "l", lwd = 2, cex = 1, col = col)
            for (j in seq(1, N, 20)) {
                lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), 
                  lwd = 1, cex = 1, col = colors()[12])
            }
            lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, 
                col = 1)
            lines(c(Obs.arg.ES[i], Obs.arg.ES[i]), c(min.plot, 
                max.plot), lwd = 1, lty = 3, cex = 1, col = col)
            for (j in 1:N) {
                if (Obs.indicator[i, j] == 1) {
                  lines(c(j, j), c(min.plot + 1.25 * delta, min.plot + 
                    1.75 * delta), lwd = 1, lty = 1, cex = 1, 
                    col = 1)
                }
            }
            lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, 
                cex = 1, col = 1)
            lines(c(1, N), c(zero.corr.line, zero.corr.line), 
                lwd = 1, lty = 1, cex = 1, col = 1)
            temp <- order(abs(obs.s2n), decreasing = T)
            arg.correl <- temp[N]
            lines(c(arg.correl, arg.correl), c(min.plot, max.plot), 
                lwd = 1, lty = 3, cex = 1, col = 3)
            leg.txt <- paste("\"", phen1, "\" ", sep = "", collapse = "")
            text(x = 1, y = min.plot, adj = c(0, 0), labels = leg.txt, 
                cex = 1)
            leg.txt <- paste("\"", phen2, "\" ", sep = "", collapse = "")
            text(x = N, y = min.plot, adj = c(1, 0), labels = leg.txt, 
                cex = 1)
            adjx <- ifelse(Obs.ES[i] > 0, 0, 1)
            leg.txt <- paste("Peak at ", Obs.arg.ES[i], sep = "", 
                collapse = "")
            text(x = Obs.arg.ES[i], y = min.plot + 1.8 * delta, 
                adj = c(adjx, 0), labels = leg.txt, cex = 1)
            leg.txt <- paste("Zero crossing at ", arg.correl, 
                sep = "", collapse = "")
            text(x = arg.correl, y = min.plot + 1.95 * delta, 
                adj = c(adjx, 0), labels = leg.txt, cex = 1)
            sub.string <- paste("ES =", signif(Obs.ES[i], digits = 3), 
                " NES =", signif(Obs.ES.norm[i], digits = 3), 
                "Nom. p-val=", signif(p.vals[i, 1], digits = 3), 
                "FWER=", signif(p.vals[i, 2], digits = 3), "FDR=", 
                signif(FDR.mean.sorted[i], digits = 3))
            temp <- density(phi[i, ], adjust = adjust.param)
            x.plot.range <- range(temp$x)
            y.plot.range <- c(-0.125 * max(temp$y), 1.5 * max(temp$y))
            plot(temp$x, temp$y, type = "l", sub = sub.string, 
                xlim = x.plot.range, ylim = y.plot.range, lwd = 2, 
                col = 2, main = "Gene Set Null Distribution", 
                xlab = "ES", ylab = "P(ES)")
            x.loc <- which.min(abs(temp$x - Obs.ES[i]))
            lines(c(Obs.ES[i], Obs.ES[i]), c(0, temp$y[x.loc]), 
                lwd = 2, lty = 1, cex = 1, col = 1)
            lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, 
                col = 1)
            leg.txt <- c("Gene Set Null Density", "Observed Gene Set ES value")
            c.vec <- c(2, 1)
            lty.vec <- c(1, 1)
            lwd.vec <- c(2, 2)
            legend(x = -0.2, y = y.plot.range[2], bty = "n", 
                bg = "white", legend = leg.txt, lty = lty.vec, 
                lwd = lwd.vec, col = c.vec, cex = 1)
            leg.txt <- paste("Neg. ES \"", phen2, "\" ", sep = "", 
                collapse = "")
            text(x = x.plot.range[1], y = -0.1 * max(temp$y), 
                adj = c(0, 0), labels = leg.txt, cex = 1)
            leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep = "", 
                collapse = "")
            text(x = x.plot.range[2], y = -0.1 * max(temp$y), 
                adj = c(1, 0), labels = leg.txt, cex = 1)
            kk <- 1
            pinko <- matrix(0, nrow = size.G[i], ncol = cols)
            pinko.gene.names <- vector(length = size.G[i], mode = "character")
            for (k in 1:rows) {
                if (Obs.indicator[i, k] == 1) {
                  pinko[kk, ] <- A[obs.index[k], ]
                  pinko.gene.names[kk] <- obs.gene.symbols[k]
                  kk <- kk + 1
                }
            }
            GSEA.HeatMapPlot(V = pinko, row.names = pinko.gene.names, 
                col.labels = class.labels, col.classes = class.phen, 
                col.names = sample.names, main = " Heat Map for Genes in Gene Set", 
                xlab = " ", ylab = " ")
            if (non.interactive.run == F) {
                if (.Platform$OS.type == "windows") {
                  savePlot(filename = gs.filename, type = "jpeg", 
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
    return(list(report1 = report.phen1, report2 = report.phen2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
