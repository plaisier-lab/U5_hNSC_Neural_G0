\name{OPAM.compute.models.overlap}
\alias{OPAM.compute.models.overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
OPAM.compute.models.overlap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OPAM.compute.models.overlap(model, model.set = NA, models.dir, thres = 0.05, size.lim = 500, results.file, produce.overlap.models = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{model.set}{
%%     ~~Describe \code{model.set} here~~
}
  \item{models.dir}{
%%     ~~Describe \code{models.dir} here~~
}
  \item{thres}{
%%     ~~Describe \code{thres} here~~
}
  \item{size.lim}{
%%     ~~Describe \code{size.lim} here~~
}
  \item{results.file}{
%%     ~~Describe \code{results.file} here~~
}
  \item{produce.overlap.models}{
%%     ~~Describe \code{produce.overlap.models} here~~
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
function (model, model.set = NA, models.dir, thres = 0.05, size.lim = 500, 
    results.file, produce.overlap.models = F) 
{
    file.list <- list.files(models.dir)
    files <- file.list[regexpr(pattern = ".mod", file.list) > 
        1]
    if (!is.na(model.set[1])) {
        model.set <- c(model, model.set)
        model.set <- paste(model.set, ".mod", sep = "")
        files <- intersect(files, model.set)
    }
    print(c("files=", files))
    max.models <- length(files)
    model.sizes <- NULL
    model.genes <- NULL
    model.names <- NULL
    num.models <- 0
    for (model.i in 1:max.models) {
        m.file <- paste(models.dir, files[model.i], sep = "")
        file.content <- readLines(m.file, n = -1)
        len <- length(file.content)
        for (i in 1:len) {
            temp <- unlist(strsplit(file.content[[i]], "\t"))
            len.param <- length(temp)
            if (len.param == 2) {
                param.vals <- temp[2]
            }
            else {
                param.vals <- paste(noquote(temp[2:len.param]), 
                  collapse = ",")
                param.vals <- paste("c(", param.vals, ")", sep = "")
            }
            assignment.string <- paste(noquote(temp[1]), " <- ", 
                param.vals, sep = "")
            eval(parse(text = assignment.string))
        }
        if (exists("msig.up.genes3")) 
            msig.up.genes <- msig.up.genes3
        m.size <- length(msig.up.genes)
        if (m.size <= size.lim) {
            model.names <- c(model.names, model.name)
            model.sizes <- c(model.sizes, m.size)
            model.genes <- rbind(model.genes, c(msig.up.genes, 
                rep(NA, size.lim - m.size)))
            num.models <- num.models + 1
        }
        rm(model.name, m.size, msig.up.genes, msig.up.genes3)
    }
    print(c("model:", model, " model.names:", model.names))
    loc <- match(model, model.names)
    if (is.na(loc)) 
        stop(paste("model:", model, " is not in models directory: ", 
            models.dir))
    m.size <- model.sizes[loc]
    m.genes <- model.genes[loc, 1:m.size]
    model.overlap.size <- vector(length = num.models, mode = "numeric")
    model.overlap.genes <- matrix(NA, nrow = num.models, ncol = size.lim)
    model.overlap.signif <- vector(length = num.models, mode = "numeric")
    for (i in 1:num.models) {
        overlap <- intersect(m.genes, model.genes[i, 1:model.sizes[i]])
        if (length(overlap) == 0) {
            model.overlap.size[i] <- 0
        }
        else {
            model.overlap.size[i] <- length(overlap)
            model.overlap.genes[i, 1:length(overlap)] <- overlap
            model.overlap.signif[i] <- (length(overlap)/m.size) * 
                (length(overlap)/model.sizes[i])
        }
    }
    save.flag <- ifelse(model.overlap.signif >= thres, 1, 0)
    total.save <- sum(save.flag)
    results.tab <- cbind(model.names, model.sizes, rep(m.size, 
        num.models), model.overlap.size, model.overlap.signif, 
        save.flag)
    ind <- order(model.overlap.signif, decreasing = T)
    results.tab <- results.tab[ind, ]
    colnames(results.tab) <- c("Name", "Size", "Model Size", 
        "Overlap", "Score", "Save Model")
    results.file <- paste(models.dir, "OVER_", model, ".txt", 
        sep = "")
    write.table(results.tab, file = results.file, quote = F, 
        row.names = F, sep = "\t")
    print(noquote(results.tab[1:(min(50, num.models)), ]))
    if (produce.overlap.models == T) {
        for (k in 1:total.save) {
            if (model == model.names[ind[k]]) 
                next
            print(paste("Creating model for model:", model.names[ind[k]], 
                sep = " "))
            gene.set <- model.overlap.genes[ind[k], 1:model.overlap.size[ind[k]]]
            gene.set.string <- paste("c('", paste(gene.set, collapse = "','"), 
                "')\n", sep = "")
            overlap.model.name <- paste("OVER_", model, "_", 
                model.names[ind[k]], sep = "")
            print(paste("Saving model: ", overlap.model.name))
            print(gene.set.string)
            print("------------")
            m.file <- paste(models.dir, overlap.model.name, ".mod", 
                sep = "")
            cat("model.creation.date", paste("'", date(), "'\n", 
                sep = ""), file = m.file, append = FALSE, sep = "\t")
            cat("model.name", paste("'", overlap.model.name, 
                "'\n", sep = ""), file = m.file, append = TRUE, 
                sep = "\t")
            cat("model.description", paste("'Overlap Model'\n", 
                sep = ""), file = m.file, append = TRUE, sep = "\t")
            cat("sample.norm.type", "'rank'\n", file = m.file, 
                append = TRUE, sep = "\t")
            cat("statistic", "'area.under.RES'\n", file = m.file, 
                append = TRUE, sep = "\t")
            cat("weight", paste(0.25, "\n", sep = ""), file = m.file, 
                append = TRUE, sep = "\t")
            cat("random.seed", paste(12345, "\n", sep = ""), 
                file = m.file, append = TRUE, sep = "\t")
            cat("nperm", paste(0, "\n", sep = ""), file = m.file, 
                append = TRUE, sep = "\t")
            cat("link.function", "'logit'\n", file = m.file, 
                append = TRUE, sep = "\t")
            cat("c1", "c('black','lightgrey')\n", file = m.file, 
                append = TRUE, sep = "\t")
            cat("msig.up.genes", gene.set.string, file = m.file, 
                append = TRUE, sep = "\t")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
