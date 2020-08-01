\name{MSIG.2D.Plot}
\alias{MSIG.2D.Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSIG.2D.Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MSIG.2D.Plot(input.ds, input.cls = "", output.2D.proj.plot, output.heatmap.plot, title = "", non.interactive.run = F, heatmap.row.norm = F, heatmap.cmap.type = 1, col = c("grey3", "cadetblue", "darkgreen", "chartreuse2", "red1", "darkred", "orange", "blue2", "lightblue", "pink", "coral"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.ds}{
%%     ~~Describe \code{input.ds} here~~
}
  \item{input.cls}{
%%     ~~Describe \code{input.cls} here~~
}
  \item{output.2D.proj.plot}{
%%     ~~Describe \code{output.2D.proj.plot} here~~
}
  \item{output.heatmap.plot}{
%%     ~~Describe \code{output.heatmap.plot} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{non.interactive.run}{
%%     ~~Describe \code{non.interactive.run} here~~
}
  \item{heatmap.row.norm}{
%%     ~~Describe \code{heatmap.row.norm} here~~
}
  \item{heatmap.cmap.type}{
%%     ~~Describe \code{heatmap.cmap.type} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
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
function (input.ds, input.cls = "", output.2D.proj.plot, output.heatmap.plot, 
    title = "", non.interactive.run = F, heatmap.row.norm = F, 
    heatmap.cmap.type = 1, col = c("grey3", "cadetblue", "darkgreen", 
        "chartreuse2", "red1", "darkred", "orange", "blue2", 
        "lightblue", "pink", "coral")) 
{
    print("Running MSIG.2D.Plot...")
    print(heatmap.row.norm)
    dataset <- MSIG.Gct2Frame(filename = input.ds)
    m <- data.matrix(dataset$ds)
    gs.names <- dataset$row.names
    gs.descs <- dataset$descs
    sample.names <- dataset$names
    dim(m)
    Ns <- length(m[1, ])
    k.proj <- length(m[, 1])
    if (input.cls != "") {
        CLS <- ReadClsFile(file = input.cls)
        class.labels <- CLS$class.v
        class.phen <- CLS$phen
    }
    else {
        class.labels <- rep(1, Ns)
        class.phen <- "Samples"
    }
    S1 <- as.real(m[1, ])
    S2 <- as.real(m[2, ])
    range.S1 <- range(S1)
    range.S2 <- range(S2)
    c0 <- col
    c1 <- colors()[match(c0, colors())]
    color <- c1[class.labels]
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.2D.proj.plot
            x11(height = 14, width = 22)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.2D.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 22)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.2D.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 22)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.2D.proj.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = 14, width = 22)
        }
    }
    nf <- layout(matrix(c(1, 2), 1, 2, byrow = T), widths = c(3, 
        2), heights = 1, respect = FALSE)
    plot(S1, S2, xlim = range.S1, ylim = range.S2, type = "n", 
        main = paste(title, " -- 2D Plot", sep = ""), sub = input.ds)
    for (j in 1:Ns) {
        if (min(class.labels) == 0) {
            color.code <- c1[class.labels[j] + 1]
        }
        else {
            color.code <- c1[class.labels[j]]
        }
        points(S1[j], S2[j], pch = 22, type = "p", cex = 2, bg = color.code, 
            col = "black")
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
            savePlot(filename = plot.filename, type = "jpeg", 
                device = dev.cur())
        }
        else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    }
    else {
        dev.off()
    }
    height <- ifelse(k.proj > 50, 15, 0.2 * k.proj + 4.8)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            plot.filename <- output.heatmap.plot
            x11(height = height, width = 15)
        }
        else if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
    }
    else {
        if (.Platform$OS.type == "unix") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
        else if (.Platform$OS.type == "windows") {
            plot.filename <- paste(output.heatmap.plot, ".pdf", 
                sep = "", collapse = "")
            pdf(file = plot.filename, height = height, width = 15)
        }
    }
    MSIG.HeatMapPlot(V = m, row.names = gs.names, col.labels = class.labels, 
        col.classes = class.phen, col.names = sample.names, main = paste(title, 
            " -- Heat Map", sep = ""), sub = " ", xlab = " ", 
        ylab = " ", row.norm = heatmap.row.norm, cmap.type = heatmap.cmap.type)
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = plot.filename, type = "jpeg", 
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
