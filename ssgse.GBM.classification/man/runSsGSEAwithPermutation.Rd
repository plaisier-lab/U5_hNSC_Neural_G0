\name{runSsGSEAwithPermutation}
\alias{runSsGSEAwithPermutation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
runSsGSEAwithPermutation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runSsGSEAwithPermutation(profile_data_file, number_perms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{profile_data_file}{
%%     ~~Describe \code{profile_data_file} here~~
}
  \item{number_perms}{
%%     ~~Describe \code{number_perms} here~~
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
function (profile_data_file, number_perms) 
{
    MOD_file = "MOD"
    profile_data_path = dirname(file.path(profile_data_file))
    profile_data_name = basename(profile_data_file)
    data <- read.table(paste(profile_data_path, profile_data_name, 
        sep = "/"), header = T, row.names = 1, skip = 2)
    size <- dim(data[-1])
    num_rows <- size[1]
    num_cols <- size[2]
    num_perm <- number_perms
    colnames <- colnames(data)
    print("Data loading finished!")
    line1 <- "#1.2"
    line2 <- paste(num_rows, num_perm, sep = "\t")
    line3 <- paste("X", 1:num_perm, sep = "", collapse = "\t")
    line3 <- paste("NAME", "Description", line3, sep = "\t")
    write.table(line1, paste(profile_data_path, "/random_profile_", 
        profile_data_name, ".gct", sep = ""), quote = F, col.name = F, 
        row.name = F)
    write.table(line2, paste(profile_data_path, "/random_profile_", 
        profile_data_name, ".gct", sep = ""), quote = F, col.name = F, 
        row.name = F, append = T)
    write.table(line3, paste(profile_data_path, "/random_profile_", 
        profile_data_name, ".gct", sep = ""), quote = F, col.name = F, 
        row.name = F, append = T)
    random_profile <- data[1]
    for (i in 1:num_perm) {
        a <- data[-1][cbind(seq(1:num_rows), sample(1:num_cols, 
            num_rows, replace = T))]
        random_profile <- cbind(random_profile, a)
        if (i\%\%100 == 0) {
            print(i)
        }
    }
    write.table(random_profile, paste(profile_data_path, "/random_profile_", 
        profile_data_name, ".gct", sep = ""), sep = "\t", quote = F, 
        col.name = F, append = T)
    print("Random profiles was genereated!")
    selected.models <- c("Proneural", "Classical", "Mesenchymal")
    OPAM.apply.model.2(input.ds = paste(profile_data_path, "/random_profile_", 
        profile_data_name, ".gct", sep = ""), models.dir = MOD_file, 
        models = selected.models, raw.score.outfile = paste(profile_data_path, 
            "/random_raw.score_", profile_data_name, ".gct", 
            sep = ""), norm.score.outfile = "", model.score.outfile = "", 
        prob.outfile = "", graphics.off = T)
    print("SsGSEA was performed on random profiles!")
    OPAM.apply.model.2(input.ds = paste(profile_data_path, profile_data_name, 
        sep = "/"), models.dir = MOD_file, models = selected.models, 
        raw.score.outfile = paste(profile_data_path, "/raw.score_", 
            profile_data_name, ".gct", sep = ""), norm.score.outfile = "", 
        model.score.outfile = "", prob.outfile = "", graphics.off = T)
    print("SsGSEA was performed on the original profiles!")
    random_result <- read.table(paste(profile_data_path, "/random_raw.score_", 
        profile_data_name, ".gct", sep = ""), header = T, row.names = 1, 
        skip = 2)
    random_result <- random_result[-1]
    random_result <- t(random_result)
    original_result <- read.table(paste(profile_data_path, "/raw.score_", 
        profile_data_name, ".gct", sep = ""), header = T, row.names = 1, 
        skip = 2)
    original_result <- original_result[-1]
    original_result <- t(original_result)
    p_result <- original_result
    for (i in 1:dim(original_result)[1]) {
        p_result[i, ] <- colSums(sweep(random_result, 2, original_result[i, 
            ]) >= 0)
    }
    write.table(cbind(original_result, p_result), paste(profile_data_path, 
        "/p_result_", profile_data_name, ".txt", sep = ""), sep = "\t", 
        quote = F)
    print("P_values for each subtype was calculated!")
    print("Finished!")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
