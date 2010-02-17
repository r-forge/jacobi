\name{jTucker3Diag}
\Rdversion{1.1}
\alias{jTucker3Diag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jTucker3Diag(a, eps = 1e-06, itmax = 100, vectors = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{itmax}{
%%     ~~Describe \code{itmax} here~~
}
  \item{vectors}{
%%     ~~Describe \code{vectors} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (a, eps = 1e-06, itmax = 100, vectors = TRUE, verbose = TRUE) 
{
    n <- dim(a)[1]
    m <- dim(a)[2]
    k <- dim(a)[3]
    nmk <- min(n, m, k)
    kn <- diag(n)
    km <- diag(m)
    kk <- diag(k)
    ossq <- 0
    itel <- 1
    repeat {
        for (i in 1:(n - 1)) for (j in (i + 1):n) {
            ai <- a[i, , ]
            aj <- a[j, , ]
            acc <- ass <- asc <- 0
            if (i <= min(m, k)) {
                acc <- acc + a[i, i, i]^2
                ass <- ass + a[j, i, i]^2
                asc <- asc + a[i, i, i] * a[j, i, i]
            }
            if (j <= min(m, k)) {
                acc <- acc + a[j, j, j]^2
                ass <- ass + a[i, j, j]^2
                asc <- asc - a[j, j, j] * a[i, j, j]
            }
            u <- eigen(matrix(c(acc, asc, asc, ass), 2, 2))$vectors[, 
                1]
            c <- u[1]
            s <- u[2]
            a[i, , ] <- c * ai + s * aj
            a[j, , ] <- c * aj - s * ai
            if (vectors) {
                ki <- kn[i, ]
                kj <- kn[j, ]
                kn[i, ] <- c * ki + s * kj
                kn[j, ] <- c * kj - s * ki
            }
        }
        for (i in 1:(m - 1)) for (j in (i + 1):m) {
            ai <- a[, i, ]
            aj <- a[, j, ]
            acc <- ass <- asc <- 0
            if (i <= min(n, k)) {
                acc <- acc + a[i, i, i]^2
                ass <- ass + a[i, j, i]^2
                asc <- asc + a[i, i, i] * a[i, j, i]
            }
            if (j <= min(n, k)) {
                acc <- acc + a[j, j, j]^2
                ass <- ass + a[j, i, j]^2
                asc <- asc - a[j, j, j] * a[j, i, j]
            }
            u <- eigen(matrix(c(acc, asc, asc, ass), 2, 2))$vectors[, 
                1]
            c <- u[1]
            s <- u[2]
            a[, i, ] <- c * ai + s * aj
            a[, j, ] <- c * aj - s * ai
            if (vectors) {
                ki <- km[i, ]
                kj <- km[j, ]
                km[i, ] <- c * ki + s * kj
                km[j, ] <- c * kj - s * ki
            }
        }
        for (i in 1:(k - 1)) for (j in (i + 1):k) {
            ai <- a[, , i]
            aj <- a[, , j]
            acc <- ass <- asc <- 0
            if (i <= min(n, m)) {
                acc <- acc + a[i, i, i]^2
                ass <- ass + a[i, i, j]^2
                asc <- asc + a[i, i, i] * a[i, i, j]
            }
            if (j <= min(n, m)) {
                acc <- acc + a[j, j, j]^2
                ass <- ass + a[j, j, i]^2
                asc <- asc - a[j, j, j] * a[j, j, i]
            }
            u <- eigen(matrix(c(acc, asc, asc, ass), 2, 2))$vectors[, 
                1]
            c <- u[1]
            s <- u[2]
            a[, , i] <- c * ai + s * aj
            a[, , j] <- c * aj - s * ai
            if (vectors) {
                ki <- kk[i, ]
                kj <- kk[j, ]
                kk[i, ] <- c * ki + s * kj
                kk[j, ] <- c * kj - s * ki
            }
        }
        nssq <- 0
        for (v in 1:nmk) nssq <- nssq + a[v, v, v]^2
        if (verbose) 
            cat("Iteration ", formatC(itel, digits = 4), "ssq ", 
                formatC(nssq, digits = 10, width = 15), "\n")
        if (((nssq - ossq) < eps) || (itel == itmax)) 
            break()
        itel <- itel + 1
        ossq <- nssq
    }
    d <- rep(0, nmk)
    for (v in 1:nmk) d[v] <- a[v, v, v]
    return(list(a = a, d = d, kn = kn, km = km, kk = kk))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
