\name{jSimDiag}
\Rdversion{1.1}
\alias{jSimDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jSimDiag(a, eps = 1e-10, itmax = 100, vectors = TRUE, verbose = FALSE)
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
function (a, eps = 1e-10, itmax = 100, vectors = TRUE, verbose = FALSE) 
{
    n <- dim(a)[1]
    kk <- diag(n)
    m <- dim(a)[3]
    itel <- 1
    saa <- sum(a^2)
    fold <- saa - sum(apply(a, 3, function(x) sum(diag(x^2))))
    repeat {
        for (i in 1:(n - 1)) for (j in (i + 1):n) {
            ad <- (a[i, i, ] - a[j, j, ])/2
            av <- a[i, j, ]
            v <- matrix(0, 2, 2)
            v[1, 1] <- sum(ad^2)
            v[1, 2] <- v[2, 1] <- sum(av * ad)
            v[2, 2] <- sum(av^2)
            u <- eigen(v)$vectors[, 2]
            c <- sqrt((1 + u[2])/2)
            s <- sign(u[1]) * sqrt((1 - u[2])/2)
            for (k in 1:m) {
                ss <- s^2
                cc <- c^2
                sc <- s * c
                ai <- a[i, , k]
                aj <- a[j, , k]
                aii <- a[i, i, k]
                ajj <- a[j, j, k]
                aij <- a[i, j, k]
                a[i, , k] <- a[, i, k] <- c * ai - s * aj
                a[j, , k] <- a[, j, k] <- s * ai + c * aj
                a[i, j, k] <- a[j, i, k] <- u[1] * (aii - ajj)/2 + 
                  u[2] * aij
                a[i, i, k] <- aii * cc + ajj * ss - 2 * sc * 
                  aij
                a[j, j, k] <- ajj * cc + aii * ss + 2 * sc * 
                  aij
            }
            if (vectors) {
                ki <- kk[, i]
                kj <- kk[, j]
                kk[, i] <- c * ki - s * kj
                kk[, j] <- s * ki + c * kj
            }
        }
        fnew <- saa - sum(apply(a, 3, function(x) sum(diag(x^2))))
        if (verbose) 
            cat("Iteration ", formatC(itel, digits = 4), "old loss ", 
                formatC(fold, width = 10), "new loss ", formatC(fnew, 
                  width = 10), "\n")
        if (((fold - fnew) < eps) || (itel == itmax)) 
            break()
        itel <- itel + 1
        fold <- fnew
    }
    return(list(a = a, d <- apply(a, 3, diag), k = kk))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
