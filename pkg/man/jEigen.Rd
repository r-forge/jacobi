\name{jEigen}
\Rdversion{1.1}
\alias{jEigen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Eigenvectors of a matrix using Jacobi plane rotations.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jEigen(a, eps1 = 1e-06, eps2 = 1e-10, itmax = 100, vectors = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{eps1}{
%%     ~~Describe \code{eps1} here~~
}
  \item{eps2}{
%%     ~~Describe \code{eps2} here~~
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
function (a, eps1 = 1e-06, eps2 = 1e-10, itmax = 100, vectors = TRUE, 
    verbose = FALSE) 
{
    n <- nrow(a)
    k <- diag(n)
    itel <- 1
    mx <- 0
    saa <- sum(a^2)
    repeat {
        for (i in 1:(n - 1)) for (j in (i + 1):n) {
            aij <- a[i, j]
            bij <- abs(aij)
            if (bij < eps1) 
                next()
            mx <- max(mx, bij)
            am <- (a[i, i] - a[j, j])/2
            u <- c(aij, -am)
            u <- u/sqrt(sum(u^2))
            c <- sqrt((1 + u[2])/2)
            s <- sign(u[1]) * sqrt((1 - u[2])/2)
            ss <- s^2
            cc <- c^2
            sc <- s * c
            ai <- a[i, ]
            aj <- a[j, ]
            aii <- a[i, i]
            ajj <- a[j, j]
            a[i, ] <- a[, i] <- c * ai - s * aj
            a[j, ] <- a[, j] <- s * ai + c * aj
            a[i, j] <- a[j, i] <- 0
            a[i, i] <- aii * cc + ajj * ss - 2 * sc * aij
            a[j, j] <- ajj * cc + aii * ss + 2 * sc * aij
            if (vectors) {
                ki <- k[, i]
                kj <- k[, j]
                k[, i] <- c * ki - s * kj
                k[, j] <- s * ki + c * kj
            }
        }
        ff <- sqrt(saa - sum(diag(a)^2))
        if (verbose) 
            cat("Iteration ", formatC(itel, digits = 4), "maxel ", 
                formatC(mx, width = 10), "loss ", formatC(ff, 
                  width = 10), "\n")
        if ((mx < eps1) || (ff < eps2) || (itel == itmax)) 
            break()
        itel <- itel + 1
        mx <- 0
    }
    d <- diag(a)
    o <- order(d, decreasing = TRUE)
    if (vectors) 
        return(list(values = d[o], vectors = k[, o]))
    else return(values = d[o])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
