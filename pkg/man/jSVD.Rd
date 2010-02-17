\name{jSVD}
\Rdversion{1.1}
\alias{jSVD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jSVD(x, eps1 = 1e-06, eps2 = 1e-06, itmax = 1000, vectors = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, eps1 = 1e-06, eps2 = 1e-06, itmax = 1000, vectors = TRUE, 
    verbose = FALSE) 
{
    n <- nrow(x)
    m <- ncol(x)
    itel <- 1
    mx <- 0
    kkk <- diag(n)
    lll <- diag(m)
    sxx <- sum(x^2)
    sxm <- sqrt(sxx/(n * m))
    repeat {
        for (i in 1:(n - 1)) {
            if (i > m) 
                next()
            for (j in (i + 1):n) {
                xi <- x[i, ]
                xj <- x[j, ]
                xij <- ifelse(j > m, 0, x[i, j])
                xjj <- ifelse(j > m, 0, x[j, j])
                xii <- x[i, i]
                xji <- x[j, i]
                mx <- max(mx, abs(xij)/sxm, abs(xji)/sxm)
                v <- matrix(0, 2, 2)
                v[1, 1] <- xij^2 + xji^2
                v[1, 2] <- v[2, 1] <- xii * xji - xjj * xij
                v[2, 2] <- xii^2 + xjj^2
                u <- eigen(v)$vectors[, 1]
                x[i, ] <- u[2] * xi + u[1] * xj
                x[j, ] <- u[2] * xj - u[1] * xi
                if (vectors) {
                  ki <- kkk[i, ]
                  kj <- kkk[j, ]
                  kkk[i, ] <- u[2] * ki + u[1] * kj
                  kkk[j, ] <- u[2] * kj - u[1] * ki
                }
            }
        }
        ff <- sqrt((sxx - sum(diag(x)^2))/sxx)
        if (verbose) 
            cat(" Left iteration ", formatC(itel, digits = 4), 
                "maxel ", formatC(mx, width = 10), "loss ", formatC(ff, 
                  width = 10), "\n")
        for (k in 1:(m - 1)) {
            if (k > n) 
                next()
            for (l in (k + 1):m) {
                xk <- x[, k]
                xl <- x[, l]
                xlk <- ifelse(l > n, 0, x[l, k])
                xll <- ifelse(l > n, 0, x[l, l])
                xkk <- x[k, k]
                xkl <- x[k, l]
                mx <- max(mx, abs(xkl)/sxm, abs(xlk)/sxm)
                v <- matrix(0, 2, 2)
                v[1, 1] <- xkl^2 + xlk^2
                v[1, 2] <- v[2, 1] <- xll * xlk - xkk * xkl
                v[2, 2] <- xkk^2 + xll^2
                u <- eigen(v)$vectors[, 1]
                x[, k] <- u[2] * xk - u[1] * xl
                x[, l] <- u[1] * xk + u[2] * xl
                if (vectors) {
                  lk <- lll[, k]
                  ll <- lll[, l]
                  lll[, k] <- u[2] * lk - u[1] * ll
                  lll[, l] <- u[1] * lk + u[2] * ll
                }
            }
        }
        ff <- sqrt((sxx - sum(diag(x)^2))/sxx)
        if (verbose) 
            cat("Right iteration ", formatC(itel, digits = 4), 
                "maxel ", formatC(mx, width = 10), "loss ", formatC(ff, 
                  width = 10), "\n")
        if ((mx < eps1) || (ff < eps2) || (itel == itmax)) 
            break()
        itel <- itel + 1
        mx <- 0
    }
    return(list(d = diag(x), u = t(kkk), v = lll))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
