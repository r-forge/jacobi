\name{jSimSVD}
\Rdversion{1.1}
\alias{jSimSVD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jSimSVD(x, eps = 1e-06, itmax = 1000, vectors = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, eps = 1e-06, itmax = 1000, vectors = TRUE, verbose = FALSE) 
{
    n <- dim(x)[1]
    m <- dim(x)[2]
    nmat <- dim(x)[3]
    itel <- 1
    kkk <- diag(n)
    lll <- diag(m)
    sxx <- sum(x^2)
    fold <- Inf
    print(dim(x))
    repeat {
        for (i in 1:(n - 1)) {
            if (i > m) 
                next()
            for (j in (i + 1):n) {
                v <- matrix(0, 2, 2)
                for (imat in 1:nmat) {
                  xij <- ifelse(j > m, 0, x[i, j, imat])
                  xjj <- ifelse(j > m, 0, x[j, j, imat])
                  xii <- x[i, i, imat]
                  xji <- x[j, i, imat]
                  v[1, 1] <- v[1, 1] + (xij^2 + xji^2)
                  v[1, 2] <- v[2, 1] <- v[1, 2] + (xii * xji - 
                    xjj * xij)
                  v[2, 2] <- v[2, 2] + (xii^2 + xjj^2)
                }
                u <- eigen(v)$vectors[, 1]
                for (imat in 1:nmat) {
                  xi <- x[i, , imat]
                  xj <- x[j, , imat]
                  x[i, , imat] <- u[2] * xi + u[1] * xj
                  x[j, , imat] <- u[2] * xj - u[1] * xi
                }
                if (vectors) {
                  ki <- kkk[i, ]
                  kj <- kkk[j, ]
                  kkk[i, ] <- u[2] * ki + u[1] * kj
                  kkk[j, ] <- u[2] * kj - u[1] * ki
                }
            }
        }
        ss <- sum(apply(x, 3, diag)^2)
        fnew <- sqrt((sxx - ss)/sxx)
        if (verbose) 
            cat(" Left iteration ", formatC(itel, digits = 4), 
                "loss ", formatC(fnew, digits = 6, width = 10), 
                "\n")
        for (k in 1:(m - 1)) {
            if (k > n) 
                next()
            for (l in (k + 1):m) {
                v <- matrix(0, 2, 2)
                for (imat in 1:nmat) {
                  xlk <- ifelse(l > n, 0, x[l, k, imat])
                  xll <- ifelse(l > n, 0, x[l, l, imat])
                  xkl <- x[k, l, imat]
                  xkk <- x[k, k, imat]
                  v[1, 1] <- v[1, 1] + (xkl^2 + xlk^2)
                  v[1, 2] <- v[2, 1] <- v[1, 2] + (xll * xlk - 
                    xkk * xkl)
                  v[2, 2] <- v[2, 2] + (xkk^2 + xll^2)
                }
                u <- eigen(v)$vectors[, 1]
                for (imat in 1:nmat) {
                  xk <- x[, k, imat]
                  xl <- x[, l, imat]
                  x[, k, imat] <- u[2] * xk - u[1] * xl
                  x[, l, imat] <- u[1] * xk + u[2] * xl
                }
                if (vectors) {
                  lk <- lll[, k]
                  ll <- lll[, l]
                  lll[, k] <- u[2] * lk - u[1] * ll
                  lll[, l] <- u[1] * lk + u[2] * ll
                }
            }
        }
        ss <- sum(apply(x, 3, diag)^2)
        fnew <- sqrt((sxx - ss)/sxx)
        if (verbose) 
            cat("Right iteration ", formatC(itel, digits = 4), 
                "loss ", formatC(fnew, digits = 6, width = 10), 
                "\n")
        if (((fold - fnew) < eps) || (itel == itmax)) 
            break()
        itel <- itel + 1
        fold <- fnew
    }
    return(list(d = apply(x, 3, diag), u = t(kkk), v = lll))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
