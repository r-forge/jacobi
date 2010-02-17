\name{jMCA}
\Rdversion{1.1}
\alias{jMCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
jMCA(burt, k, eps = 1e-06, itmax = 500, verbose = TRUE, vectors = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{burt}{
%%     ~~Describe \code{burt} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{itmax}{
%%     ~~Describe \code{itmax} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{vectors}{
%%     ~~Describe \code{vectors} here~~
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
function (burt, k, eps = 1e-06, itmax = 500, verbose = TRUE, 
    vectors = TRUE) 
{
    m <- length(k)
    burt <- m * m * burt/sum(burt)
    sk <- sum(k)
    db <- diag(burt)
    ll <- kk <- ww <- diag(sk)
    itel <- 1
    ossq <- 0
    klw <- 1 + cumsum(c(0, k))[1:m]
    kup <- cumsum(k)
    ind <- lapply(1:m, function(i) klw[i]:kup[i])
    sburt <- burt/sqrt(outer(db, db))
    for (i in 1:m) kk[ind[[i]], ind[[i]]] <- t(svd(sburt[ind[[i]], 
        ])$u)
    kbk <- kk \%*\% sburt \%*\% t(kk)
    for (i in 1:m) for (j in 1:m) ww[ind[[i]], ind[[j]]] <- ifelse(outer(1:k[i], 
        1:k[j], "=="), 1, 0)
    repeat {
        for (l in 1:m) {
            if (k[l] == 2) 
                next()
            li <- ind[[l]]
            for (i in (klw[l] + 1):(kup[l] - 1)) for (j in (i + 
                1):kup[l]) {
                bi <- kbk[i, -li]
                bj <- kbk[j, -li]
                wi <- ww[i, -li]
                wj <- ww[j, -li]
                acc <- sum(wi * bi^2) + sum(wj * bj^2)
                acs <- sum((wi - wj) * bi * bj)
                ass <- sum(wi * bj^2) + sum(wj * bi^2)
                u <- eigen(matrix(c(acc, acs, acs, ass), 2, 2))$vectors[, 
                  1]
                c <- u[1]
                s <- u[2]
                kbk[-li, i] <- kbk[i, -li] <- c * bi + s * bj
                kbk[-li, j] <- kbk[j, -li] <- c * bj - s * bi
                if (vectors) {
                  ki <- kk[i, li]
                  kj <- kk[j, li]
                  kk[i, li] <- c * ki + s * kj
                  kk[j, li] <- c * kj - s * ki
                }
            }
        }
        nssq <- sum(ww * kbk^2) - sum(diag(kbk)^2)
        if (verbose) 
            cat("Iteration ", formatC(itel, digits = 4), "ssq ", 
                formatC(nssq, digits = 10, width = 15), "\n")
        if (((nssq - ossq) < eps) || (itel == itmax)) 
            break()
        itel <- itel + 1
        ossq <- nssq
    }
    kl <- unlist(sapply(k, function(i) 1:i))
    pp <- ifelse(outer(1:sk, order(kl), "=="), 1, 0)
    pkbkp <- t(pp) \%*\% kbk \%*\% pp
    pk <- t(pp) \%*\% kk
    km <- as.vector(table(kl))
    nm <- length(km)
    klw <- 1 + cumsum(c(0, km))[1:nm]
    kup <- cumsum(km)
    for (i in 1:length(km)) {
        if (km[i] == 1) 
            next()
        ind <- klw[i]:kup[i]
        ll[ind, ind] <- eigen(pkbkp[ind, ind])$vectors
    }
    lpkbkpl <- t(ll) \%*\% pkbkp \%*\% ll
    lpk <- t(ll) \%*\% pk
    return(list(kbk = kbk, pkbkp = pkbkp, lpkbkpl = lpkbkpl, 
        kk = t(kk), pp = pp, ll = ll, kpl = t(lpk)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
