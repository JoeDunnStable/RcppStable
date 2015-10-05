\name{StableDistribution}

\title{Fitting Stable Distribution Function to Sample Data}

\alias{stable_fit}
\alias{dstable.quick}

\concept{Stable Distribution}

\description{
    \code{stable_fit} estimates the four parameters of the stable distribution
    that best fit the supplied sample data.  It uses either McCulloch's method,
    which is based on quantiles or a combination of McCulloch's followed by
    a maximum likelihood estimation.  For large datasets \code{dstable.quick}
    is used, which is faster, but less accurate than \code{dstable}
}


\usage{
stable_fit(y,type="q",quick=T,pm=0)

dstable.quick(x, alpha, beta, gamma=1, delta=0, pm=0,
        log=F,
        tol = 64 * .Machine$double.eps, zeta.tol = NULL,
        subdivisions = 1000)
}

\arguments{
  \item{y}{the sample data for stable_fit}

  \item{type}{the method used to obtain the parameters.  The choices are "q",
  which indicates McCulloch's quantile based method, or "mle", which indicates that
  a maximum likelihood method is to be used with an iniital quess from McCulloch's
  method.}

  \item{alpha, beta, gamma, delta}{
    value of the index parameter \code{alpha} with \code{alpha = (0,2]};
    skewness parameter \code{beta}, in the range [-1, 1];
    scale parameter \code{gamma}; and location (or \sQuote{shift})
    parameter \code{delta}.
  }

  \item{pm}{
    parameterization, an integer in \code{0, 1, 2}; by default \code{pm=0},
    the \sQuote{S0} parameterization.}

  \item{quick}{A boolean indicating whether \code{dstable.quick} is to be used in the
    maximum likelihood estimation.}

  \item{log}{logical; if TRUE, probabilities p are given as log(p).}

  \item{tol}{numerical tolerance,
    \describe{
      \item{dstable.quick():}{used for numerical integration,}
    }
  }
  \item{zeta.tol}{(\code{dstable}) numerical tolerance for checking if
    \code{x} is close to \eqn{\zeta(\alpha,\beta)}.  The default,
    \code{NULL} depends itself on \eqn{(\alpha,\beta)}.
    \cr
    As it is \emph{experimental} and not guaranteed to remain in the
    future, its use is not recommended in production code.  Rather e-mail the
    package maintainer about it.}
  \item{subdivisions}{maximal number of intervals for integration, see
    \code{\link{integrate}}.}
}

\value{
    Function \code{dstable.quick} returns a numeric vectors with the density.

    Function \code{stable.fit} returns a small database containing the fitted
    parameters
}

\details{
    \bold{Skew Stable Distribution:}
    \cr\cr
    The function \code{dstable.quick} selects a subset of the input sample, \code{y}
    and uses \code{dstable} to calculate the log of the skew stable distribution probability.
    These probabilites are then compared to the corresponding log probabilities for the Student
    t distribution with \code{df=alpha} and a natural cubic spline is fit through the subsample      points and used to estimate the log probabiliites for the rest of the sample.
}
\seealso{
  the \code{\link{dstable}()} function for details and
}
\author{
  Joseph Dunn.
}

\examples{
## stable -

## Plot stable random number series
   set.seed(1953)
   r <- rstable(n = 1000, alpha = 1.9, beta = 0.3)
   plot(r, type = "l", main = "stable: alpha=1.9 beta=0.3",
        col = "steelblue")
   grid()

## Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white",
        col = "steelblue")
   x <- seq(-5, 5, 0.25)
   lines(x, dstable(x, alpha = 1.9, beta = 0.3, tol= 1e-3), lwd = 2)

## Plot df and compare with true df:
   plot(ecdf(r), do.points=TRUE, col = "steelblue",
        main = "Probabilities:  ecdf(rstable(1000,..)) and true  cdf F()")
   rug(r)
   lines(x, pstable(q = x, alpha = 1.9, beta = 0.3),
         col="#0000FF88", lwd= 2.5)

## Switching  sign(beta)  <==> Mirror the distribution around  x == delta:
curve(dstable(x, alpha=1.2, beta =  .8, gamma = 3, delta = 2), -10, 10)
curve(dstable(x, alpha=1.2, beta = -.8, gamma = 3, delta = 2),
      add=TRUE, col=2)
## or the same
curve(dstable(2*2-x, alpha=1.2, beta = +.8, gamma = 3, delta = 2),
      add=TRUE, col=adjustcolor("gray",0.2), lwd=5)
abline(v = 2, col = "gray", lty=2, lwd=2)
axis(1, at = 2, label = expression(delta == 2))

## Compute quantiles:
   x. <- -4:4
   px <- pstable(x., alpha = 1.9, beta = 0.3)
  (qs <- qstable(px, alpha = 1.9, beta = 0.3))
   stopifnot(all.equal(as.vector(qs), x., tol = 1e-5))
}

\keyword{distribution}
