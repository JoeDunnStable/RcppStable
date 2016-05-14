\name{stable_fit}

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
        log=FALSE,
        tol = 64 * .Machine$double.eps, zeta.tol = NULL,
        subdivisions = 1000)
}

\arguments{
  \item{y}{the sample data for stable_fit}

  \item{type}{the method used to obtain the parameters.  The choices are "q",
  which indicates McCulloch's quantile based method, or "mle", which indicates that
  a maximum likelihood method is to be used with an iniital quess from McCulloch's
  method.}

  \item{x}{numeric vector of quantiles for dstable.quick}

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

    Function \code{stable.fit} returns a list with two items: \code{parameters}, a small          database containing the fitted parameters, and \code{fit_mle}, the output of stats4::mle      when the \code{type} is "mle" or \code{NULL} when the \code{type} is "q".
}

\details{
    \bold{Skew Stable Distribution:}
    \cr\cr
    The function \code{dstable.quick} uses \code{dstable} to evaluate (1) all points \code{y}
    close to the mode, (2) a few points at the extreme tail of \code{y}, and (3) a                sample of points away from the mode or the extreme tail. The log probabilites for points      in category (3) are then compared to the corresponding log probabilities for the Student
    t distribution with \code{df=alpha} and a natural cubic spline if fit using the ratios so     calculated as knots.  The cubic spline is then used to estimate the log probabiliites for     the rest of the points \code{y} that are not in categories (1) and (2).
}

\seealso{
  the \code{\link{dstable}()} function for details and references.
}
\author{
  Joseph Dunn.
}

\examples{
## stable_fit

## Generate series of stable random numbers with known alpha and beta.
   set.seed(1953)
   r <- rstable(n = 10000, alpha = 1.9, beta = 0.3)

## Calculate the fit using McCulloch's method
   fit_q<-stable_fit(r,type="q")
   print(fit_q)

## Calculate the fit using McCulloch's method for the initial quess, followed
## by the Maximum Likelihood Method.
   fit_mle<-stable_fit(r,type="mle")
   print(fit_mle)
}

\keyword{distribution}
