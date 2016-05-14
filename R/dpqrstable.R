

pi2<-pi/2

##' @title omega() according to Lambert & Lindsey (1999), p.412
##' @param gamma [dpqr]stable()'s scale parameter, > 0         -- of length 1
##' @param alpha [dpqr]stable()'s "main" parameter, in [0, 2]  -- of length 1
##' @return omega(.) = tan(pi/2 alpha) if alpha != 1 ...
.om <- function(gamma,alpha) {
  if(alpha != round(alpha)) # non-integer usual case
    tan(pi2*alpha)# not tanpi2() !
  else if(alpha == 1)
    (2/pi)*log(gamma)
  else 0 # for alpha = 0 or = 2
}

dstable<-function (x, alpha, beta, gamma = 1, delta = 0, pm = 0, log = FALSE,
          tol = 64 * .Machine$double.eps, zeta.tol = NULL, subdivisions = 1000) {
  ## Original implemented by Diethelm Wuertz;
  ## Changes for efficiency and accuracy by Martin Maechler

  ## Description:
  ##	 Returns density for stable DF

  ## Details:
  ##	 The function uses the approach of J.P. Nolan for general
  ##	 stable distributions. Nolan derived expressions in form
  ##	 of integrals based on the charcteristic function for
  ##	 standardized stable random variables. These integrals
  ##	 can be numerically evaluated.

  ## Arguments:
  ##	 alpha = index of stability, in the range (0,2]
  ##	 beta  = skewness, in the range [-1, 1]
  ##	 gamma = scale, in the range (0, infinity)
  ##	 delta = location, in the range (-infinity, +infinity)
  ##	 pm = type of parmeterization
  ##   lower.tail = TRUE for lower tail, FALSE for upper tail
  ##   log_flag = return log of density

  ## Note: S+ compatibility no longer considered (explicitly)

  ## Parameter Check:
  ## NB: (gamma, delta) can be *vector*s (vectorized along x)
  stopifnot(0 < alpha, alpha <= 2, length(alpha) == 1, -1 <=
              beta, beta <= 1, length(beta) == 1, 0 <= gamma, length(pm) ==
              1, pm %in% 0:2, tol > 0, subdivisions > 0)
  ## not an official argument {no doc!}:
  verbose <- getOption("dstable.debug", default = FALSE)
  ## Parameterizations:
  if (pm == 1) {
    delta <- delta + beta * gamma * .om(gamma, alpha)
  }
  else if (pm == 2) {
    delta <- delta - alpha^(-1/alpha) * gamma * stableMode(alpha,
                                                           beta)
    gamma <- alpha^(-1/alpha) * gamma
  } ## else pm == 0

  ## Shift and Scale:
  x <- (x - delta)/gamma
  ## Special Cases:
  if (alpha == 2) {
    ans<-dnorm(x, mean = 0, sd = sqrt(2), log = log)
  }
  else if (alpha == 1 && beta == 0) {
    ans<-dcauchy(x, log = log)
  }
  else {
    if (is.null(zeta.tol))
      zeta.tol<-0.
    ans<-as.vector(sdstable(x,alpha, beta, log,
                         tol, zeta.tol, subdivisions, verbose))
  }
  if (log)
    ans - log(gamma)
  else
    ans/gamma
} ## {dstable}

dstable.quick<-function (x, alpha, beta, gamma = 1, delta = 0, pm = 0, log = FALSE,
                   tol = 64 * .Machine$double.eps, zeta.tol = NULL, subdivisions = 1000) {
  ## Like dstable but with approximations to speed calculation of large sized x
  ## Description:
  ##	 Returns density for stable DF

  ## Arguments:
  ##	 alpha = index of stability, in the range (0,2]
  ##	 beta  = skewness, in the range [-1, 1]
  ##	 gamma = scale, in the range (0, infinity)
  ##	 delta = location, in the range (-infinity, +infinity)
  ##	 pm = type of parmeterization
  ##   lower.tail = TRUE for lower tail, FALSE for upper tail
  ##   log_flag = return log of density

  ## Parameter Check:
  ## NB: (gamma, delta) can be *vector*s (vectorized along x)
  stopifnot(0 < alpha, alpha <= 2, length(alpha) == 1, -1 <=
              beta, beta <= 1, length(beta) == 1, 0 <= gamma, length(pm) ==
              1, pm %in% 0:2, tol > 0, subdivisions > 0)
  ## not an official argument {no doc!}:
  verbose <- getOption("dstable.debug", default = FALSE)
  ## Parameterizations:
  if (pm == 1) {
    delta <- delta + beta * gamma * .om(gamma, alpha)
  }
  else if (pm == 2) {
    delta <- delta - alpha^(-1/alpha) * gamma * stableMode(alpha,
                                                           beta)
    gamma <- alpha^(-1/alpha) * gamma
  } ## else pm == 0

  ## Shift and Scale:
  x <- (x - delta)/gamma
  ## Special Cases:
  if (alpha == 2) {
    ans<-dnorm(x, mean = 0, sd = sqrt(2), log = log)
  }
  else if (alpha == 1 && beta == 0) {
    ans<-dcauchy(x, log = log)
  }
  else {
    if (is.null(zeta.tol))
      zeta.tol<-0.
    ans<-as.vector(sdstable_quick(x,alpha, beta, log,
                            tol, zeta.tol, subdivisions, verbose))
  }
  if (log)
    ans - log(gamma)
  else
    ans/gamma
} ## {dstable}

### ------------------------------------------------------------------------------


pstable <- function(q, alpha, beta, gamma = 1, delta = 0, pm = 0,
                    lower.tail = TRUE, log.p = FALSE,
                    tol = 64*.Machine$double.eps, subdivisions = 1000)
{
  ## A function implemented by Diethelm Wuertz

  ## Description:
  ##	 Returns probability for stable DF

  x <- q
  ## Parameter Check:
  ## NB: (gamma, delta) can be *vector*s (vectorized along x)
  stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
             -1 <= beta, beta	<= 1, length(beta) == 1,
             0 <= gamma, length(pm) == 1, pm %in% 0:2,
             tol > 0, subdivisions > 0)
  ## not an official argument {no doc!}:
  verbose <- getOption("pstable.debug", default = FALSE)
  ## Parameterizations:
  if (pm == 1) {
    delta <- delta + beta*gamma * .om(gamma,alpha)
  } else if (pm == 2) {
    delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
    gamma <- alpha^(-1/alpha) * gamma
  } ## else pm == 0

  ## Shift and Scale:
  x <- (x - delta) / gamma

  ## Return directly
  ## ------  first, special cases:
  if (alpha == 2) {
    pnorm(x, mean = 0, sd = sqrt(2), lower.tail=lower.tail, log.p=log.p)
  } else if (alpha == 1 && beta == 0) {
    pcauchy(x, lower.tail=lower.tail, log.p=log.p)
  } else {
    as.vector(spstable(x,alpha, beta, lower.tail, log.p,
                  tol, subdivisions, verbose))
  }
}## {pstable}

## ------------------------------------------------------------------------------
qstable <- function(p, alpha, beta, gamma = 1, delta = 0, pm = 0,
                    lower.tail = TRUE, log.p = FALSE,
                    tol = 64*.Machine$double.eps, maxiter = 1000, trace = 0,
                    integ.tol = 64*.Machine$double.eps, subdivisions = 200)
{
  ## A function implemented by Diethelm Wuertz

  ## Description:
  ##	 Returns quantiles for stable DF

  ## Parameter Check:
  ## NB: (gamma, delta) can be *vector*s (vectorized along x)
  stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
             -1 <= beta, beta	<= 1, length(beta) == 1,
             0 <= gamma, length(pm) == 1, pm %in% 0:2,
             tol > 0, subdivisions > 0)
  ## not an official argument {no doc!}:
  verbose <- getOption("qstable.debug", default = FALSE)

  ## Parameterizations:
  if (pm == 1) {
    delta <- delta + beta*gamma * .om(gamma,alpha)
  } else if (pm == 2) {
    delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
    gamma <- alpha^(-1/alpha) * gamma
  } ## else pm == 0

  ## Special Cases:
  if (alpha == 2)
    result<-qnorm(p, mean = 0, sd = sqrt(2), lower.tail=lower.tail, log.p=log.p)
  else if (alpha == 1 && beta == 0)
    result<-qcauchy(p, lower.tail=lower.tail, log.p=log.p)
  ## General Case
  else { ## -------------- 0 < alpha < 2 ---------------
    ## Calculate:
    result<-as.vector(sqstable(p, alpha, beta, lower.tail, log.p,
                   tol, integ.tol, subdivisions, verbose))
  }

  ## Result:
  result * gamma + delta
}

## ------------------------------------------------------------------------------

rstable <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 0)
{
  ## Description:
  ##	 Returns random variates for stable DF

  ## slightly amended along  copula::rstable1

  ## Parameter Check:
  ## NB: (gamma, delta) can be *vector*s (vectorized along x)
  stopifnot( 0 < alpha, alpha <= 2, length(alpha) == 1,
             -1 <= beta, beta	<= 1, length(beta) == 1,
             0 <= gamma, length(pm) == 1, pm %in% 0:2)

  ## Parameterizations:
  if (pm == 1) {
    delta <- delta + beta*gamma * .om(gamma,alpha)
  } else if (pm == 2) {
    delta <- delta - alpha^(-1/alpha)*gamma*stableMode(alpha, beta)
    gamma <- alpha^(-1/alpha) * gamma
  } ## else pm == 0

  ## Calculate uniform and exponential distributed random numbers:
  u1<-runif(n)
  u2<-runif(n)
  result<-srstable(alpha,beta,u1,u2)
  ## Result:
  result * gamma + delta
}

##' dstable() for very small alpha > 0
##' ok only for  x > zeta := - beta * tan(pi/2 *alpha)
dstable.smallA <- function(x, alpha, beta, log=FALSE) {
  r <- log(alpha) + log1p(beta) - (1 + log(2*x + pi*alpha*beta))
  if(log) r else exp(r)
}

## 1e-17: seems "good", but not "optimized" at all -- hidden for now
.alpha.small.dstable <- 1e-17

##' @title C_alpha - the tail constant
##' @param alpha numeric vector of stable tail parameters, in [0,2]
##' @return
##' @author Martin Maechler
C.stable.tail <- function(alpha, log = FALSE) {
  stopifnot(0 <= alpha, alpha <= 2)
  r <- alpha
  i0 <- alpha == 0
  r[i0] <- if(log) -log(2) else 0.5
  al <- alpha[!i0]
  r[!i0] <-
    if(log) lgamma(al)-log(pi)+ log(sin(al*pi2))
  else gamma(al)/pi * sin(al*pi2)
  if(any(a2 <- alpha == 2)) r[a2] <- if(log) -Inf else 0
  r
}

##' According to Nolan's  "tail.pdf" paper, where he takes *derivatives*
##' of the tail approximation 1-F(x) ~ (1+b) C_a x^{-a}  to prove
##' that    f(x) ~  a(1+b) C_a x^{-(1+a)} ...
##'
##' @title tail approximation density for dstable()
##' @param x
##' @param alpha
##' @param beta
##' @param log if true, return  log(f(.))
##' @return
##' @author Martin Maechler
dPareto <- function(x, alpha, beta, log = FALSE) {
  if(any(neg <- x < 0)) { ## left tail
    x   [neg] <- -x	  [neg]
    beta <- rep(beta, length.out=length(x))
    beta[neg] <- -beta[neg]
  }
  if(log)
    log(alpha)+ log1p(beta)+ C.stable.tail(alpha, log=TRUE) -(1+alpha)*log(x)
  else
    alpha*(1+beta)* C.stable.tail(alpha)* x^(-(1+alpha))
}

pPareto <- function(x, alpha, beta, lower.tail = TRUE, log.p = FALSE) {
  neg <- x < 0  ## left tail
  beta <- ifelse(neg,-beta, beta)
  if(log.p) {
    ifelse(lower.tail != neg,
      log1p((1+beta)* C.stable.tail(alpha)* abs(x)^(-alpha)),
      log1p(beta)+ C.stable.tail(alpha, log=TRUE) - alpha*log(abs(x)))
  } else {
    iF <- (1+beta)* C.stable.tail(alpha)* abs(x)^(-alpha)
    ifelse(lower.tail != neg, 1-iF, iF)
  }
}

stable_fit<-function(y,type="q",quick=T,pm=0) {
  stopifnot(length(pm) == 1, pm %in% 0:2)
  results<-stable_fit_cpp(y, type, quick);
  ## Parameterizations:
  if (pm == 1) {
    for (i in 1:nrow(results)) {
      om <- .om(results$gamma[i],results$alpha[i])
    }
    results$delta <- with(results, delta - beta * gamma * om)
    results$pm <- rep(1,nrow(results))
  }
  else if (pm == 2) {
    for (i in 1:nrow(results)) {
      mode <- stableMode(results$alpha[i],results$beta[i])
    }
      results$delta <- with(results,delta + gamma * mode)
      results$gamma <- with(results, alpha^(1/alpha) * gamma)
      results$pm <- rep(2, nrow(results))
  }
  ## else pm == 0
  results
}


