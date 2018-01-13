z_error_integrand <- function(t, n, alpha, beta, x) {
  arg <- 1i*t*x
  term <- 1
  s <- rep(1,length(t))
  if (n>1) {
    for (k in 1:(n-1)) {
      term <- term * arg / k
      s <- s + term
    }
  }
  Re((exp(1i*t*x)-s)*exp(-t^alpha*exp(1i*pi/2*alpha*beta)))
}

z_error<- function(n, alpha, beta, x) {
  if (alpha < 1) {
    tmp<-integrate(z_error_integrand, 0, Inf, n, alpha, beta, x, stop.on.error=F)
    abs(tmp$value)/pi
  } else {
    x^(-1-alpha) * z_error(n, 1/alpha, 1-(2-alpha)*(1+beta), x^(-alpha))
  }
}

z_error_comp <- function(alpha, beta, x) {
  theta0<-atan(beta*tan(pi/2*alpha))/alpha
  k_alpha <- ifelse(alpha>1, alpha-2, alpha)
  betaB <- alpha * theta0/(pi/2*k_alpha)
  gammaB <- cos(alpha*theta0)^(-1/alpha)
  xB <- x/gammaB
  d_stable <- stablecpp::dstable(x, alpha, beta, pm=1)
  cat(c("theta0 = ",format(theta0), "\n"))
  cat(c("k_alpha = ", format(k_alpha), "\n"))
  cat(c("betaB = ", format(betaB), "\n"))
  cat(c("gammaB = ", format(gammaB), "\n"))
  cat(c("xB = ", format(xB), "\n"))
  cat(c("dstable = ",format(d_stable),"\n"))
  if (alpha > 1) {
    alpha_star = 1/alpha
    cat(c("alpha* = ", format(alpha_star), "\n"))
    betaB_star = 1-(2-alpha)*(1+betaB)
    cat(c("betaB* = ", format(betaB_star),"\n"))
    cat(c("xB* = ", format(xB^-alpha),"\n"))
    cat(c("cos_ab = ", format(cos(pi/2*alpha_star*betaB_star)^(1/alpha)),"\n"))
  }
  z_estimate <- rep(0,0)
  error <- rep(0,0)
  error2 <- rep(0,0)
  error_bound <- rep(0,0)
  for (n in 1:50) {
    z_estimate <- c(z_estimate, z_2_5_1(n, alpha, betaB, xB)/gammaB)
    error <- c(error,z_error(n, alpha, betaB, xB)/gammaB)
    error2 <- c(error2, z_estimate[n]-d_stable)
    error_bound <- c(error_bound,z_error_bound(n, alpha, betaB, xB)/gammaB)
  }
  data.frame(n=1:length(error), z_estimate=z_estimate, error=error,
             error2=error2, error_bound=error_bound, ratio=error/error_bound)
}

z_error_bound <- function(n, alpha, beta, x) {
  if (alpha<1) {
#    cos_ab<-cos(pi/2*alpha*beta)^(1/alpha)
    cos_ab <- 1
    1/(alpha*pi)*gamma((n+1)/alpha)/gamma(n+1)*cos_ab^(-(n+1))*x^n
  } else
    x^(-1-alpha)*z_error_bound(n, 1/alpha, 1-(2-alpha)*(1+beta), x^(-alpha))
}

z_2_5_1 <- function(n, alpha, beta, x) {
  if (alpha < 1) {
    k <- 0:(n-1)
    (1/(pi*alpha))*sum(gamma((k+1)/alpha)/gamma(k+1)*sin(pi/2*(k+1)*(1-beta))*x^k)
  } else {
    x^(-1-alpha) * z_2_5_1(n, 1/alpha, 1-(2-alpha)*(1+beta), x^(-alpha))
  }
}
