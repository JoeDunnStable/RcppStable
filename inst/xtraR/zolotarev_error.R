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
  d_stable <- RcppStable::dstable(x, alpha, beta, pm=1)
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
  n<-50
  z_asymptotic <- z_2_5_1(n, alpha, betaB, xB)/gammaB
  z_convergent <- z_2_4_6(n, alpha, betaB, xB)/gammaB
  error_asymptotic <- z_asymptotic-d_stable
  error_convergent <- z_convergent- d_stable
  data.frame(n=1:n, z_asymptotic=z_asymptotic, z_convergent=z_convergent,
             error_asymptotic=error_asymptotic, error_convergent=error_convergent)
}

z_error_bound <- function(n, alpha, beta, x) {
  if (alpha<1) {
#    cos_ab<-cos(pi/2*alpha*beta)^(1/alpha)
    cos_ab <- 1
    1/(alpha*pi)*gamma((n+1)/alpha)/gamma(n+1)*cos_ab^(-(n+1))*x^n
  } else
    x^(-1-alpha)*z_error_bound(n, 1/alpha, 1-(2-alpha)*(1+beta), x^(-alpha))
}

z_2_5_1 <- function(n, alpha, betaB, xB) {
  # asymptotic series for x positive
  if (alpha < 1) {
    # for xB small, beta != 1
    k <- 0:(n-1)
    (1/(pi*alpha))*cumsum(gamma((k+1)/alpha)/gamma(k+1)*sin(pi/2*(k+1)*(1-betaB))*xB^k)
  } else {
    # for x large, beta != -1
    xB^(-1-alpha) * z_2_5_1(n, 1/alpha, 1-(2-alpha)*(1+betaB), xB^(-alpha))
  }
}

z_2_4_6 <- function(n, alpha, betaB, xB) {
  # convergent series for positive x
  if (alpha < 1) {
    # for x large
    k<-1:n
    theta<-betaB
    rho<-(1+theta)/2
    (1/pi)*cumsum((-1)^(k-1)*gamma(k*alpha+1)/gamma(k+1)*sin(pi*k*rho*alpha)*xB^(-k*alpha-1))
  } else {
    # alpha > 1 and x small,
    xB^(-1-alpha) * z_2_4_6(n, 1/alpha, 1-(2-alpha)*(1+betaB), xB^-alpha)
  }
}
