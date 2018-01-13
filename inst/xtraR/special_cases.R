dstable_holtsmark<-function(x){
  # dstable for alpha=3/2 and beta = 0
  # x[4*x^6/3^6>1]<-NA
  (1/pi)*gamma(5/3)*genhypergeo(U=c(5/12,11/12),L=c(1/3,1/2,5/6),z=-2^2*x^6/3^6)-
    x^2/(3*pi)*genhypergeo(U=c(3/4,1,5/4),L=c(2/3,5/6,7/6,4/3),z=-2^2*x^6/3^6)+
    7*x^4/(3^4*pi)*gamma(4/3)*genhypergeo(U=c(13/12,19/12),L=c(7/6,3/2,5/3),z=-2^2*x^6/3^6,debug=TRUE)
}

dstable_taleb<-function(x){
  require(Bessel)
  ## dstable for alpha=3/2 and beta=1
  ## From Taleb with slight modification (x=-x)
  arg<-x^2/(3^(4/3)*2^(2/3))
  ret<- -2^(1/3)*exp((x)^3/27)*(3^(1/3)*(x)*AiryA(arg,deriv=0)+
                           3*2^(1/3)*AiryA(arg,deriv=1))/3^(5/3)
  ret
}

dstable_fresnel<-function(x){
  require(pracma)
  # dstable for alpha=1/2 and beta=0
  abs(x)^(-3/2)/sqrt(2*pi)*(sin(1/(4*abs(x)))*(1/2-fresnelS(sqrt(1/(2*pi*abs(x)))))
                          +cos(1/(4*abs(x)))*(1/2-fresnelC(sqrt(1/(2*pi*abs(x))))))
}

dstable_bessel<-function(x){
  #dstable for alpha=1/3 and beta=1
  ret<-rep(0,length(x))
  ret[x>0]<-(1/pi)*2^(3/2)/3^(7/4)*x[x>0]^(-3/2)*
              besselK(2^(5/2)/3^(9/4)*x[x>0]^(-1/2),nu=1/3)
  ret
}

lommelS<-function(x,mu,nu){
  require(hypergeo)
  ret<-genhypergeo(U=1,L=c(.5*(mu-nu+3),.5*(mu+nu+3)),-x^2/4)*x^(mu+1)/((mu+1)^2-nu^2)+
       2^(mu+nu-1)*gamma(nu)*gamma(.5*(mu+nu+1))/(gamma(.5*(-mu+nu+1)))*genhypergeo(U=NULL,L=1-nu,-x^2/4)*x^-nu+
       2^(mu-nu-1)*gamma(.5*(mu-nu+1))*gamma(-nu)/(gamma(.5*(-mu-nu+1)))*genhypergeo(U=NULL,L=1+nu,-x^2/4)*x^nu
  ret
}

dstable_lommel<-function(x){
  # dstable for alpha=1/3 and beta-0
  ret<-rep(NA,length(x))
  sel<-x!=0
  x<-x[sel]
  ret[sel]<-2*exp(-1i*pi/4)/(3^(3/2)*pi)*abs(x)^-(3/2)*lommelS(2*exp(1i*pi/4)/3^(3/2)*abs(x)^-(1/2),mu=0,nu=1/3)
  Re(ret)
}

dstable_garoni_frankel<-function(x){
  # dstable for alpha=4/3 and beta=0
  require(hypergeo)
  ret<-rep(NA,length(x))
  sel<-abs(3^3*x^4/2^8)<1
  x<-x[sel]
  ret[sel]<-3^(5/4)/(2^(5/2)*sqrt(pi))*gamma(7/12)*gamma(11/12)/(gamma(1/2)*gamma(3/4))*
        genhypergeo(U=c(7/12,11/12),L=c(1/2,3/4),3^3*x^4/2^8)-
       3^(11/4)*x^2/(2^(13/2)*sqrt(pi))*gamma(13/12)*gamma(17/12)/(gamma(3/2)*gamma(5/4))*
        genhypergeo(U=c(13/12,17/12),L=c(3/2,5/4),3^3*x^4/2^8)
  ret
}

dstable_whittaker1<-function(x){
  require(fAsianOptions)
  # dstable for alpha=2/3 and beta =0
  ret<-rep(NA,length(x))
  sel<-abs(4/27*x^-2) < pi
  x<-x[sel]
  ret[sel]<-1/(2*sqrt(3*pi))*abs(x)^-1*exp(2/27*x^-2)*whittakerW(4/27*x^-2,-.5,1/6)
  Re(ret)
}

dstable_whittaker2<-function(x){
  require(fAsianOptions)
  # dstable for alpha=2/3 and beta =1
  ret<-rep(NA,length(x))
  ret[x<=0]<-0
  sel<-x>0 & abs(32/27*x^-2) < pi
  x<-x[sel]
  ret[sel]<-sqrt(3/pi)*abs(x)^-1*exp(-16/27*x^-2)*whittakerW(32/27*x^-2,.5,1/6)
  Re(ret)
}

dstable_zolotarev<- function(x0, alpha, beta) {
  require(Rmpfr)
  # pdf using zolotarev expansion
  # alpha and beta are Nolan's pm = 0'
  prec<-getPrec(x0)
  print(prec)
  pi.<-Const("pi",prec)
  eps<-mpfr(2,prec)^(1-prec)
  if (alpha!=1) {
    zeta <- -beta*tan(alpha*pi./2)
    if (x0<zeta) {
      beta <- -beta
      zeta <- -zeta
      x0 <- -x0
    }
    x <- x0 - zeta
    theta0 <- atan(-zeta)/alpha
    theta <- theta0/(pi./2)
    rho <- (1+theta)/2
    gammaB <- 1/cos(alpha*theta0)^(1/alpha) # = Zolotarev lambda^(1/alpha)
  } else { # alpha = 1
    # Zolotarev 1.1.8
    betaB<-beta
    gammaB<-2/pi.
    x<-x0
  }
  xB <- x/gammaB
  if (alpha < 1) {
    # Zolotarev 2.4.8
    series <- 0
    abs_series <-0
    for (n in 1:1000000) {
      term <-(-1)^(n-1)*exp(lgamma(n*alpha+1)-lgamma(mpfr(n+1,prec))) * sin(pi.*n*rho*alpha)*xB^(-n*alpha-1)/(pi.*gammaB)
      series <- series + term
      abs_series <- abs_series + abs(term)
      if (!is.finite(term) || abs(term) < 2*eps*abs(series)) break
    }
    series_error <- abs(term) + abs_series*eps
    print(c(n,c(series,series_error)))
    k_alpha<-alpha
    betaB <- alpha * theta / k_alpha
    # Zolotarev Theorem 2.5.6 with Zolotarev mu = -zeta
    x_adj <- x
    while (abs(x_error<-x_adj-zeta*x_adj^(1-alpha)-x) > 2*eps*abs(x)) {
      x_adj <- x_adj - x_error/(1-(1-alpha)*zeta*x^(-alpha))
    }
    series2<-mpfr(0,prec)
    for (n in 1:50) {
      term<-A(n, alpha, zeta, prec) * x_adj^(-n*alpha-1)/pi.
      if (n > 1 && abs(term) > abs(old_term)) break
      series2<-series2+term
      old_term<-term
      if (abs(term) < eps*abs(series2)) break
    }
    print(paste("series-series2",format(series-series2)))
    # Zolotarev 2.5.
    term <- gamma(1/alpha)/gamma(1)*sin(pi./2*(1-betaB))/(pi.*alpha*gammaB)
    tail <- term
    old_term <- term
    for (k in 1:100) {
      term <- exp(lgamma((k+1)/alpha)-lgamma(k+1))*sin(pi./2*(k+1)*(1-betaB))*xB^k/(pi.*alpha*gammaB)
      if (abs(term) > abs(old_term)) break
      tail <- tail + term
      old_term <- term
    }
    tail_error <- abs(term)
    print(c(k,as.double(tail_error)))
    if (series_error < tail_error) {
      estimated_error<<-series_error
      ret <- series
    } else {
      estimated_error<<-tail_error
      ret<- tail
    }
  } else if (alpha == 1) {
    # Zolotarev 2.4.7
    series_error<-1e300
    # Zolotarev 2.5.23
    log_x<-log(xB)
    tail<-mpfr(0,prec)
    num_small_terms<-0
    for (n in 1:50) {
      fac<-mpfr(0,prec)
      for (l in 0:n){
        m<-l:n
        term0<-chooseMpfr(mpfr(n,prec),m)*chooseMpfr(mpfr(m,prec),l)
        term0<-term0 * (1-2*((m-l)%%2)) * gamma_at_integers[m-l+1,1+n]
        term0<-term0 * betaB^m * (pi./2*(1+betaB))^(n-m) * sin(pi./2*(n-m))
        r_l_n<-sum(term0)
#        print(c(l,n,as.double(r_l_n)))
        fac<-fac+r_l_n*(log_x)^l
      }
      term<-fac * xB^(-n-1) / (factorialMpfr(n, precBits=prec)*pi.*gammaB)
      print(c(n,as.double(term)))
      if (n>1 && abs(term) > abs(old_term)) break
      tail<-tail+term
      if (abs(term) < abs(tail) * eps) {
        if (num_small_terms > 0) break
        num_small_terms<-num_small_terms+1
      } else {
        num_small_terms<-0
        old_term<-term
      }
    }
    tail_error<-max(abs(term),abs(tail) * eps)
    # Zolotarev Theorem 2.5.5

    x_adj<-xB
    while (abs(x_error<-x_adj+betaB*log(x_adj)-xB)> 2*eps*abs(xB)) {
      x_adj <- x_adj - x_error/(1-betaB/x_adj)
    }
    tail2<-mpfr(0,prec)
    num_small_terms<-0
    for (n in 1:50) {
      m<-0:floor((n-1)/2)
      term0<-(-1)^(n+m-1) * chooseMpfr(mpfr(n,prec),2*m+1)
      term0<-term0 * (pi./2*(1+betaB))^(2*m+1) * betaB^(n-2*m-1)
      term0<-term0 * gamma_at_integers[n-2*m,1+n]
      d_n <- sum(term0) / factorialMpfr(n, precBits=prec)
      term <- d_n * x_adj^(-n-1)/(pi.* gammaB)
      if (n>1 && abs(term) > abs(old_term)) break
      tail2<-tail2+term
      if (abs(term) < abs(tail2) * eps) {
        if (num_small_terms > 0) break
        num_small_terms<-num_small_terms+1
      } else {
        num_small_terms<-0
        old_term<-term
      }
    } # for n
    print(c(tail, tail2, 1-tail/tail2))
    if (series_error < tail_error) {
      ret<-series
      estimated_error<<-series_error
    } else {
      ret<-tail
      estimated_error<<-tail_error
    }
  } else {
    # alpha > 1
    # Zolotarev 2.4.6
    k_alpha <- alpha - 2
    betaB <-alpha * theta/ k_alpha
    series<-0
    abs_series <-0
    for (n in 1:1000000) {
      term<- (-1)^(n-1)*exp(lgamma(n/alpha+1) - lgamma(n+1)) * sin(pi.*n*rho) * xB ^ (n-1) /pi./gammaB
      series <- series + term
      abs_series <- abs_series + abs(term)
      if (!is.finite(term) || abs(term) < 2*eps * abs(series) ) break
    }
    series_error <- abs(term) + abs_series*eps
    print(c(n,series_error))
    # Zolotarev 2.5.4
    term<-(alpha/(gammaB*pi.))*gamma(alpha)/gamma(1)*sin(pi./2*(2-alpha)*(1+betaB))*xB^(-alpha-1)
    tail<-term
    old_term<-term
    for (k in 2:100) {
      term<-(alpha/(gammaB*pi.))*exp(lgamma(k*alpha)-lgamma(k))*sin(pi.*k/2*(2-alpha)*(1+betaB))*xB^(-alpha*k-1)
      if (abs(term) > abs(old_term)) break
      tail<-tail+term
    }
    tail_error<-abs(term)
    print(c(k,tail_error))
    if (series_error < tail_error) {
      ret <- series
      estimated_error <- series_error
    } else {
      ret <- tail
      estimated_error <- tail_error
    }
  }
  ret
}

A<-function(n, alpha, zeta, prec) {
  pi.<-Const("pi", prec)
  A_n<-mpfr(0,prec)
  for (k in 1:n) {
    m<-0:k
    term0<-sum(chooseMpfr(mpfr(k,prec),m)*zeta^m*sin(pi./2*(m-alpha*k)))*(-1)^k
    term0<-term0*gamma(alpha*k+n-k+1)
    term0<-term0/gamma(mpfr(k+1,prec))
    term0<-term0/gamma(mpfr(n-k+1,prec))
    term0<-term0*zeta^(n-k)
    A_n <- A_n + term0
  }
  A_n
}

gamma_derivative_at_integers <-function(nn) {
  require(Rmpfr)
  n <- seqMpfr(mpfr(1,166), nn)
  # First do x = 1
  # http://dlmf.nist.gov/5.4.E12 http://dlmf.nist.gov/5.15.E2
  polygamma <- c(-Const("gamma",166), (2*(n%%2)-1)*gamma(n+1)*zeta(n+1))
  n<-c(mpfr(0,166),n)
  dim(polygamma)<-c(nn+1,1)
  # http://dlmf.nist.gov/5.15.E5
  for (i in 2:(nn+1)) {
    polygamma<-cbind(polygamma,polygamma[,i-1]+(1-2*(n%%2))*factorial(n)*(i-1)^(-n-1))
  }
  n<-c(mpfr(0,166),n+1)
  lgamma_taylor<-rbind(lgamma(seqMpfr(mpfr(1,166),nn+1)),polygamma)/matrix(factorial(n),nn+2,nn+1)
  gamma_taylor<-t(array(mpfr(0,166),dim=c(nn+1,nn+2)))
  gamma_taylor[1,] <- exp(lgamma_taylor[1,])
  # j is power of x in lgamma_taylor its coefficients in row j+1 of lgamma_taylor
  for (j in 1:(nn+1)) {
    tmp <- gamma_taylor # pick up the zero order term(1) in exp
    # k i spower of x in gamma_taylor its coefficents in row k+1 of gamma_taylor
    for (k in 1:(nn+1)) {
      for (i in 1:(nn+1)){
        kk <- k - i*j
        if (kk < 0) break
        tmp[k+1,] <- tmp[k+1,] + gamma_taylor[kk+1,] * lgamma_taylor[j+1,]^i/factorial(i)
      }
    }
    gamma_taylor <- tmp
  }
  gamma_taylor*matrix(factorial(n),nn+2,nn+1)
}

#gamma_at_integers<-gamma_derivative_at_integers(50)

zolotarev_I2_test <- function(x, beta) {
  A<-x^.75
  t<-(1:1000)/10 -1i * A/x
  df <-data.frame(t=t)
  df$t3 <- exp(-1i * beta * t * log(t))
  df$abs_t3 <- abs(df$t3)
  df$d<-Im(t*log(t))
  df$d2<-Re(t)*Im(log(t))
  df$d3<--.5*x^-.25*log(Re(t)^2+x^-.5)
  df$d_bnd<-x^-.25+.5*x^-.25*log(Re(t)^2+x^-.5)
  df$d_bnd2<-1+.5*log(Re(t)^2+1)
  df
}

plot_tail<-function (x, k) {
  require(ggplot2)
  require(plyr)
  df0<-data.frame(k=k)
  f<-function(df) {
    v<-x^.75+1:200/10
    data.frame(v=v,f=(v)^df$k * (pi + abs(log(v)))^df$k * exp(-v))
  }
  f1<-function(df) {
    v<-x^.75+1:200/10
    data.frame(v=v, f=(v)^df$k * (pi + log(v))^df$k, e=exp(v/2))
  }
  df<-ddply(df0,.(k),f)
  df1 <<- ddply(df0,.(k),f1)
  qplot(x=v, y=f/exp(-v/2), data=df, color=as.factor(k), geom="line",size=I(1))

}

integrand_for_bn<-function(uu,betaB,n) {
  u<-tan(uu)
  (1/gamma(n+1))*exp(-betaB*u*log(u))*(u^(n-1))*sin((1+betaB)*u*pi/2)/cos(uu)^2
}

Q_integrand<-function(u, beta, x){
  ret<-rep(NA, length(u))
  eps<-50*.Machine$double.eps
  sel<-(u>eps)
  u<-u[sel]
  ret[sel]<-exp(-x*u-beta*(2/pi)*u*log(u))*sin((1+beta)*u)/(u*pi)
  ret[!sel]<-(1+beta)/2
  ret

}

q_integrand<-function(u, beta, x) {
  ret<-rep(NA, length(u))
  eps<-50*.Machine$double.eps
  sel<-u>eps
  u<-u[sel]
  ret[sel]<-  exp(-x*u-beta*(2/pi)*u*log(u))*sin(u*(1+beta))/(pi)
  ret[!sel]<-0
  ret
}


