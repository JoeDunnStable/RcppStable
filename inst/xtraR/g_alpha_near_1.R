require(plyr)
require(ggplot2)

eps_diff<-function(l,r){
  stopifnot(length(r)==length(l))
  out <- r  #just for size
  sel_eq <- r==l
  out[sel_eq]<-0
  sel_inf <-!sel_eq & (abs(r)==Inf | abs(l)==Inf)
  out[sel_inf]<-Inf
  sel<-!(sel_eq | sel_inf)
  r<-r[sel]
  l<-l[sel]
  out[sel]<- abs(r-l)/pmax(1, abs(l),abs(r))/.Machine$double.eps
  round(out,1)
}

g_far_1 <- function(x, alpha, beta) {
  # S0 parameterization
  df0<-expand.grid(x=x, alpha=alpha, beta=beta)
  g_df_1<- function(df0) {
    x<-df0$x
    alpha <- df0$alpha
    alpha_m_1 <- alpha - 1
    beta <- df0$beta
    pi2<-pi/2
    tol<-8*.Machine$double.eps

    if (alpha_m_1 != 0) {
      zeta <- -beta*tan(pi2*alpha)
      x_m_zet <- abs(x-zeta)
      if (x < zeta) {
        x <- -x
        beta <- -beta
        zeta <- -zeta
      }
      theta0<-atan(-zeta)/alpha
      cat0<-1/sqrt(1+zeta^2)
      th_min<-0;
      th_max<-(pi2)+theta0;
      add_l<-theta0-pi2;
      add_r<-max(0,pi-alpha*(theta0+(pi2)));
      if (alpha<1) {
        if (beta==1){
          add_l<-0
        } else if (beta==-1){
          add_l<-pi
        }
      } else if (alpha>1) {
        if (beta==-1){
          add_r<-0;
        }
      }
      if (x_m_zet < max(1,abs(zeta))) {
        type <- "left"
        g <- function (th_l){
          # Similar to g except the input variable is th_l = th+theta0
          ret <- th_l  # just for size
          i_g0 <- ((alpha < 1 & th_l==0) | (alpha > 1 & abs(th_l-th_max)<tol))
          i_g_inf <- ((alpha < 1 & abs(th_l-th_max)<tol) | (alpha >1 & th_l==0))
          if ((alpha < 1 & beta==1) || (alpha > 1 & beta == -1)){
            g0 <- (cat0*(x_m_zet/alpha)^alpha)^(1/(alpha_m_1))*abs(1-alpha)
          } else{
            g0 <- 0
          }
          ret[i_g0] <- g0
          ret[i_g_inf] <- Inf
          sel <- !(i_g0 | i_g_inf)
          th_l<-th_l[sel]
          catt_m_t <- pmax(0,-sin((alpha-1)*th_l+add_l))
          costh <- pmax(0,sin(th_l-add_l))
          att <- alpha*th_l;
          if (abs(zeta) < 1 || abs(x_m_zet+zeta) > .1 * abs(zeta)) {
            x_sin <- x_m_zet/sin(att)
            pow1 <- x_sin^alpha;
            pow2 <- (cat0*costh*pow1)^(1/(alpha-1))
          } else {
            ln_pow1 <- alpha * (log(abs(zeta))+log1p(-x_m_zet/zeta-1)-log(sin(att)))
            ln_pow2 <- (log(cat0) + log(costh)+ln_pow1)/(alpha-1)
            pow2 = exp(ln_pow2);
          }
          ret[sel]<-pow2*catt_m_t;
          ret
        } # g left
      } else {
        type <- "right"
        g<-function(th_r) {
          # Similar to g except the input variable is th_r = pi/2 - th
          ret <- th_r # just for size
          i_g_inf <-((alpha>1 & abs(th_r-th_max)<tol) | (alpha<1 & th_r==0) )
          i_g0<-((alpha>1 & th_r==0) | (alpha <1 & abs(th_r-th_max)<tol))
          ret[i_g_inf] <- Inf
          if ((alpha<1 && beta == 1) | (alpha > 1 && beta == -1)){
            g0 <- (cat0*(x_m_zet/alpha)^alpha)^(1/(alpha_m_1))*abs(1-alpha)
          } else{
            g0 <- 0
          }
          ret[i_g0] <- g0;
          sel <-!(i_g0 | i_g_inf)
          th_r <- th_r[sel]
          att <- alpha*th_r+add_r
          catt_m_t <- pmax(0.,sin((alpha-1)*th_r+add_r))
          if (x_m_zet < 1e100) {
            costh <- pmax(0,sin(th_r));
            pow1 <- (x_m_zet/sin(att))^alpha
            pow2 <- (cat0 * costh * pow1)^(1/(alpha-1))
          } else {
            ln_costh <- log(sin(th_r))
            ln_pow1 <- alpha*(log(x_m_zet) - log(sin(att)));
            pow2 <- exp((log(cat0) + ln_costh + ln_pow1)/(alpha - 1));
          }
          ret[sel]<- pow2*catt_m_t;
          ret
        } # g right
      }

    } else {
      # alpha == 1
      type <- "right"
      th_min <- 0
      th_max <- pi
      if (x < 0) {
        x <- -x
        beta <- -beta
      }
      g<-function(th_r) {
        i_g0 <- (beta>0 & abs(th_r -th_max) < tol) |
          (beta<0 & th_r == 0)
        i_g_inf <- (beta > 0 & th_r == 0) |
          (beta<0 & abs(th_r - th_max) < tol)
        ln_g <- th_r
        ea <- -pi2*x/beta
        h <- pi2/beta + pi2 -th_r
        h_tanh <- h / tan(th_r)
        ln_h2b <- log(h * beta / pi2)
        ln_costh <- log(sin(th_r))
        residual <- 0
        ln_g[i_g0] <- -Inf
        ln_g[i_g_inf] <- Inf
        sel<-!(i_g0 | i_g_inf)
        ln_g[sel]<-ea+h_tanh[sel]+ln_h2b[sel]-ln_costh[sel]+residual
        exp(ln_g)
      } #g right for alpha == 1
    }
    if (abs(th_max-th_min) > 4*.Machine$double.eps) {
      g_low <- g(th_min)
      g_hi <- g(th_max)
      if (min(g_low,g_hi) < 1 && 1 < max(g_low,g_hi)) {
        g2 <- function(th) { min(g(th),.Machine$double.xmax) - 1}
        th2 <- uniroot(g2, c(th_min, th_max),tol=64*.Machine$double.eps)$root
        th <-  sort(c(th_min+(0:400)/400 * (th_max-th_min),
                      th_min+(50:100)/100*(th2-th_min),
                      th2+(1:50)/100*min(th2-th_min,th_max-th2)))
      } else {
        th<-th_min +(0:500)/500*(th_max-th_min)
      }
    data.frame(th=th, type=type, g=g(th))
    } else {
      data.frame()
    }
  } # g_df_1
  df <- ddply(df0, .(x, alpha, beta), g_df_1)
  df$exp_m_g <- exp(-df$g)
  df$g_exp_m_g <- ifelse(df$exp_m_g==0,0,df$g * df$exp_m_g)
  df$g_1_alpha_g_exp_g <- ifelse(df$exp_m_g==0,0,df$g*(1-df$alpha*df$g)*df$exp_m_g)
  df
} # g for alpha far from 1

g_mpfr <- function(x, alpha, beta, prec=120) {
  require(Rmpfr)
  # S0 parameterization
  df0<-expand.grid(x=x, alpha=alpha, beta=beta)
  g_df_1<- function(df0) {
    x<-mpfr(df0$x,prec=prec)
    alpha <- mpfr(df0$alpha,prec=prec)
    alpha_m_1 <- alpha - 1
    beta <- mpfr(df0$beta,prec=prec)
    pi2<-Const("pi", prec=prec)/2
    tol<-2^-(prec-8)

    if (alpha_m_1 != 0) {
      zeta <- -beta*tan(pi2*alpha)
      x_m_zet <- abs(x-zeta)
      if (x < zeta) {
        x <- -x
        beta <- -beta
        zeta <- -zeta
      }
      theta0<-atan(-zeta)/alpha
      cat0<-1/sqrt(1+zeta^2)
      th_min<-mpfr(0,prec=prec);
      th_max<-(pi2)+theta0;
      add_l<-theta0-pi2;
      add_r<-max(mpfr(0,prec=prec),2*pi2-alpha*(theta0+(pi2)));
      if (alpha<1) {
        if (beta==1){
          add_l<-0
        } else if (beta==-1){
          add_l<-pi
        }
      } else if (alpha>1) {
        if (beta==-1){
          add_r<-0;
        }
      }
      if (x_m_zet < max(mpfr(1,prec=prec),abs(zeta))) {
        type <- "left"
        g <- function (th_l){
          # Similar to g except the input variable is th_l = th+theta0
          ret <- th_l  # just for size
          i_g0 <- ((alpha < 1 & th_l==0) | (alpha > 1 & abs(th_l-th_max)<tol))
          i_g_inf <- ((alpha < 1 & abs(th_l-th_max)<tol) | (alpha >1 & th_l==0))
          if ((alpha < 1 & beta==1) || (alpha > 1 & beta == -1)){
            g0 <- (cat0*(x_m_zet/alpha)^alpha)^(1/(alpha_m_1))*abs(1-alpha)
          } else{
            g0 <- 0
          }
          ret[i_g0] <- g0
          ret[i_g_inf] <- Inf
          sel <- !(i_g0 | i_g_inf)
          th_l<-th_l[sel]
          catt_m_t <- pmax(0,-sin((alpha-1)*th_l+add_l))
          costh <- pmax(0,sin(th_l-add_l))
          att <- alpha*th_l;
          if (abs(zeta) < 1 || abs(x_m_zet+zeta) > .1 * abs(zeta)) {
            x_sin <- x_m_zet/sin(att)
            pow1 <- x_sin^alpha;
            pow2 <- (cat0*costh*pow1)^(1/(alpha-1))
          } else {
            ln_pow1 <- alpha * (log(abs(zeta))+log1p(-x_m_zet/zeta-1)-log(sin(att)))
            ln_pow2 <- (log(cat0) + log(costh)+ln_pow1)/(alpha-1)
            pow2 = exp(ln_pow2);
          }
          ret[sel]<-pow2*catt_m_t;
          ret
        } # g left
      } else {
        type <- "right"
        g<-function(th_r) {
          # Similar to g except the input variable is th_r = pi/2 - th
          ret <- th_r # just for size
          i_g_inf <-((alpha>1 & abs(th_r-th_max)<tol) | (alpha<1 & th_r==0) )
          i_g0<-((alpha>1 & th_r==0) | (alpha <1 & abs(th_r-th_max)<tol))
          ret[i_g_inf] <- Inf
          if ((alpha<1 && beta == 1) | (alpha > 1 && beta == -1)){
            g0 <- (cat0*(x_m_zet/alpha)^alpha)^(1/(alpha_m_1))*abs(1-alpha)
          } else{
            g0 <- 0
          }
          ret[i_g0] <- g0;
          sel <-!(i_g0 | i_g_inf)
          th_r <- th_r[sel]
          att <- alpha*th_r+add_r
          catt_m_t <- pmax(0.,sin((alpha-1)*th_r+add_r))
          if (x_m_zet < 1e100) {
            costh <- pmax(0,sin(th_r));
            pow1 <- (x_m_zet/sin(att))^alpha
            pow2 <- (cat0 * costh * pow1)^(1/(alpha-1))
          } else {
            ln_costh <- log(sin(th_r))
            ln_pow1 <- alpha*(log(x_m_zet) - log(sin(att)));
            pow2 <- exp((log(cat0) + ln_costh + ln_pow1)/(alpha - 1));
          }
          ret[sel]<- pow2*catt_m_t;
          ret
        } # g right
      }

    } else {
      # alpha == 1
      type <- "right"
      th_min <- 0
      th_max <- 2*pi2
      if (x < 0) {
        x <- -x
        beta <- -beta
      }
      g<-function(th_r) {
        i_g0 <- (beta>0 & abs(th_r -th_max) < tol) |
          (beta<0 & th_r == 0)
        i_g_inf <- (beta > 0 & th_r == 0) |
          (beta<0 & abs(th_r - th_max) < tol)
        ln_g <- th_r
        ea <- -pi2*x/beta
        h <- pi2/beta + pi2 -th_r
        h_tanh <- h / tan(th_r)
        ln_h2b <- log(h * beta / pi2)
        ln_costh <- log(sin(th_r))
        residual <- 0
        ln_g[i_g0] <- -Inf
        ln_g[i_g_inf] <- Inf
        sel<-!(i_g0 | i_g_inf)
        ln_g[sel]<-ea+h_tanh[sel]+ln_h2b[sel]-ln_costh[sel]+residual
        exp(ln_g)
      } #g right for alpha == 1
    }
    if (abs(th_max-th_min) > tol) {
      g_low <- g(th_min)
      g_hi <- g(th_max)
      if (min(g_low,g_hi) < 1 && 1 < max(g_low,g_hi)) {
        g2 <- function(th) {
          g_raw<-g(th)
          # print(c(th,g_raw))
          min(g(th)-1,mpfr(.Machine$double.xmax,prec=prec))
          }
        th2 <- unirootR(g2, c(th_min, th_max),tol=2^-(prec-8))$root
        th <-  sort(c(th_min+(0:400)/400 * (th_max-th_min),
                      th_min+(50:100)/100*(th2-th_min),
                      th2+(1:50)/100*min(th2-th_min,th_max-th2)))
      } else {
        th<-th_min +(0:500)/500*(th_max-th_min)
      }
      data.frame(th=as.double(th), type=type, g=as.double(g(th)))
    } else {
      data.frame()
    }
  } # g_df_1
  df <- ddply(df0, .(x, alpha, beta), g_df_1)
  df$exp_m_g <- exp(-df$g)
  df$g_exp_m_g <- ifelse(df$exp_m_g==0,0,df$g * df$exp_m_g)
  df$g_1_alpha_g_exp_g <- ifelse(df$exp_m_g==0,0,df$g*(1-df$alpha*df$g)*df$exp_m_g)
  df
} # g for alpha far from 1

g_alpha_near_1 <- function(x, alpha_m_1, beta) {
  # S0 parameterization
  df0<-expand.grid(x=x, alpha_m_1=alpha_m_1, beta=beta)
  g_df_1 <-function(df0) {
    x<-df0$x
    alpha_m_1<-df0$alpha_m_1
    alpha <- 1 + alpha_m_1
    beta<-df0$beta
    pi2<-pi/2
    tol <- 8*.Machine$double.eps
    if (alpha_m_1 != 0) {
      zeta <- beta/tan(pi2*alpha_m_1)
      if (x < zeta) {
        x <- -x
        beta <- -beta
        zeta <- -zeta
      }
      if (sign(x) != sign(zeta)) {
        add_r <- (-alpha_m_1*pi2-atan(1/zeta))
        g <- function(th_r) {
          i_g0 <- (alpha_m_1 < 0 & abs(th_r - th_max) < tol)|
            (alpha_m_1 > 0 & abs(th_r) < tol)
          i_g_inf <- (alpha_m_1 < 0 & abs(th_r) < tol)|
            (alpha_m_1 > 0 & abs(th_r - th_max) < tol)
          ea <- (alpha/alpha_m_1) * log1p(-x/zeta)
          residual <- -1/(2*alpha_m_1)*log1p(zeta^-2)
          if ((alpha<1 && beta == 1) | (alpha > 1 && beta == -1)){
            ln_g0 <- ea + residual - alpha/alpha_m_1*log1p(alpha_m_1) +log(abs(zeta*alpha_m_1))
          } else{
            ln_g0 <- -Inf
          }
          ln_g <- th_r
          h_tanh <- -(alpha/alpha_m_1)*log1p(pmax(-1,-2 * sin((alpha_m_1*th_r+add_r)/2)^2 +
                                               1/tan(th_r)*sin(alpha_m_1*th_r + add_r)))
          ln_h2b <- log(abs(zeta)*sin(alpha_m_1*th_r+add_r))
          ln_costh <- log(sin(th_r))
          ln_g[i_g0] <- ln_g0
          ln_g[i_g_inf] <- Inf
          sel<-!(i_g0 | i_g_inf)
          ln_g[sel]<-ea+h_tanh[sel]+ln_h2b[sel]-ln_costh[sel]+residual
          data.frame(type=rep("right",length(th_r)),g=exp(ln_g), ea=ea, h_tanh=h_tanh, ln_h2b=ln_h2b,
                     ln_costh=ln_costh, residual=residual)
        } # g
      } else {
        add_l<--1/alpha * (alpha_m_1 * pi2-atan(1/beta * tan(alpha_m_1*pi2)))
        g <- function(th_l){
          i_g_inf <- (alpha_m_1 < 0 & abs(th_l - th_max) < tol)|
            (alpha_m_1 > 0 & abs(th_l) < tol)
          i_g0 <- (alpha_m_1 < 0 & abs(th_l) < tol)|
            (alpha_m_1 > 0 & abs(th_l - th_max) < tol)
          ea <- (alpha/alpha_m_1) * log1p(-x/zeta)
          residual <- -1/(2*alpha_m_1)*log1p(zeta^-2)
          if ((alpha<1 && beta == 1) | (alpha > 1 && beta == -1)){
            ln_g0 <- ea + residual - alpha/alpha_m_1*log1p(alpha_m_1) +log(abs(zeta*alpha_m_1))
          } else{
            ln_g0 <- -Inf
          }
          ln_g <- th_l
          h_tanh <- -(alpha/alpha_m_1)*log1p(pmax(-1,-2 * sin((alpha_m_1*th_l+add_l)/2)^2 +
                                               1/tan(th_l-add_l)*sin(alpha_m_1*th_l + add_l)))
          ln_h2b <- log(abs(zeta)*-sin(alpha_m_1*th_l+add_l))
          ln_costh <- log(sin(th_l-add_l))
          ln_g[i_g0] <- ln_g0
          ln_g[i_g_inf] <- Inf
          sel<-!(i_g0 | i_g_inf)
          ln_g[sel]<-ea+h_tanh[sel]+ln_h2b[sel]-ln_costh[sel]+residual
          data.frame(type=rep("left",length(th_l)),g=exp(ln_g), ea=ea, h_tanh=h_tanh, ln_h2b=ln_h2b,
                     ln_costh=ln_costh, residual=residual)
        } # g
      }
    } else { # alpha_m_1 == 0
      if (x < 0) {
        x <- -x
        beta <- -beta
      }
      g<-function(th_r) {
        i_g0 <- (beta>0 & abs(th_r -th_max) < tol) |
                (beta<0 & abs(th_r) < tol)
        i_g_inf <- (beta > 0 & abs(th_r) < tol) |
                (beta<0 & abs(th_r - th_max) < tol)
        ln_g <- th_r
        ea <- -pi2*x/beta
        h <- pi2/beta + pi2 -th_r
        h_tanh <- h / tan(th_r)
        ln_h2b <- log(h * beta / pi2)
        ln_costh <- log(sin(th_r))
        residual <- 0
        ln_g[i_g0] <- -Inf
        ln_g[i_g_inf] <- Inf
        sel<-!(i_g0 | i_g_inf)
        ln_g[sel]<-ea+h_tanh[sel]+ln_h2b[sel]-ln_costh[sel]+residual
        data.frame(type=rep("right",length(th_r)),g=exp(ln_g), ea=ea, h_tanh=h_tanh, ln_h2b=ln_h2b,
                     ln_costh=ln_costh, residual=residual)
      }
    }
    if (alpha_m_1 == 0 || abs(x) < abs(zeta)) {
      th_min<-0
      th_max<-pi+(-alpha_m_1*pi2+atan(tan(alpha_m_1*pi2)/beta))/alpha
      g2 <- function(th) { min(g(th)$g,.Machine$double.xmax) - 1}
      th2 <- uniroot(g2, c(0, th_max),tol=64*.Machine$double.eps)$root
      th <-  sort(c(th_min+(0:400)/400 * (th_max-th_min),
                    th_min+(50:100)/100*(th2-th_min),
                    th2+(1:50)/100*min(th2-th_min,th_max-th2)))
      df<- cbind(data.frame(th=th), g(th))
      df$exp_m_g <- exp(-df$g)
      df$g_exp_m_g <- ifelse(df$exp_m_g==0,0,df$g * df$exp_m_g)
    } else {
      df = data.frame()
    }
    df
  } #g_df_1
  ddply(df0, .(x, alpha_m_1, beta), g_df_1)
}

check_alpha_near_1 <- function(alpha_m_1, beta) {
  require(Rmpfr)
  df<-expand.grid(alpha_m_1=alpha_m_1, beta=beta)
  check_one<-function(alpha_m_1, beta) {
    # Check alternate formulae near alpha = 1
    # mpfr default precision is 120 bits, or about 36 decimal places
    near<-abs(alpha_m_1)<.5*abs(beta)
    prec=mpfr_default_prec(120)
    pi2<-pi/2
    big_pi2 <- Const("pi")/2
    alpha <- 1+alpha_m_1
    big_alpha_m_1<-mpfr(alpha_m_1,prec=prec)
    big_alpha<- 1+big_alpha_m_1
    big_beta<-mpfr(beta, prec=prec)
    big_zeta <- -big_beta * tan(big_pi2*big_alpha)
    zeta_far<--beta*tan(pi2*alpha)
    zeta_near<-beta/tan(pi2*alpha_m_1)
    zeta<-ifelse(near,zeta_near,zeta_far)
    err_zeta<-eps_diff(zeta, as.double(big_zeta))
    big_theta0 <- atan(-big_zeta)/big_alpha
    theta0_far <- atan(-zeta_far)/alpha
    theta0_near <- ifelse(zeta_near < 0,1,-1)*
                   (pi2-(pi2*alpha_m_1 + atan(1/abs(zeta_near)))/alpha)
    theta0<-ifelse(near, theta0_near, theta0_far)
    err_theta0<-eps_diff(theta0, as.double(big_theta0))
    big_th_max <- big_pi2 + big_theta0
    th_max_far<-pi2+theta0_far
    th_max_near<-ifelse(zeta_near<0,
                    pi-(pi2*alpha_m_1+atan(1/abs(zeta_near)))/alpha,
                    (pi2*alpha_m_1+atan(1/abs(zeta_near)))/alpha)
    th_max<-ifelse(near, th_max_near, th_max_far)
    err_th_max<-eps_diff(th_max, as.double(big_th_max))
    big_add_l<- big_theta0 - big_pi2
    add_l_far<-theta0_far-pi2
    add_l_far[alpha_m_1<0 & beta==1]<-0
    add_l_near<-ifelse(zeta_near<0,-(pi2*alpha_m_1 + atan(1/abs(zeta_near)))/alpha,
                   (pi2*alpha_m_1 + atan(1/abs(zeta_near)))/alpha-pi)
    add_l<-ifelse(near,add_l_near, add_l_far)
    err_add_l<-eps_diff(add_l, as.double(big_add_l))
    big_add_r<-2*big_pi2-big_alpha*(big_th_max)
    add_r_far<-pi-alpha*(th_max_far)
    add_r_far[alpha_m_1>0 & beta==-1]<-0
    add_r_near<-ifelse(zeta_near<0,
                   -alpha_m_1*pi+(pi2*alpha_m_1+atan(1/abs(zeta_near))),
                   pi-(pi2*alpha_m_1+atan(1/abs(zeta_near))))
    add_r<-ifelse(near,add_r_near, add_r_far)
    err_add_r<-eps_diff(add_r, as.double(big_add_r))

    data.frame(alpha_m_1=alpha_m_1, beta=beta,
               big_zeta=as.double(big_zeta), zeta=zeta, err_zeta=err_zeta,
               big_theta0=as.double(big_theta0), theta0=theta0,err_theta0=err_theta0,
               big_th_max=as.double(big_th_max), th_max=th_max, err_th_max=err_th_max,
               big_add_l=as.double(big_add_l), add_l=add_l, err_add_l=err_add_l,
               big_add_r=as.double(big_add_r), add_r=add_r, err_add_r=err_add_r)
  }
  check_one(df$alpha_m_1, df$beta)
}
tmp<-c(2^(-2-8*(6:0)),.5,1-2^(-2-8*(0:6)))
alpha_m_1 <- c(-rev(tmp), tmp)
beta <- (-4:4)/4
df_chk<-check_alpha_near_1(alpha_m_1, beta)

tmp<-2^(-8*(6:2))
alpha_m_1<-c(-rev(tmp),tmp)
x <- -10
beta<-.5
df_1 <- g_alpha_near_1(x, alpha_m_1, beta)
qplot(x=th, y=g_exp_m_g, data=df_1[df_1$g_exp_m_g>1e-9,],
      color = as.factor(alpha_m_1), geom="line",
      main="dstable(-10, alpha, .5)")

x <- 10
df_2 <- g_alpha_near_1(x, alpha_m_1, beta)
qplot(x=th, y=g_exp_m_g, data=df_2[df_2$g_exp_m_g>1e-9,],
      color = as.factor(alpha_m_1), geom="line",
     main="dstable(10, alpha, .5)")

alpha <- 1+1/128
beta <- 1
x <--.1

df_3 <- g_far_1(x, alpha, beta)
df_4<-g_mpfr(x,alpha,beta)
df_5<-g_alpha_near_1(x, alpha-1, beta)
err_far<-eps_diff(df_3$g,df_4$g)
err_near<-eps_diff(df_5$g,df_4$g)
list(near=sum(err_far>err_near),far=sum(err_near>err_far))

qplot(x=th, y=g_exp_m_g, data=df_3[df_3$g_exp_m_g > 1e-9, ],
      color = as.factor(alpha), geom="line",
      main=paste("Integrand for dstable(x = ", x,", alpha, beta =", beta, ")"))
