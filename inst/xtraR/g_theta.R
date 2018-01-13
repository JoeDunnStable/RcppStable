require(plyr)

g_df<- function(x, alpha, beta) {
  df0<-expand.grid(x=x, alpha=alpha, beta=beta)
  g_df_1 <-function(df0) {
    x<-df0$x
    alpha<-df0$alpha
    beta<-df0$beta
    pi2<-pi/2
    if (alpha != 1) {
      zeta <- -beta * tan(pi2 * alpha)
      theta0 <- atan(-zeta)/alpha
      x.m.zet <- abs(x - zeta)
      if(x < zeta) theta0 <- -theta0
      
      g <- function(th) {
        a_1 <- alpha-1
        at0 <- alpha * theta0
        cat0 <- cos(alpha * theta0)
        r <- th
        ## g(-pi/2) or g(pi/2) could become  NaN --> work around
        i.bnd <- abs(pi2 -sign(a_1)*th) < 64*.Machine$double.eps
        r[i.bnd] <- 0
        th <- th[io <- !i.bnd]
        att <- at0 + alpha*th ## = alpha*(theta0 + theta)
        r[io] <- (cat0 * cos(th) * (x.m.zet/sin(att))^alpha)^(1/a_1) * cos(att-th)
        r
      } # g
    } else {
      zeta <- 0
      if (x < 0) {
        x <- -x
        beta<--beta
      }
      theta0 <- pi2
      i2b <- 1/(2*beta)
      p2b <- pi*i2b # = pi/(2 beta)
      ea <- -p2b*x
      if(is.infinite(ea)) return else 0
      
      ##' g() is strictly monotone;
      ##'  for beta > 0: increasing from g(-pi2) = 0   to  g(+pi2) = Inf
      ##'  for beta < 0: decreasing from g(-pi2) = Inf to  g(+pi2) = 0
      t0 <- -sign(beta)*pi2# g(t0) == 0  mathematically, but not always numerically
      g <- function(th) {
        r <- th
        r[i <- abs(th-t0) < 1e-10] <- 0
        th <- th[!i]
        h <- p2b+ th # == g'/beta where g' := pi/2 + beta*th = pi/2* (1 + beta*u)
        r[!i] <- (h/p2b) * exp(ea + h*tan(th)) / cos(th)
        r
      }
      
    }
    
    g2 <- function(th) { min(g(th),.Machine$double.xmax) - 1}
    
    
    th <- -theta0 + (0:500)/500 * (pi2+theta0)
    th2 <- uniroot(g2, c(-theta0, pi2), tol=64*.Machine$double.eps)$root
    th <- sort(c(th,th2))
    df<- data.frame(th=th, g=g(th))
    df$exp_m_g <- exp(-df$g)
    df$g_exp_m_g <- df$g * df$exp_m_g
    df
  } #g_df_1
  ddply(df0, .(x, alpha, beta), g_df_1)        
}

alpha <- 1.5
beta <- .5
zeta <- -beta*tan(pi/2*alpha)
x<- zeta +c(.1, 1, 10)
df_1.5<-g_df(x, alpha, beta)

require(ggplot2)
qplot(x=th, y=g_exp_m_g, data=df_1.5, geom="line", color=as.factor(x),
      main="Integrand for dstable(alpha = 1.5, beta = .5)")

alpha <- 1
zeta <- 0
x<- zeta +c(-10, -1, 0, 1, 10)
df_1<-g_df(x, alpha, beta)

require(ggplot2)
qplot(x=th, y=g_exp_m_g, data=df_1, geom="line", color=as.factor(x),
      main="Integrand for dstable(alpha = 1, beta = .5)")

alpha <- .5
zeta <- -beta*tan(pi/2*alpha)
x<- zeta +c(.01, .1, 1, 10, 100)
df_.5<-g_df(x, alpha, beta)

require(ggplot2)
qplot(x=th, y=g_exp_m_g, data=df_.5, geom="line", color=as.factor(x),
      main="Integrand for dstable(alpha = .5, beta = .5)")


alpha <- .99
beta <- 0
x<- +c( .1, 1, 10, 100)
df_pk1<-g_df(x, alpha, beta)

qplot(x=th, y=g_exp_m_g, data=df_pk1, geom="line", color=as.factor(x),
     main="Integrand for dstable(alpha = .99, beta = 0)")

alpha <- 1
beta <- .01
x<- +c( 0, 1, 10, 100)
df_pk2<-g_df(x, alpha, beta)

qplot(x=th, y=g_exp_m_g, data=df_pk2, geom="line", color=as.factor(x),
      main="Integrand for dstable(alpha = 1, beta = .01)")
