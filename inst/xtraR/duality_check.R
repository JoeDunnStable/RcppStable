require(stablecpp)
require(plyr)
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
  out
}

map_plus<-function(alpha, beta) {
  pi2<-pi/2
  alpha_prime<-1/alpha
  Q<-1-atan(beta*tan(pi2*alpha))/pi2
  if (alpha==2 || beta==-1)
    beta_star <- 1
  else {
    beta_star<-1/(tan(pi2*alpha_prime)*tan(pi2*Q*alpha_prime))
  }
  D<-(1+(beta*tan(pi2*alpha))^2)^(1/(2*alpha))/sin(pi2*Q*alpha_prime)
  data.frame(alpha_prime=alpha_prime, beta_star=beta_star, D=D)
}

map_minus<-function(alpha, beta) {
  out<-map_plus(alpha,-beta)
  out$beta_star <- -out$beta_star
  out
}

duality_check<-function(y, alpha, beta) {
  #Zolotarev Theorem 2.3.2
  stopifnot(alpha>1,-1<=beta, beta<=1)
  y<-y[y!=0]
  min_ln<-log(.Machine$double.xmin)
  df_plus <- map_plus(alpha, beta)
  df_minus <- map_minus(alpha, beta)
  ln_pdf1  <- y  #for size only
  ln_pdf2 <- y   #for size only
  alpha_prime<- y
  beta_star<- y
  sel<-y>0
  ln_pdf1[sel]<-with(df_plus,dstable(y[sel], alpha_prime, beta_star, pm=1, log=T))
  ln_pdf2[sel]<-with(df_plus,log(D) + (-1-alpha_prime)*log(y[sel]) +
              dstable(D*y[sel]^(-alpha_prime),alpha, beta, pm=1, log=T))
  alpha_prime[sel] <- with(df_plus, alpha_prime)
  beta_star[sel]<-with(df_plus, beta_star)
  sel<-y<0
  ln_pdf1[sel]<-with(df_minus,dstable(y[sel], alpha_prime, beta_star, pm=1, log=T))
  ln_pdf2[sel]<-with(df_minus,log(D) + (-1-alpha_prime)*log(abs(y[sel])) +
                       dstable(-D*abs(y[sel])^(-alpha_prime),alpha, beta, pm=1, log=T))
  alpha_prime[sel] <- with(df_minus, alpha_prime)
  beta_star[sel]<-with(df_minus, beta_star)
  ln_pdf1[ln_pdf1<min_ln]<--Inf
  ln_pdf2[ln_pdf2<min_ln]<--Inf
  data.frame(y=y,alpha_prime=alpha_prime, beta_star=beta_star,ln_pdf1=ln_pdf1,
             alpha=alpha, beta=beta, ln_pdf2=ln_pdf2, eps_diff=eps_diff(ln_pdf1,ln_pdf2))
}

alpha<-1+c(1/128,(1:12)/12)
beta<-(-16:16)/16
df_plus<-ddply(expand.grid(alpha=alpha,beta=beta),.(alpha,beta),
               function(df){map_plus(df$alpha,df$beta)})
df_minus<-ddply(expand.grid(alpha=alpha,beta=beta),.(alpha,beta),
                function(df){map_minus(df$alpha,df$beta)})
require(ggplot2)
qplot(x=beta, y=beta_star, data=df_plus, color=as.factor(alpha), geom="line",
      main="Duality mapping for positive x")
qplot(x=beta, y=beta_star, data=df_minus, color=as.factor(alpha), geom="line",
      main="Duality mapping for negative x")
x<-c(-2^(32:-32),2^(-32:32))
df_in<-expand.grid(alpha=alpha, beta=beta, x=x)
require(plyr)
df<-ddply(df_in, .(alpha, beta), function(df) {duality_check(df$x, df$alpha[1], df$beta[1])})
fivenum(df$eps_diff)
