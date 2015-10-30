require(ggplot2)

show_g<-function(x,alpha,beta){
  cat(sprintf("show_g: alpha = %20.15g\n", alpha))
  cat(sprintf("show_g: beta = %20.15g\n", beta))
  if (alpha!= 1){
    tan.pa2 = tan(pi*alpha/2)
    zeta = -beta*tan.pa2
    cat(sprintf("show_g: zeta         = %20.15g\n",zeta))
    theta0<-min(max(-pi/2,atan(-zeta)/alpha),pi/2)
    if (x<zeta) theta0<--theta0
    cat(sprintf("show_g: theta0       = %20.15g\n",theta0))
    x_m_zeta<-abs(x-zeta)
    cat(sprintf("show_g: x_m_zeta     = %20.15g\n",x_m_zeta))
    at0<-alpha*theta0
    cat(sprintf("show_g: alpha*theta0 = %20.15g\n",at0))
    cat0<-cos(at0)
    c2 <- (alpha/(pi * abs(alpha-1) * x_m_zeta))
    cat(sprintf("show_g: c2 =           %20.15g\n",c2))
    d_th<-(theta0+pi/2)/1000
    th<-c(-theta0,-theta0+(10^-(10:1))*d_th,seq(-theta0+d_th,pi/2-d_th,length.out=999),pi/2-(10^-(1:10))*d_th,pi/2)
    th_l<-c(0,(10^-(10:1))*d_th,d_th*(1:999),(pi/2+theta0)-(10^-(1:10))*d_th,pi/2+theta0)
    th_r<-c(pi/2+theta0,pi/2+theta0-(10^-(10:1)*d_th),d_th*(999:1),(10^-(1:10))*d_th,0)
    att=alpha*th+at0
    att_l<-alpha*th_l
    att_r<-alpha*(pi/2+theta0-th_r)
    x_sin<-x_m_zeta/sin(att)                   #  ln(x_m_zeta)-ln(alpha)                                          -ln(th)
    x_sin_l<-x_m_zeta/sin(att_l)
    x_sin_r<-x_m_zeta/sin(att_r)
    pow1<-(x_sin)^alpha            #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
    pow1_l<-(x_sin_l)^alpha            #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
    pow1_r<-(x_sin_r)^alpha            #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
    pow2<-(cat0*cos(th)*pow1)^(1/(alpha-1))    #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
    pow2_l<-(cat0*cos(th_l-theta0)*pow1_l)^(1/(alpha-1))    #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
    pow2_r<-(cat0*cos(pi/2-th_r)*pow1_r)^(1/(alpha-1))    #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
    g<-pow2*cos(att-th)                        #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1))+ln(1-alpha) 0*ln(th)
    g_l<-pow2_l*cos((alpha-1)*th_l+theta0)                        #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1))+ln(1-alpha) 0*ln(th)
    g_r<-pow2_r*cos((alpha-1)*pi/2 +alpha*theta0-(alpha-1)*th_r)                        #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1))+ln(1-alpha) 0*ln(th)
    g0<-((cat0*(x_m_zeta/alpha)^alpha))^(1/(alpha-1))*(1-alpha)
    if (beta==1)
      cat(sprintf("show_g: g0 = %20.15g",g0))
    g_exp_m_g<-g*exp(-g)
    g_exp_m_g[g==Inf]<-0
    ln_cat0<-log(cat0)
    ln_costh<-log(cos(th))
    ln_x_sin<-log(x_m_zeta)-log(sin(att))
    ln_pow1<-alpha*ln_x_sin
    ln_pow2<-(ln_cat0+ln_costh+ln_pow1)/(alpha-1)
    ln_catt_m_t<-log(cos(att-th))
    ln_g<-ln_pow2+ln_catt_m_t
      df<<-data.frame(th=th,
               th_l=th_l,
               th_r=th_r,
               cat0=rep(cat0,length(th)),
               costh=cos(th),
               x_sin=x_m_zeta/sin(att),
               x_sin_l=x_sin_l,
               x_sin_r=x_sin_r,
               pow1=pow1,
               pow1_l=pow1_l,
               pow1_r=pow1_r,
               pow2=pow2,
               pow2_l=pow2_l,
               pow2_r=pow2_r,
               catt_m_t=cos(att-th),
               g=g,
               g_l=g_l,
               g_r=g_r,
               g_exp_m_g=g_exp_m_g,
               ln_dth=log(th-th[1]),
               ln_cat0=rep(ln_cat0,length(th)),
               ln_costh=ln_costh,
               ln_x_sin=ln_x_sin,
               ln_pow1=ln_pow1,
               ln_pow2=ln_pow2,
               ln_catt_m_t=ln_catt_m_t,
               ln_g=ln_g)
    show(qplot(x=th,y=g_exp_m_g,data=df,geom="line"))
#       geom_line(aes(x=th,y=exp(ln_g)*exp(-exp(ln_g))),color="red"))

  }
}
