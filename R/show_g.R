require(ggplot2)

show_g<-function(x,alpha,beta){
  cat(sprintf("show_g: alpha = %20.15g\n", alpha))
  cat(sprintf("show_g: beta = %20.15g\n", beta))
  if (alpha!= 1){
    tan.pa2 = tan(pi*alpha/2)
    zeta = -beta*tan.pa2
    cat(sprintf("show_g: zeta = %20.15g\n",zeta))
    theta0<-min(max(-pi/2,atan(-zeta)/alpha),pi/2)
    cat(sprintf("show_g: theta0 = %20.15g\n",theta0))
    x_m_zeta<-abs(x-zeta)
    cat(sprintf("show_g: x_m_zeta = %20.15g\n",x_m_zeta))
    at0<-alpha*theta0
    cat(sprintf("show_g: alpha*theta0 = %20.15g\n",at0))
    cat0<-cos(at0)
    th=seq(from=-theta0,to=pi/2, length.out=200)
    att=alpha*th+at0
    x_sin<-x_m_zeta/sin(att)
    pow1=(x_m_zeta/sin(att))^alpha
    pow2=(cat0*cos(th)*pow1)^(1/(alpha-1))
    data.frame(th=th,
               cat0=rep(cat0,length(th)),
               costh=cos(th),
               x_sin=x_m_zeta/sin(att),
               pow1=pow1,
               pow2=pow2,
               catt_m_t=cos(att-th),
               g=pow2*cos(att-th))
  }
}
