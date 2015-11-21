require(ggplot2)

show_g<-function(x,alpha,beta){
  cat(sprintf("show_g: alpha = %20.15g\n", alpha))
  cat(sprintf("show_g: beta = %20.15g\n", beta))
  if (alpha!= 1){
    tan.pa2 = tan(pi*alpha/2)
    zeta = -beta*tan.pa2
    cat(sprintf("show_g: zeta         = %20.15g\n",zeta))
    if (alpha<1){
      if (beta==1)
        theta0=pi/2
      else if (beta==-1)
        theta0=-pi/2
      else
        theta0<-min(max(-pi/2,atan(-zeta)/alpha),pi/2)
    } else {
      if (beta==1)
        theta0=pi/2-pi/alpha
      else if (beta==-1)
        theta0=pi/alpha-pi/2
      else
        theta0<-min(max(-pi/2,atan(-zeta)/alpha),pi/2)
    }
    if (x<zeta) theta0<--theta0
    cat(sprintf("show_g: theta0       = %20.15g\n",theta0))
    x_m_zeta<-abs(x-zeta)
    cat(sprintf("show_g: x_m_zeta     = %20.15g\n",x_m_zeta))
    at0<-alpha*theta0
    cat(sprintf("show_g: alpha*theta0 = %20.15g\n",at0))
    cat0<-cos(at0)
    c2 <- (alpha/(pi * abs(alpha-1) * x_m_zeta))
    cat(sprintf("show_g: c2 =           %20.15g\n",c2))
    add_l<-(theta0-pi/2)
    add_r<-(pi-alpha*(theta0+pi/2))
    fun_g<-function(th) {
      costh<-cos(th)
      att<-alpha*th+at0
      x_sin<-x_m_zeta/sin(att)                #  ln(x_m_zeta)-ln(alpha)                                          -ln(th)
      pow1<-(x_sin)^alpha                     #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
      pow2<-(cat0*costh*pow1)^(1/(alpha-1))   #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
      catt_m_t<-cos(att-th)
      if (abs(beta)==1){
        g0<-((cat0*(x_m_zeta/alpha)^alpha))^(1/(alpha-1))*abs(1-alpha)
        if (alpha<1)(
          if (th==-theta0)
            r<-g0
          else
            r<-pow2*catt_m_t
        ) else if (theta0>0 && th==pi/2)
          r<-g0
        else
          r<-pow2*catt_m_t
      } else
        r<-pow2*catt_m_t
      r-1
    }
    fun_g_l<-function(th_l) {
      costh<-sin(th_l-add_l)
      att<-alpha*th_l
      x_sin<-x_m_zeta/sin(att)                #  ln(x_m_zeta)-ln(alpha)                                          -ln(th)
      pow1<-(x_sin)^alpha                     #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
      pow2<-(cat0*costh*pow1)^(1/(alpha-1))   #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
      catt_m_t<--sin((alpha-1)*th_l+add_l)
      if (abs(beta)==1){
        g0<-((cat0*(x_m_zeta/alpha)^alpha))^(1/(alpha-1))*abs(1-alpha)
        if (alpha<1)(
          if (th_l==0)
            r<-g0
          else
            r<-pow2*catt_m_t
        ) else if (theta0>0 && th_l==pi/2+theta0)
          r<-g0
        else
          r<-pow2*catt_m_t
      } else
        r<-pow2*catt_m_t
      r-1
    }
    fun_g_r<-function(th_r) {
      costh<-sin(th_r)
      att<-alpha*th_r+add_r
      x_sin<-x_m_zeta/sin(att)                #  ln(x_m_zeta)-ln(alpha)                                          -ln(th)
      pow1<-(x_sin)^alpha                     #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
      pow2<-(cat0*costh*pow1)^(1/(alpha-1))   #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
      catt_m_t<-sin((alpha-1)*th_r+add_r)
      if (abs(beta)==1){
        g0<-((cat0*(x_m_zeta/alpha)^alpha))^(1/(alpha-1))*abs(1-alpha)
        if (alpha<1)(
          if (th_r==pi/2+theta0)
            r<-g0
          else
            r<-pow2*catt_m_t
        ) else if (theta0>0 && th_r==0)
          r<-g0
        else
          r<-pow2*catt_m_t
      } else
        r<-pow2*catt_m_t
      r-1
    }
    f_mth0<-fun_g(-theta0)
    f_pi2<-fun_g(pi/2)
    if (f_mth0*f_pi2<0){
      thr<-uniroot(fun_g,c(-theta0,pi/2),tol=64*.Machine$double.eps)
      th2<-thr$root
      cat(sprintf("show_g: th2 =           %20.15g, g(th2) = %20.15g\n",th2,thr$f.root+1))
    }else{
      if (is.finite(f_mth0)){
        th2<--theta0
        cat(sprintf("show_g: th2 =           %20.15g, g(th2) = %20.15g\n",th2,f_mth0+1))
      } else {
        th2<-pi/2
        cat(sprintf("show_g: th2 =           %20.15g, g(th2) = %20.15g\n",th2,f_pi2+1))
      }
    }
#    f_l_lo<-fun_g_l(0)
#    f_l_hi<-fun_g_l(pi/2+theta0)
#    if (f_l_lo*f_l_hi<0){
#      thr<-uniroot(fun_g_l,c(0,pi/2+theta0),tol=64*.Machine$double.eps)
#      th2_l<-thr$root
#      cat(sprintf("show_g: th2_l-theta0 =    %20.15g, g(th2_l) = %20.15g\n",th2_l-theta0,thr$f.root+1))
#    }else{
#      if (is.finite(f_l_lo)){
#        th2_l<-0
#        cat(sprintf("show_g: th2_l-theta0 =   %20.15g, g(th2_l) = %20.15g\n",th2_l-theta0,f_l_lo+1))
#      } else {
#        th2_l<-pi/2+theta0
#        cat(sprintf("show_g: th2_l-theta0 =   %20.15g, g(th2_l) = %20.15g\n",th2_l-theta0,f_l_hi+1))
#      }
#    }
#    f_r_hi<-fun_g_r(pi/2+theta0)
#    f_r_lo<-fun_g_r(0)
#    if (f_r_hi*f_r_lo<0){
#      thr<-uniroot(fun_g_r,c(0,pi/2+theta0),tol=64*.Machine$double.eps)
#      th2_r<-thr$root
#      cat(sprintf("show_g: pi/2-th2_r =   %20.15g, g(th2) = %20.15g\n",pi/2-th2_r,thr$f.root+1))
#    }else
#      {
#      if (is.finite(f_r_lo)){
#        th2_r<-0
#        cat(sprintf("show_g: pi/2-th2_r =   %20.15g, g(th2_r) = %20.15g\n",pi/2-th2_r,f_r_lo+1))
#      } else {
#        th2_r<-pi/2+theta0
#        cat(sprintf("show_g: pi/2-th2_r =   %20.15g, g(th2_r) = %20.15g\n",pi/2-th2_r,f_r_hi+1))
#      }
#    }
    d_th<-(theta0+pi/2)/1000
    th<-c(-theta0,-theta0+(10^-(10:1))*d_th,seq(-theta0+d_th,pi/2-d_th,length.out=999),pi/2-(10^-(1:10))*d_th,pi/2)
    th<-sort(c(th,th2))
    th_l<-c(0,(10^-(10:1))*d_th,d_th*(1:999),(pi/2+theta0)-(10^-(1:10))*d_th,pi/2+theta0)
    th_l<-sort(c(th_l,th2+theta0))
    th_r<-c(pi/2+theta0,pi/2+theta0-(10^-(10:1)*d_th),d_th*(999:1),(10^-(1:10))*d_th,0)
    th_r<-sort(c(th_r,pi/2-th2),decreasing=T)
    costh<-cos(th)
    cat(sprintf("show_g: theta0-pi/2 = %20.15g\n",add_l))
    costh_l<-sin(th_l-add_l)
    costh_r<-sin(th_r)
    att=alpha*th+at0
    att_l<-alpha*th_l
    cat(sprintf("show_g: pi-alpha*(theta0+pi/2) = %20.15g\n",add_r))
    att_r<-alpha*th_r+add_r     #alpha*th_r+(pi-alpha*(pi/2+theta0))
    x_sin<-x_m_zeta/sin(att)                   #  ln(x_m_zeta)-ln(alpha)                                          -ln(th)
    x_sin_l<-x_m_zeta/sin(att_l)
    x_sin_r<-x_m_zeta/sin(att_r)
    pow1<-(x_sin)^alpha            #  alpha*(ln(x_m_zet)-ln(alpha))                             -alpha*ln(th)
    pow1_l<-(x_sin_l)^alpha
    pow1_r<-(x_sin_r)^alpha
    pow2<-(cat0*costh*pow1)^(1/(alpha-1))    #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1)               -ln(th)
    pow2_l<-(cat0*costh_l*pow1_l)^(1/(alpha-1))
    pow2_r<-(cat0*costh_r*pow1_r)^(1/(alpha-1))
    catt_m_t<-cos(att-th)
    catt_m_t_l<--sin((alpha-1)*th_l+add_l)
    catt_m_t_r<-sin((alpha-1)*th_r+add_r)
    g<-pow2*catt_m_t                        #  (ln(cat0)+alpha*(ln(x_m_zet)-ln(alpha))/(alpha-1))+ln(1-alpha) 0*ln(th)
    g_l<-pow2_l*catt_m_t_l
    g_r<-pow2_r*catt_m_t_r
    if (abs(beta)==1){
      g0<-((cat0*(x_m_zeta/alpha)^alpha))^(1/(alpha-1))*abs(1-alpha)
      cat(sprintf("show_g: g0 = %20.15g",g0))
      if (alpha<1){
        sel<-th==-theta0
        g[sel]<-g_l[sel]<-g_r[sel]<-g0
      } else if (theta0>0){
        sel<-th==pi/2
        g[sel]<-g_l[sel]<-g_r[sel]<-g0
      }
    }
    g_tst<-1+vapply(th,fun_g,0.)
    g_l_tst<-1+vapply(th_l,fun_g_l,0.)
    g_r_tst<-1+vapply(th_r,fun_g_r,0.)
    if (x_m_zeta<=.01*min(1,abs(zeta)))
      g_cpp<-vapply(th_l,stablecpp:::test_g,0.,x=x,alpha=alpha,beta=beta)
    else
      g_cpp<-vapply(th_r,stablecpp:::test_g,0.,x=x,alpha=alpha,beta=beta)
    g_exp_m_g<-g*exp(-g)
    g_exp_m_g[g==Inf]<-0
    ln_cat0<-log(cat0)
    ln_costh<-log(cos(th))
    ln_x_sin<-log(x_m_zeta)-log(sin(att))
    ln_pow1<-alpha*ln_x_sin
    ln_pow2<-(ln_cat0+ln_costh+ln_pow1)/(alpha-1)
    ln_catt_m_t<-log(cos(att-th))
    ln_g<-ln_pow2+ln_catt_m_t
    df_out<-data.frame(
               th=th,
               th_l=th_l,
               th_r=th_r,
               cat0=rep(cat0,length(th)),
               costh=costh,
               costh_l=costh_l,
               costh_r=costh_r,
               x_sin=x_m_zeta/sin(att),
               x_sin_l=x_sin_l,
               x_sin_r=x_sin_r,
               pow1=pow1,
               pow1_l=pow1_l,
               pow1_r=pow1_r,
               pow2=pow2,
               pow2_l=pow2_l,
               pow2_r=pow2_r,
               catt_m_t=catt_m_t,
               catt_m_t_l=catt_m_t_l,
               catt_m_t_r=catt_m_t_r,
               g=g,
               g_tst=g_tst,
               g_l=g_l,
               g_l_tst=g_l_tst,
               g_r=g_r,
               g_r_tst=g_r_tst,
               g_cpp=g_cpp,
               g_exp_m_g=g_exp_m_g,
               ln_dth=log(th-th[1]),
               ln_cat0=rep(ln_cat0,length(th)),
               ln_costh=ln_costh,
               ln_x_sin=ln_x_sin,
               ln_pow1=ln_pow1,
               ln_pow2=ln_pow2,
               ln_catt_m_t=ln_catt_m_t,
               ln_g=ln_g)

  } else if (alpha==1){
    x_in<-x
    beta_in<-beta
    if (x<0) {
      x<--x
      beta<--beta
    }
    i2b <- 1/(2*beta)
    cat(sprintf("show_g: 1/(2*beta) = %20.15g\n",i2b))
    p2b <- pi*i2b ## = pi/(2 beta)
    cat(sprintf("show_g: pi/(2*beta) = %20.15g\n",p2b))
    ea <- -p2b*x
    cat(sprintf("show_g: ea = %20.15g\n",ea))
    cc<-pi/2*abs(i2b)
    cat(sprintf("show_g: cc = %20.15g\n",cc))
    u0 <- -((beta>=0)-(beta<=0))
    fun_g<-function(u) {

      h<-p2b+u*pi/2
      h2b<-h/p2b
      if(u==-1){
        tanth<--Inf
        if (beta == 1)
           h_tanth<--1
        else
           h_tanth<--Inf
      } else if (u==1) {
        tanth<-Inf
        if (beta == -1)
          h_tanth<--1
        else
          h_tanth<-Inf
      } else {
        tanth<-tanpi(u/2)
        h_tanth<-h*tanth
      }
      exp_ea_p_h_tan_th<-exp(ea+h_tanth)
      costh<-cospi(u/2)
      if (u==-1) {
        if (beta>0) {
          if (beta==1)
            r<-exp_ea_p_h_tan_th*2/pi
          else
            r<-0
        } else {
          r<-Inf
        }
      } else if (u==1){
        if (beta<0){
          if (beta==-1)
            r<-exp_ea_p_h_tan_th*2/pi
          else
            r<-0
        } else
          r <- Inf
      } else if (exp_ea_p_h_tan_th ==0)
        r<-0
      else
        r<-h2b*exp(ea+h_tanth)/costh
      r-1
    }
    fun_g_r<-function(u_r) {

      h<-p2b+pi/2-u_r*pi/2
      h2b<-h/p2b
      if(u_r==2){
        tanth<--Inf
        if (beta == 1)
          h_tanth<--1
        else
          h_tanth<--Inf
      } else if (u_r==0) {
        tanth<-Inf
        if (beta == -1)
          h_tanth<--1
        else
          h_tanth<-Inf
      } else {
        tanth<-cospi(u_r/2)/sinpi(u_r/2)
        h_tanth<-h*tanth
      }
      exp_ea_p_h_tan_th<-exp(ea+h_tanth)
      costh<-sinpi(u_r/2)
      if (u_r==2) {
        if (beta>0) {
          if (beta==1)
            r<-exp_ea_p_h_tan_th*2/pi
          else
            r<-0
        } else {
          r<-Inf
        }
      } else if (u_r==0){
        if (beta<0){
          if (beta==-1)
            r<-exp_ea_p_h_tan_th*2/pi
          else
            r<-0
        } else
          r <- Inf
      } else if (exp_ea_p_h_tan_th ==0)
        r<-0
      else
        r<-h2b*exp(ea+h_tanth)/costh
      r-1
    }
    f_m1<-fun_g(-1)
    f_p1<-fun_g(+1)
    if (f_m1*f_p1<0){
      ur<-uniroot(fun_g,c(-1,1),tol=64*.Machine$double.eps)
      u2<-ur$root
      cat(sprintf("show_g: th2 = %20.15g, g(u2) = %20.15g\n",u2*pi/2,ur$f.root+1))
    }else{
      if (is.finite(f_m1)){
        u2<--1
        cat(sprintf("show_g: th2 = %20.15g, g(u2) = %20.15g\n",u2*pi/2,f_m1+1))
      } else {
        u2<-1
        cat(sprintf("show_g: th2 = %20.15g, g(u2) = %20.15g\n",u2*pi/2,f_p1+1))
      }
    }
    f_r_lo<-fun_g_r(0)
    f_r_hi<-fun_g_r(2)
    if (f_r_lo*f_r_hi<0){
      ur<-uniroot(fun_g_r,c(0,2),tol=64*.Machine$double.eps)
      u2_r<-ur$root
      cat(sprintf("show_g: pi/2*(1-u2_r) = %20.15g, g(u2) = %20.15g\n",pi/2*(1-u2_r),ur$f.root+1))
    }else{
      if (is.finite(f_r_lo)){
        u2_r<-0
        cat(sprintf("show_g: pi/2*(1-u2_r) = %20.15g, g(u2) = %20.15g\n",pi/2*(1-u2_r),f_r_lo+1))
      } else {
        u2_r<-2
        cat(sprintf("show_g: pi/2*(1-u2_r) = %20.15g, g(u2) = %20.15g\n",pi/2*(1-u2_4),f_r_hi+1))
      }
    }
    d_u<-2/1000
    u<-c(-1,-1+(10^-(10:1))*d_u,seq(-1+d_u,1-d_u,length.out=999),1-(10^-(1:10))*d_u,1)
    u<-sort(c(u,u2))
    u_l<-c(0,(10^-(10:1))*d_u,seq(d_u,2-d_u,length.out=999),2-(10^-(1:10))*d_u,2)
    u_r<-rev(u_l)
    u_l<-sort(c(u_l,u2+1))
    u_r<-sort(c(u_r,1-u2),decreasing=T)
    g_cpp<-vapply(u_r,stablecpp:::test_g,0.,x=x_in,alpha=alpha,beta=beta_in)
    n<-length(u)
    th = u*pi/2
    th_l<-u_l*pi/2
    th_r<-u_r*pi/2
    h = p2b + th
    h_l<-th_l+(p2b-pi/2)
    h_r<--th_r+(p2b+pi/2)
    h_2b<-(h/p2b)
    h_2b_l<-(h_l/p2b)
    h_2b_r<-(h_r/p2b)
    tanth<-tanth_l<-tanth_r<-rep(0,n)
    sel<-abs(u)!=1
    tanth[sel]<-tanpi(u[sel]/2)
    tanth_l[sel]<--cospi(u_l[sel]/2)/sinpi(u_l[sel]/2)
    tanth_r[sel]<-cospi(u_r[sel]/2)/sinpi(u_r[sel]/2)
    tanth[1]<-tanth_l[1]<-tanth_r[1]<--Inf
    tanth[n]<-tanth_l[n]<-tanth_r[n]<-Inf
    h_tan_th<-h*tanth
    h_tan_th_l<-h_l*tanth_l
    h_tan_th_r<-h_r*tanth_r
    if (abs(beta)==1){
      if (beta==1) {
        h_tan_th[1]<-h_tan_th_l[1]<-h_tan_th_r[1]<--1
      } else {
        h_tan_th[n]<-h_tan_th_l[n]<-h_tan_th_r[n]<--1
      }
    }
    costh<-cospi(u/2)
    costh_l<-sinpi(u_l/2)
    costh_r<-sinpi(u_r/2)
    exp_ea_p_h_tan_th<-exp(ea + h_tan_th)
    exp_ea_p_h_tan_th_l<-exp(ea + h_tan_th_l)
    exp_ea_p_h_tan_th_r<-exp(ea + h_tan_th_r)
    sel<-exp_ea_p_h_tan_th != 0
    g<-rep(0,n)
    g_l<-rep(0,n)
    g_r<-rep(0,n)
    g[sel] <- h_2b[sel] * exp_ea_p_h_tan_th[sel] / costh[sel]
    g_l[sel] <- h_2b_l[sel] * exp_ea_p_h_tan_th_l[sel] / costh_l[sel]
    g_r[sel] <- h_2b_r[sel] * exp_ea_p_h_tan_th_r[sel] / costh_r[sel]
    if (abs(beta)==1){
      if (beta==1) {
        g[1]<-exp_ea_p_h_tan_th[1]*2/pi
        g_l[1]<-exp_ea_p_h_tan_th_l[1]*2/pi
        g_r[1]<-exp_ea_p_h_tan_th_r[1]*2/pi
      } else {
        g[n]<-exp_ea_p_h_tan_th[n]*2/pi
        g_l[n]<-exp_ea_p_h_tan_th_l[n]*2/pi
        g_r[n]<-exp_ea_p_h_tan_th_r[n]*2/pi
      }
    }
    g_tst<-1+vapply(u,fun_g,0.)
    g_r_tst<-1+vapply(u_r,fun_g_r,0.)
    g_exp_m_g<-g*exp(-g)
    g_exp_m_g_l<-g_l*exp(-g_l)
    g_exp_m_g_r<-g_r*exp(-g_r)
    g_exp_m_g[g==Inf]<-0
    g_exp_m_g_l[g_l==Inf]<-0
    g_exp_m_g_r[g_r==Inf]<-0
    df_out<-data.frame(
              u=u,
              u_l=u_l,
              u_r=u_r,
              th=th,
              th_l=th_l,
              th_r=th_r,
              h=h,
              h_l=h_l,
              h_r=h_r,
              h_2b=h_2b,
              h_2b_l=h_2b_l,
              h_2b_r=h_2b_r,
              h_tan_th=h_tan_th,
              h_tan_th_l=h_tan_th_l,
              h_tan_th_r=h_tan_th_r,
              exp_ea_p_h_tan_th=exp_ea_p_h_tan_th,
              exp_ea_p_h_tan_th_l=exp_ea_p_h_tan_th_l,
              exp_ea_p_h_tan_th_r=exp_ea_p_h_tan_th_r,
              costh=costh,
              costh_l=costh_l,
              costh_r=costh_r,
              g=g,
              g_tst=g_tst,
              g_l=g_l,
              g_r=g_r,
              g_r_tst=g_r_tst,
              g_cpp=g_cpp,
              g_exp_m_g=g_exp_m_g,
              g_exp_m_g_l=g_exp_m_g_l,
              g_exp_m_g_r=g_exp_m_g_r)
  }
  show(qplot(x=th,y=g_exp_m_g,data=df_out,geom="line",xlim=c(-pi/2,pi/2)))
  df_out
}
