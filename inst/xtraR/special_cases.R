dstable_holtsmark<-function(x){
  # dstable for alpha=3/2 and beta = 0
  x[4*x^6/3^6>1]<-NA
  (1/pi)*gamma(5/3)*hypergeo::genhypergeo(U=c(5/12,11/12),L=c(1/3,1/2,5/6),z=-2^2*x^6/3^6)-
    x^2/(3*pi)*hypergeo::genhypergeo(U=c(3/4,1,5/4),L=c(2/3,5/6,7/6,4/3),z=-2^2*x^6/3^6)+
    7*x^4/(3^4*pi)*gamma(4/3)*hypergeo::genhypergeo(U=c(13/12,19/12),L=c(7/6,3/2,5/3),z=-2^2*x^6/3^6)
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
  ret<-genhypergeo(U=1,L=c(.5*(mu-nu+3),.5*(mu+nu+3)),x^2/4)*x^(mu+1)/((mu+1)^2-nu^2)+
       2^(mu+nu-1)*gamma(nu)*gamma(.5*(mu+nu+1))/(gamma(.5*(-mu+nu+1)))*genhypergeo(U=NULL,L=1-nu,-x^2/4)*x^-nu+
       2^(mu-nu-1)*gamma(.5*(mu-nu+1))*gamma(-nu)/(gamma(.5*(-mu-nu+1)))*genhypergeo(U=NULL,L=1+nu,-x^2/4)*x^nu

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
