require(stabledist)
require(RcppStable)
require(reshape2)

check_stableMode<-function(){
  write("Test of stableMode for beta close to 1\n","")
  alpha <- seq(.1, 2, by = 0.1)
  beta <- c(.99,.99999,.99999999,.99999999999)
  ans<-data.frame()
  for (a in alpha)
    for (b in beta){
    ans<-rbind(ans,data.frame(alpha=a,beta=b,mode=stableMode(a,b)))
  }
  print(acast(ans,alpha ~ beta,value.var="mode"))
  ans
}
##
##	alpha	       0.99	  0.99999    0.99999999 0.99999999999
##	0.0    0.000000e+00  0.000000e+00  0.000000e+00	 0.000000e+00
##	0.2   -3.214142e-01 -3.246759e-01 -3.246787e-01 -3.246788e-01
##	0.4   -6.105318e-01 -6.158562e-01 -6.158616e-01 -6.158616e-01
##	0.6   -6.550106e-01 -6.594746e-01 -6.594790e-01 -6.594790e-01
##	0.8   -5.558811e-01 -5.590032e-01 -5.590063e-01 -5.590063e-01
##	1.0   -4.271033e-01 -4.293078e-01 -4.293099e-01 -4.293099e-01
##	1.2   -3.074015e-01 -3.090820e-01 -3.090804e-01 -3.090804e-01
##	1.4   -2.050956e-01 -2.063979e-01 -2.063951e-01 -2.063951e-01
##	1.6   -1.199623e-01 -1.208875e-01 -1.208853e-01 -1.208853e-01
##	1.8   -5.098617e-02 -5.145758e-02 -5.145639e-02 -5.145639e-02
##	2.0   -7.487432e-05 -7.487432e-05 -7.487432e-05 -7.487432e-05


compare_dstable<-function(xs,alphas,betas) {
  write("Comparing RcppStable::dstable to stabledist::dstable","")
  n<-length(xs)
  returnNA <-function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha\n","")
  for (a in alphas) {
    write(a,"")
    for (b in betas) {
      v_cpp<-tryCatch(RcppStable::dstable(xs,a,b,pm=0),error=returnNA)
      v_r<-tryCatch(stabledist::dstable(xs,a,b,pm=0),error=returnNA)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      x=xs,v_cpp=v_cpp,v_r=v_r))
    }
  }
  df_out$diff<-with(df_out,abs(v_cpp-v_r))
  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write(noquote(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp)),sep=" "))),"")
  write("The worst 100 varinaces","")
  print(head(df_out,100))
  df_out
}

compare_pstable<-function(xs,alphas,betas) {
  write("Comparing RcppStable::pstable to stabledist::pstable\n","")
  n<-length(xs)
  returnNA <- function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha\n","")
    for (a in alphas) {
    write(a,"")
    for (b in betas) {
      v_cpp<-tryCatch(RcppStable::pstable(xs,a,b,pm=0),error=returnNA)
      v_r<-tryCatch(stabledist::pstable(xs,a,b,pm=0),error=returnNA)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      x=xs,v_cpp=v_cpp,v_r=v_r))
    }
  }
  write(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp))),
              sep=" "),"")
  df_out$diff<-with(df_out,abs(v_cpp-v_r))
  df_out<-df_out[order(df_out$diff, decreasing=T),]
  write("The worst 100 varinaces\n","")
  print(head(df_out,100))
  df_out
}

compare_qstable<-function(ps,alphas,betas) {
  write("Comparing RcppStable::qstable to stabledist::qstable\n","")
  n<-length(ps)
  returnNA <- function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha","")
  for (a in alphas) {
    print(a)
    for (b in betas) {
      v_cpp<-tryCatch(RcppStable::qstable(ps,a,b,pm=0),error=returnNA)
      v_r<-tryCatch(stabledist::qstable(ps,a,b,pm=0),error=returnNA)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      p=ps,v_cpp=v_cpp,v_r=v_r))
    }
  }
  write(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp))),
                     sep=" "),"")
  df_out$diff <- with(df_out,abs(v_cpp-v_r))
  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write("The worst 100 varinaces\n","")
  print(head(df_out,100))
  df_out
}

smode_out<-check_stableMode()

xs<-c(-1000,-100,seq(from=-10,to=10,by=.5),100,1000)
ps<-c(.0001,.001,.01,1:9/10,.99,.999,.9999)
alphas<-c(.1,.5,1,1.5,2)
betas<-c(-1,-.5,0,.5,1)
df_d<-compare_dstable(xs,alphas,betas)
df_p<-compare_pstable(xs,alphas,betas)
df_q<-compare_qstable(ps,alphas,betas)

