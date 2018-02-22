require(RcppStable)
require(reshape2)

compare_dstable.quick<-function(alphas,betas) {
  write("Comparing RcppStable::dstable to RcppStable::dstable.quick","")
  n<-2000
  returnNA <-function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha\n","")
  for (a in alphas) {
    write(a,"")
    for (b in betas) {
      xs<-RcppStable::qstable((.5+0:1999)/2000,a,b,pm=0)
      v_quick<-RcppStable::dstable.quick(xs,a,b,pm=0,log=T)
      v_exact<-RcppStable::dstable(xs,a,b,pm=0,log=T)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      x=xs,v_quick=v_quick,v_exact=v_exact))
      cat(sprintf("alpha = %6g, beta = %6g, loglik_quick = %20g, loglik_exact = %20g\n",
                   a, b, sum(v_quick), sum(v_exact)))
    }
  }
  df_out$diff<-with(df_out,abs(v_quick-v_exact))

  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write(noquote(paste("Number of missing v_quick =",
              format(sum(is.na(df_out$v_quick)),sep=" "))),"")
  write("The worst 100 varinaces","")
  print(head(df_out,100))
  df_out
}

alphas<-c(.1,.5,1,1.5,1.99,2)
betas<-c(-1,-.5,0,.5,1)
df_d<-compare_dstable.quick(alphas,betas)

