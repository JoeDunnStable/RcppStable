require(stablecpp)
require(reshape2)

compare_dstable.quick<-function(xs,alphas,betas) {
  write("Comparing stablecpp::dstable to stablecpp::dstable.quick","")
  n<-length(xs)
  returnNA <-function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha\n","")
  for (a in alphas) {
    write(a,"")
    for (b in betas) {
      v_quick<-stablecpp::dstable.quick(xs,a,b,pm=0,zeta.tol=5e-5)
      v_exact<-stablecpp::dstable(xs,a,b,pm=0,zeta.tol=5e-5)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      x=xs,v_quick=v_quick,v_exact=v_exact))
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

xs<-seq(from=-1000,to=1000)/10
alphas<-c(.1,.5,1,1.5,1.99,2)
betas<-c(-1,-.5,0,.5,1)
df_d<-compare_dstable.quick(xs,alphas,betas)

