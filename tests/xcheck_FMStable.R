require(FMStable)
require(RcppStable)
require(reshape2)


compare_dstable<-function(xs,alphas) {
  write("Comparing RcppStable::dstable to FMStable::dEstable","")
  n<-length(xs)
  FM<-function(x,a,b,pm){
    if (b==-1) {
      b=1
      x<--x
    }
    obj<-FMStable::setParam(alpha=a,location=0,logscale=0,pm=pm)
    dEstable(x,obj)
  }
  returnNA <-function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha\n","")
  for (a in alphas) {
    write(a,"")
    for (b in c(-1,1)) {
      v_cpp<-tryCatch(RcppStable::dstable(xs,a,b,pm=0),error=returnNA)
      v_fm<-tryCatch(FM(xs,a,b,pm=0),error=returnNA)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      x=xs,v_cpp=v_cpp,v_fm=v_fm))
    }
  }
  df_out$diff<-with(df_out,abs(v_cpp-v_fm))
  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write(noquote(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp)),sep=" "))),"")
  write("The worst 100 varinaces","")
  print(head(df_out,100))
  df_out
}

compare_pstable<-function(xs,alphas) {
  write("Comparing RcppStable::pstable to FMStable::pEstable\n","")
  n<-length(xs)
  FM<-function(x,a,b,pm){
    if (b==-1) {
      x<--x
    }
    obj<-FMStable::setParam(alpha=a,location=0,logscale=0,pm=pm)
    r<-pEstable(x,obj)
    if (b==-1)
      1-r
    else
      r
  }
  returnNA <- function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha\n","")
    for (a in alphas) {
    write(a,"")
    for (b in c(-1,1)) {
      v_cpp<-tryCatch(RcppStable::pstable(xs,a,b,pm=0),error=returnNA)
      v_fm<-tryCatch(FM(xs,a,b,pm=0),error=returnNA)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      x=xs,v_cpp=v_cpp,v_fm=v_fm))
    }
  }
  write(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp))),
              sep=" "),"")
  df_out$diff<-with(df_out,abs(v_cpp-v_fm))
  df_out<-df_out[order(df_out$diff, decreasing=T),]
  write("The worst 100 varinaces\n","")
  print(head(df_out,100))
  df_out
}

compare_qstable<-function(ps,alphas) {
  write("Comparing RcppStable::qstable to FMStable::qEstable\n","")
  n<-length(ps)
  FM<-function(p,a,b,pm){
    if (b==-1)
      p=1-p
    obj<-FMStable::setParam(alpha=a,location=0,logscale=0,pm=pm)
    r<-FMStable::qEstable(p,obj)
    if (b==1)
      r
    else
      -r
  }
  returnNA <- function(e) rep(NA,n)
  df_out<-data.frame()
  write("alpha","")
  for (a in alphas) {
    print(a)
    for (b in c(-1,1)) {
      v_cpp<-tryCatch(RcppStable::qstable(ps,a,b,pm=0),error=returnNA)
      v_fm<-tryCatch(FM(ps,a,b,pm=0),error=returnNA)
      df_out<-rbind(df_out,data.frame(alpha=rep(a,n),beta=rep(b,n),
                                      p=ps,v_cpp=v_cpp,v_fm=v_fm))
    }
  }
  write(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp))),
                     sep=" "),"")
  df_out$diff <- with(df_out,abs(v_cpp-v_fm))
  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write("The worst 100 varinaces\n","")
  print(head(df_out,100))
  df_out
}

xs<-c(-1e6,-1e6,-1e4,-1000,-100,seq(from=-10,to=10,by=.5),100,1000,1e4,1e5,1e6)
ps<-c(.000001,.00001,.0001,.001,.01,1:9/10,.99,.999,.9999,.99999,.999999)
alphas<-c(.1,.5,1,1.5,2)
df_d<-compare_dstable(xs,alphas)
df_p<-compare_pstable(xs,alphas)
df_q<-compare_qstable(ps,alphas)

