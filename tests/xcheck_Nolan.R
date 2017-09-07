require(stablecpp)
require(reshape2)
require(plyr)

cdir<-getwd()
if (basename(cdir)=="tests") {
  load("../../stablecpp/inst/extdata/Nolan_v3.14.02.RData")
} else {
  load ("inst/extdata/Nolan_v3.14.02.RData")
}

compare_dstable<-function() {
  write("Comparing stablecpp::dstable to to Nolan STABLE v3.14.02","")
  df_out<-ddply(df_Nolan_3.14.02,.(alpha,beta),
                 function(df) data.frame(alpha=df$alpha,beta=df$beta,x=df$x,
                                         v_cpp=dstable(df$x,df$alpha[1],df$beta[1]),
                                         v_nolan=df$dstable))
  df_out$diff<-with(df_out,abs(v_cpp-v_nolan))
  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write(noquote(paste("Number of missing v_cpp =",
              format(sum(is.na(df_out$v_cpp)),sep=" "))),"")
  write("The worst 100 varinaces","")
  print(head(df_out,100))
  df_out
}

compare_pstable<-function() {
  write("Comparing stablecpp::pstable to pstable from Nolan STABLE v3.14.02","")
  df_out<-ddply(df_Nolan_3.14.02,.(alpha,beta),
                 function(df) data.frame(alpha=df$alpha,beta=df$beta,x=df$x,
                                         v_cpp=pstable(df$x,df$alpha[1],df$beta[1]),
                                         v_nolan=df$pstable))

  df_out$diff<-with(df_out,abs(v_cpp-v_nolan))
  df_out<-df_out[order(df_out$diff,decreasing=T),]
  write(noquote(paste("Number of missing v_cpp =",
                      format(sum(is.na(df_out$v_cpp)),sep=" "))),"")
  write("The worst 100 varinaces","")
  print(head(df_out,100))
  df_out
}

df_d<-compare_dstable()
df_p<-compare_pstable()

