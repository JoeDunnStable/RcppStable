## Graphs the results of a run of stable_fit
graph_stable_fit<-function(l_in,subhead=""){
  require(ggplot2)
  parameters<-l_in$parameters
  fit_mle<-l_in$fit_mle
  label_names<-c("alpha","beta","gamma","delta","two_ll_n","pm","n")
  df_label<-data.frame(x=rep(.2,7),
                       y=c(.05+(9:3)/10),
                       value=as.double(parameters[1,label_names]),
                       confint=rep(NA,7))

  label_var<-c("alpha","beta","gamma","delta","2*loglik/n","pm","n")
  df_label$label<-c(sprintf("%11s =%8.4f",label_var[1:5],df_label$value[1:5]),
                      sprintf("%11s =%8.0f",label_var[6:7],df_label$value[6:7]))

  df_gph<-cbind(stablecpp:::stable_table,
                ghp_value=with(stablecpp:::stable_table,ifelse(type=="q_skew", value,log(value))))

  type=ifelse(nrow(parameters)==2,"mle","q")
  if (type=="mle")
    title<-paste("Fit Using McCulloch's Method for Initial Fit",
                 "and then Maximum Likelihoods",subhead, sep="\n")
  else
    title<-paste("Fit Using McCulloch's Method",subhead,sep="\n")

  gph<-ggplot(data=df_gph)+
        labs(title=title)+
        geom_contour(aes(x=alpha,y=beta,z=ghp_value,color=type))

  if (type=="mle") {
    df_label2<-data.frame(x=rep(.2,7),
                          y=-(.05+(3:9)/10),
                          value=as.double(parameters[2,label_names]),
                          confint=c(sqrt(diag(stats4::vcov(fit_mle)))*qnorm(.975),NA,NA,NA))
    df_label2$label<-c(sprintf("%11s =%8.4f +/-%8.4f",
                              label_var[1:4],df_label2$value[1:4],df_label2$confint[1:4]),
                              sprintf("%11s =%8.4f",label_var[5],df_label2$value[5]),
                              sprintf("%11s =%8.0f",label_var[6:7],df_label2$value[6:7]))

    df_label<-rbind(df_label,df_label2)
    shape <- stats4::vcov(fit_mle)[c("alpha","beta"),c("alpha","beta")]
    center <- stats4::coef(fit_mle)[c("alpha","beta")]
    radius <- sqrt(2 * qf(.95, 2, parameters[2,"n"]-2))
    angles <- (0:100) * 2 * pi/100
    unit.circle <- cbind(cos(angles), sin(angles))
##    conf_ellipse <- as.data.frame(t(center + radius * t(unit.circle %*% chol(shape))))
  }
  gph<-gph+
##    geom_polygon(mapping=aes(x=alpha,y=beta),data=conf_ellipse,fill="grey")+
    geom_point(mapping=aes(x=alpha,y=beta,shape=method),data=parameters,
               color=I("black"))+
    geom_text(data=df_label,mapping=aes(x=x,y=y,label=label),
              color=I("blue"),hjust=0,family="mono")
  show(gph)

}

require(stablecpp)

alphas=c(.1,.5,1,1.5,2)
betas<-c(-1,-.5,0,.5,1)
for (alpha in alphas){
  for (beta in betas){
    print(c(alpha,beta))
    n=10000
    set.seed(100)
    xtst<-rstable(n,alpha,beta)
    q<-quantile(xtst,p=c(.05,.25,.5,.75,.95))
    q_kurt<-(q[5]-q[1])/(q[4]-q[2])
    q_skew<-(q[5]+q[1]-2*q[3])/(q[5]-q[1])
    q_scale<-q[4]-q[2]
    q_location<-q[3]
    convergence=NA
    input_parameters<-data.frame(alpha=alpha,beta=beta,gamma=1,delta=0,pm=0,
                                 two_ll_n=2*sum(dstable.quick(xtst,alpha,beta,log=T))/n,
                                 n=n,method="input",q_kurt,q_skew,q_scale,q_location,
                                 convergence)
    row.names(input_parameters)<-NULL
    sf_out<-stable_fit(xtst,type="q_mle")
    print(rbind(input_parameters,sf_out$parameters))
  }
}
